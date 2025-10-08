import os, sys, glob
from tqdm import tqdm
import numpy as np
import urllib
import requests
import aiohttp
import warnings
import contextlib

import pyvo

import multiprocess as mp

from astropy.table import Column
from astropy.io import fits, ascii
from astropy.table import Table, vstack, hstack, QTable
from astropy import cosmology
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import constants as const
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.visualization import quantity_support
from astropy.nddata import Cutout2D

from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry

from astroquery.ipac.irsa import Irsa

from scipy import ndimage
from scipy.interpolate import interp2d, RegularGridInterpolator, interp1d


import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as path_effects

import pandas as pd

# 物理定数
c = 2.9979e10  # [cm/s]
pc = 3.085677857083127e18  # [cm]
Jy = 1.0e-23   # [erg/s/cm^2/Hz]

## Plotting stuff
mpl.rcParams['font.size'] = 14
mpl.rcParams['axes.labelpad'] = 7
mpl.rcParams['xtick.major.pad'] = 7
mpl.rcParams['ytick.major.pad'] = 7
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['xtick.minor.top'] = True
mpl.rcParams['xtick.minor.bottom'] = True
mpl.rcParams['ytick.minor.left'] = True
mpl.rcParams['ytick.minor.right'] = True
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['hatch.linewidth'] = 1

def_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
cosmo = cosmology.FlatLambdaCDM(H0=70,Om0=0.3)

def metalicity(MH):
    grid = np.array([-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.3, 0.5])
    tags = ['m40','m35','m30','m25','m20','m15','m10','m05','m00','p03','p05']
    idx = int(np.argmin((grid - MH)**2))
    print(f"metalicity {tags[idx]}")
    return 'phoenix' + tags[idx]  # 例: 'phoenixm05'

def phoenix_temperature_grid():
    # あなたの元の定義をそのまま使用
    return np.hstack((
        np.linspace(2000, 7000, 51),
        np.linspace(7200, 12000, 25),
        np.linspace(12500, 20000, 16),
        np.linspace(21000, 70000, 50)
    ))

def temperature_bracket(Teff):
    Tgrid = phoenix_temperature_grid()
    # 端ならその温度を両側にして w=0 or 1 にする
    if Teff <= Tgrid.min():
        return float(Tgrid.min()), float(Tgrid[min(1, len(Tgrid)-1)]), 0.0
    if Teff >= Tgrid.max():
        return float(Tgrid[max(0, len(Tgrid)-2)]), float(Tgrid.max()), 1.0
    # 内部なら下側インデックスを探す
    ihi = int(np.searchsorted(Tgrid, Teff))
    ilo = ihi - 1
    Tlo, Thi = float(Tgrid[ilo]), float(Tgrid[ihi])
    w = (Teff - Tlo) / (Thi - Tlo)
    #print(f"lower temperature: {Tlo}, upper temperature: {Thi}, weight: {w}")
    return Tlo, Thi, float(w)

def phoenix_filename(base_dir, metal_tag, T):
    # 例: base/phoenixm05/phoenixm05_5800.fits
    return os.path.join(base_dir, metal_tag, f"{metal_tag}_{int(round(T))}.fits")

def pick_logg_column(filename, logg_target):
    with fits.open(filename) as hdul:
        cols = hdul[1].columns.names
        # 例: 'g40','g45','g50' ... のような列名から数値を抽出
        cand = []
        for nm in cols:
            nm_l = nm.lower()
            if nm_l.startswith('g'):  # 'g45' 等
                s = nm_l[1:]
                if s.replace('.','',1).isdigit() or s.isdigit():
                    try:
                        val = float(s) / 10.0
                        cand.append((nm, val))
                    except Exception:
                        pass
        if not cand:
            # フォーマットが違う場合は元のあなたの方法にフォールバック
            data = hdul[1].data
            avail = []
            vals = []
            for i, nm in enumerate(cols):
                try:
                    if np.max(data[nm]) > 0:
                        avail.append(nm)
                        vals.append(i * 0.5)  # 仮の尺度（並び順依存）
                except Exception:
                    pass
            if not avail:
                raise RuntimeError("No usable logg columns in file.")
            vals = np.array(vals)
            idx = int(np.argmin((vals - logg_target)**2))
            return avail[idx]
        # 近い logg を選ぶ
        names, vals = zip(*cand)
        idx = int(np.argmin((np.array(vals) - logg_target)**2))
        return names[idx]

# ===== コア：温度線形補間つき SED → f_mJy =====

def get_stellar_flux(property):
    """
    property = [Teff[K], logg[cm/s^2のlog10], [M/H], distance_pc, radius_Rsun]
    戻り値: wavelength[Å], f_mJy[同次元]
    """
    teff   = float(property[0])
    logg   = float(property[1])
    mh     = float(property[2])
    dist_pc = float(property[3])   # ここでは『距離[pc]』として扱います
    radius = float(property[4])    # [R_sun]

    base_dir = "/mnt/data_cat3/yuka/data/grp/redcat/trds/grid/phoenix/"
    metal = metalicity(mh)

    # 温度の左右を取得
    Tlo, Thi, w = temperature_bracket(teff)
    fn_lo = phoenix_filename(base_dir, metal, Tlo)
    #print(f"using files:\n  {fn_lo}")
    fn_hi = phoenix_filename(base_dir, metal, Thi)
    #print(f"using files:\n  {fn_hi}")

    # 各テンプレートから最も近い logg 列を選ぶ（ファイルごとに判定）
    col_lo = pick_logg_column(fn_lo, logg)
    #print(f"logg low {col_lo}")
    col_hi = pick_logg_column(fn_hi, logg)
    #print(f"logg low {col_hi}")

    # スペクトル読み込み
    with fits.open(fn_lo) as h0, fits.open(fn_hi) as h1:
        wav_lo = np.array(h0[1].data['WAVELENGTH'], dtype=float)  # [Å]
        wav_hi = np.array(h1[1].data['WAVELENGTH'], dtype=float)
        # 波長グリッドは同一想定。違っていたら共通グリッドへ補間
        if wav_lo.shape != wav_hi.shape or not np.allclose(wav_lo, wav_hi, rtol=0, atol=1e-6):
            # 共通波長（lo側）に合わせる
            from scipy.interpolate import interp1d
            F_hi_raw = np.array(h1[1].data[col_hi], dtype=float)
            interp_hi = interp1d(wav_hi, F_hi_raw, bounds_error=False, fill_value="extrapolate")
            F_lo = np.array(h0[1].data[col_lo], dtype=float)
            F_hi = interp_hi(wav_lo)
            wavelength = wav_lo
        else:
            wavelength = wav_lo
            F_lo = np.array(h0[1].data[col_lo], dtype=float)
            F_hi = np.array(h1[1].data[col_hi], dtype=float)

    # 温度線形補間（表面フラックス F_lambda と仮定：sr^-1 ではない）
    F_surf = (1.0 - w) * F_lo + w * F_hi      # [erg/s/cm^2/Å] at stellar surface

    # 幾何学スケール (R/D)^2 で観測フラックスへ
    R_sun_cm = 6.957e10
    distance_cm = dist_pc * pc
    solid = (radius * R_sun_cm / distance_cm) ** 2
    f_lambda = F_surf * solid                  # [erg/s/cm^2/Å]

    # f_lambda → f_nu → mJy
    cspeed_A = c * 1e8
    f_nu = f_lambda * (wavelength**2) / cspeed_A        # [erg/s/cm^2/Hz]
    f_mJy = (f_nu / Jy) * 1e3    
    return wavelength, f_mJy

def bin_spherex_flux(lam, flux , lam_bins_width_angstrom):

    xgrid = np.arange(0.8,5 , lam_bins_width_angstrom/1e4/2)
    dxgrid = lam_bins_width_angstrom/1e4
    
    LAM_bin = []
    FLUX_bin = []
    FLUXERR_bin = []
    for xx in xgrid:
    
        sel = np.where( (lam >= (xx-dxgrid/2)) & (lam <= (xx+dxgrid/2)) )[0]
    
        if len(sel) > 2:
            lam_bin = np.nanmedian(lam[sel])
            f = flux[sel]
            f = f[~np.isnan(f)]
            msk = sigma_clip(f , sigma=3).mask
            f = f[~msk]
    
            if len(f) > 1:
                #mean, median, stddev = sigma_clipped_stats(f, sigma=3)        
                flux_bin = np.median(f)
                #fluxerr_bin = np.median(stddev)
                fluxerr_bin = bootstrap_median_std(f , num_bootstrap_samples=200)
                
                LAM_bin.append(lam_bin)
                FLUX_bin.append(flux_bin)
                FLUXERR_bin.append(fluxerr_bin)
    
    LAM_bin = np.asarray(LAM_bin)
    FLUX_bin = np.asarray(FLUX_bin)
    FLUXERR_bin = np.asarray(FLUXERR_bin)

    return(LAM_bin, FLUX_bin , FLUXERR_bin)

def get_spherex_cutouts(ra, dec, size_arcsec, output_path, chunk_size, n_processes):
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    img_table = get_spherex_exposures(ra, dec, size_arcsec)
    idx = np.arange(len(img_table))

    def _worker(pars):
        uri_full, ra_, dec_, size_, save_dir = pars
        try:
            return get_spherex_cutout_full(uri_full, ra_, dec_, size_, save_dir, overwrite=True)
        except Exception as e:
            print("cutout failed:", e, uri_full)
            return None

    # タスク作成（uri_full を使うのは OK：フルフレームから ZODI を自前切り出し）
    tasks = [(img_table[i]["uri_full"], ra, dec, size_arcsec, output_path) for i in idx]

    # 並列は I/O 帯域に合わせて（例：4）
    with mp.Pool(processes=n_processes) as pool:
        list(tqdm(pool.imap_unordered(_worker, tasks), total=len(tasks)))

    return img_table


def _pixscale_arcsec(header):
    w = WCS(header)
    try:
        cd = w.pixel_scale_matrix  # deg/pix
        scale_deg = float(np.hypot(cd[0,0], cd[1,1]))
        return scale_deg * 3600.0
    except Exception:
        sx = abs(header.get("CDELT1", 0))*3600.0
        sy = abs(header.get("CDELT2", 0))*3600.0
        return (sx+sy)/2 if (sx>0 and sy>0) else 6.15


def get_spherex_cutout_full(uri, ra_deg, dec_deg, size_arcsec, save_path, overwrite):
    out_name = os.path.basename(uri.split("?")[0]).replace(".fits", "-cutout.fits")
    out_path = os.path.join(save_path, out_name)

    with fits.open(uri, memmap=True) as hdul:
        # 拡張名は大文字想定（小文字で存在判定すると失敗しがち）
        img_hdu = hdul['IMAGE'] if 'IMAGE' in hdul else hdul[1]
        if 'ZODI' not in hdul:
            raise RuntimeError("ZODI extension not found in full frame.")
        zodi_hdu = hdul['ZODI']

        wave_src = hdul['WCS-WAVE']
        wave_hdu = fits.BinTableHDU(data=wave_src.data, header=wave_src.header, name='WCS-WAVE')
        # 親座標系であることと cutout の原点をメモ

        w = WCS(img_hdu.header, fobj=hdul, relax=True)
        coord = SkyCoord(ra_deg*u.deg, dec_deg*u.deg, frame="icrs")
        x, y = w.world_to_pixel(coord)

        size_pix = max(int(round(size_arcsec/6.15)), 3)
        cut_img  = Cutout2D(img_hdu.data,  (x, y), size_pix, wcs=w)
        cut_zodi = Cutout2D(zodi_hdu.data, (x, y), size_pix, wcs=w)

        xc_1b = float(x) + 1.0
        yc_1b = float(y) + 1.0
        wave_hdu.header['WVT_XC']  = (xc_1b, 'Center X in parent frame (1-based, subpixel)')
        wave_hdu.header['WVT_YC']  = (yc_1b, 'Center Y in parent frame (1-based, subpixel)')
        wave_hdu.header['WVT_IXC'] = (int(round(float(xc_1b))) + 1, 'Nearest integer X (1-based)')
        wave_hdu.header['WVT_IYC'] = (int(round(float(yc_1b))) + 1, 'Nearest integer Y (1-based)')

    h0 = fits.PrimaryHDU()
    h1 = fits.ImageHDU(data=cut_img.data,  header=cut_img.wcs.to_header(),  name='IMAGE')
    hZ = fits.ImageHDU(data=cut_zodi.data, header=cut_zodi.wcs.to_header(), name='ZODI')
    fits.HDUList([h0, h1, hZ, wave_hdu]).writeto(out_path, overwrite=overwrite)
    return out_path

    
def get_spherex_cutout_helper(uri, uri_zodi, save_path, uri_full, overwrite):
    '''
    Creates cutouts for a given access URI from the table 
    obtained by `get_spherex_exposures()`.

    access_uri: Directly `access_uri` column from table.
    size_arcsec: Cutout size in arcseconds
    save_path: Path where to save the cutouts
    overwrite: Overwrite existing cutouts (True / False)
    '''

    out_name = uri.split("?")[0].split("/")[-1].replace(".fits","-cutout.fits")
    #uri = "https://irsa.ipac.caltech.edu/ibe/cutout?ra={}&dec={}&size={}&path={}".format(ra,dec,size_deg,uri_path)
    #print(uri)

    zodi_name = uri_zodi.split("?")[0].split("/")[-1].replace(".fits","-cutout_zodi.fits")

    try:
        with fits.open(uri) as hdul:
            filename = os.path.join(save_path , out_name)
            hdul.writeto(filename , overwrite=overwrite)
        
        with fits.open(uri_zodi) as hdul:
            filename = os.path.join(save_path , zodi_name)
            hdul.writeto(filename , overwrite=overwrite)
        return(True)

    except Exception as e:
        print(e)
        return(False)

def split_by_chunk_size(items, chunk_size):
    """
    Splits a list into multiple lists, each with a maximum length of chunk_size.

    Parameters:
    items (list): The list to be split.
    chunk_size (int): The maximum length of each sublist.

    Returns:
    list of lists: A list containing sublists of specified length.
    """
    if chunk_size <= 0:
        raise ValueError("Chunk size must be greater than 0")
    
    return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]

def get_spherex_exposures(ra,dec,size_arcsec):
    '''
    Collects all SPHEREx exposures for a given coordinate (ra,dec)
    '''

    ## Tip:
    # To check all IRSA collections:
    #[print(c[0]) for c in Irsa.list_collections()]

    coord = SkyCoord(ra , dec, unit = "deg" , frame = "icrs")
    size_degree = size_arcsec / 3600 # in deg
    
    # Define the service endpoint for IRSA's Table Access Protocol (TAP)
    # so that we can query SPHEREx metadata tables.
    service = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
    
    # Define a query that will search the appropriate SPHEREx metadata tables
    # for spectral images that cover the chosen coordinates and match the
    # specified bandpass. Return the cutout data access URL and the time of observation.
    # Sort by observation time.
    query = f"""
    SELECT
        'https://irsa.ipac.caltech.edu/' || a.uri || '?center={coord.ra.degree},{coord.dec.degree}d&size={size_degree}' AS uri,
        'https://irsa.ipac.caltech.edu/' || a.uri || '?center={coord.ra.degree},{coord.dec.degree}d&size={size_degree}&ext=ZODI' AS uri_zodi,
        p.time_bounds_lower,
        'https://irsa.ipac.caltech.edu/' || a.uri AS uri_full
    FROM spherex.artifact a
    JOIN spherex.plane p ON a.planeid = p.planeid
    WHERE 1 = CONTAINS(POINT('ICRS', {coord.ra.degree}, {coord.dec.degree}), p.poly)
    ORDER BY p.time_bounds_lower
    """
    
    # Execute the query and return as an astropy Table.
    result = service.search(query)
    img_table = result.to_table()
    
    
    '''## Need to query all spherex collections separately
    img_table_allsky = Irsa.query_sia(pos=(coord, 1 * u.arcsec), collection='spherex_qr', maxrec=10000000).to_table()
    img_table_deep = Irsa.query_sia(pos=(coord, 1 * u.arcsec), collection='spherex_qr_deep', maxrec=10000000).to_table()
    #img_table_cal = Irsa.query_sia(pos=(coord, 1 * u.arcsec), collection='spherex_qr_cal', maxrec=10000000).to_table()

    ## Then combine them
    img_table = vstack([img_table_allsky,img_table_deep])'''

    
    print("Number of SPHEREx exposure: {}".format(len(img_table)))

    return(img_table)

def bootstrap_median_std(X, num_bootstrap_samples=1000):
    medians = []
    n = len(X)
    for _ in range(num_bootstrap_samples):
        sample = np.random.choice(X, size=n, replace=True)
        medians.append(np.median(sample))
    return np.std(medians)


def measure_spherex_flux(fns , aperture_radius_px, n_processes, chunk_size):
    '''
    Measured the flux (Aperture) at the center of the cutouts
    '''

    pool = mp.Pool(processes=n_processes)
    pars = [(fn,aperture_radius_px) for fn in fns]
    results = pool.map(measure_spherex_flux_helper, pars, chunksize=chunk_size)
    
    LAM = np.asarray( [res[0] for res in results] )
    FLUX = np.asarray( [res[1] for res in results] )
    FLUX_SUB = np.asarray( [res[2] for res in results] )

    return(LAM , FLUX, FLUX_SUB)

def get_spherex_lambda(hdul, x, y):
    '''
    Get the wavelength and width from an hdul at pixel position (x,y)
    '''

    img = hdul['IMAGE'].data
    hdr = hdul['IMAGE'].header

    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f): # suppress print()
        spectral_wcs = WCS(hdr, fobj=hdul, key="W")
    spectral_wcs.sip = None
            
    ## Get Lambda
    lam, bp = spectral_wcs.pixel_to_world(img.shape[0]//2 , img.shape[1]//2)
    lam = lam.value
    dlam = bp.value

    return(lam , dlam)
    

def mjsr_to_jypixel(value_mjsr, pixel_size_arcsec):
    """
    Convert surface brightness from MJy/sr to Jy/pixel.

    Parameters
    ----------
    value_mjsr : float or array-like
        Value(s) in MJy/sr.
    pixel_size_arcsec : float
        Pixel size in arcseconds (assumed square pixels).

    Returns
    -------
    float or ndarray
        Equivalent value(s) in Jy/pixel.
    """
    # Constants
    arcsec_to_rad = np.pi / (180.0 * 3600.0)  # 1 arcsec in radians

    # Convert pixel size from arcsec^2 to steradian
    pixel_area_sr = (pixel_size_arcsec * arcsec_to_rad)**2

    # 1 MJy/sr = 1e6 Jy/sr
    value_jy_sr = value_mjsr * 1e6

    # Multiply by pixel solid angle
    value_jy_pixel = value_jy_sr * pixel_area_sr

    return value_jy_pixel

def get_lambda(data, header):
    x_axis = data[0][0]
    y_axis = data[0][1]
    values = data[0][2][:,:,0]


    x_center = header["WVT_XC"]
    y_center = header["WVT_YC"]

    f = RegularGridInterpolator((y_axis, x_axis), values, method='linear',
                                bounds_error=False, fill_value=np.nan)
        
    lam = float(f((y_center + 1, x_center+1)))
    return(lam)

def measure_spherex_flux_helper(pars):

    fn = pars[0]
    aperture_radius_px = pars[1]

    ## Load
    with fits.open(fn) as hdul:
        img = hdul['IMAGE'].data
        hdr = hdul['IMAGE'].header

        zodi = hdul['ZODI'].data

        wave_hdu = hdul['WCS-WAVE']
        lam = get_lambda(wave_hdu.data, wave_hdu.header)

        #with open(os.devnull, "w") as f, contextlib.redirect_stdout(f): # suppress print()
            #spectral_wcs = WCS(hdr, fobj=hdul, key="W")
        #spectral_wcs.sip = None

    ## Convert flux MJy/sr to Jy/px
    img = mjsr_to_jypixel(value_mjsr = img, pixel_size_arcsec = 6.15)
    zodi = mjsr_to_jypixel(value_mjsr = zodi, pixel_size_arcsec = 6.15)
    #img = img * 2.350443e-5 # MJy/sr -> Jy/arcsec2
    #img = img * (6.15**2)# Jy/arcsec2 -> Jy/px
    img = img * 1e3 # Jy/px -> mJy/px
    zodi = zodi * 1e3
    
    ## Measure flux (this is super simple currently)
    #flux = np.nansum(img)

    ## Measure aperture flux
    positions = [(img.shape[0]/2.0, img.shape[1]/2.0)]
    #mean, median, stddev = sigma_clipped_stats(img[~np.isnan(img)] , sigma=3 , maxiters=5)
    aperture = CircularAperture(positions, r=aperture_radius_px)
    #aphot = aperture_photometry(img - median , aperture, subpixels=10)
    #aphot_zodi = aperture_photometry(img - median - zodi , aperture, subpixels=10)
    aphot = aperture_photometry(img , aperture, subpixels=10)
    aphot_zodi = aperture_photometry(img - zodi , aperture, subpixels=10)
    flux = aphot["aperture_sum"][0]
    flux_sub = aphot_zodi["aperture_sum"][0]
    #print("Flux: {:2.3f} mJy | Lambda: {:2.3f} um | File: {}".format(flux, lam, fn))
        

    return(lam, flux, flux_sub)


def spherex_photo_pipeline(coord,obj_name, redshift, cutout_size_arcsec, aperture_radius_px, main_path,
                           CREATECUTOUTS, OVERWRITE, MAKEFIGURE, n_processes, chunk_size, lam_bins_width_ang, propertys=None):

    ## Output path:
    spherex_main = os.path.join(main_path , "spherex_data" , obj_name)# "./data/{}".format(obj_name)
    
    ## Create Cutouts
    RUNCUTOUTS = False
    if CREATECUTOUTS:
        if (os.path.exists(spherex_main)) & (OVERWRITE):
            RUNCUTOUTS = True
        elif not os.path.exists(spherex_main):
            RUNCUTOUTS = True
        else:
            print("Directory {} exists. Because OVERWRITE = False, I do not create new cutouts.".format(spherex_main))
            RUNCUTOUTS = False
    else:
        RUNCUTOUTS = False
    
    if RUNCUTOUTS:
        results = get_spherex_cutouts(ra = coord.ra.degree, dec = coord.dec.degree,
                                      size_arcsec = cutout_size_arcsec,
                                      output_path = spherex_main,
                                      chunk_size = chunk_size,
                                      n_processes = n_processes
                                      )

    
    ## MEASURE APERTURE PHOTOMETRY ======
    
    ## Load images
    fns = glob.glob(os.path.join(spherex_main , "*.fits"))
    print("Number of images: {}".format(len(fns)))
    
    ## Measure Aperture photometry
    print("Measuring photometry . . .")
    LAM, FLUX, FLUX_SUB = measure_spherex_flux(fns ,
                                     aperture_radius_px = aperture_radius_px,
                                     n_processes = n_processes,
                                     chunk_size = chunk_size)
    
    ## Create He mask
    he_mask = (LAM > 1.05) & (LAM < 1.15)
    
    ## Apply mask
    LAM = LAM[~he_mask]
    FLUX = FLUX[~he_mask]
    FLUX_SUB = FLUX_SUB[~he_mask]


    ## BIN THE PHOTOMETRY ========
    print("Binning Data")
    LAM_bin, FLUX_bin , FLUXERR_bin = bin_spherex_flux(lam = LAM,
                                                       flux = FLUX ,
                                                       lam_bins_width_angstrom = lam_bins_width_ang)
    
    LAM_bin, FLUX_SUB_bin , FLUXERR_SUB_bin = bin_spherex_flux(lam = LAM,
                                                       flux = FLUX_SUB ,
                                                       lam_bins_width_angstrom = lam_bins_width_ang)


    ## SAVE THE PHOTOMETRY ========
    path_spectrum = "./spectra/{}/".format(obj_name)
    if not os.path.exists(path_spectrum):
        print("Path {} does not exist. Creating it".format(path_spectrum))
        os.mkdir(path_spectrum)
    
    ## Get theoretical stellar spectrum(if possible) =====
    if propertys is not None:
        print("Getting theoretical stellar spectrum")
        wavelength, f_mJy = get_stellar_flux(propertys)
        #wavelength, f_mJy_old = get_stellar_flux_old(propertys)
    
    ## Create object that will be saved as pickel
    tab1 = Table([LAM,FLUX, FLUX_SUB], names=["lam_int","flux_int", "flux_sub"])
    tab2 = Table([LAM_bin,FLUX_bin, FLUX_SUB_bin, FLUXERR_SUB_bin], names=["lam_bin","flux_bin", "flux_sub_bin", "fluxerr_sub_bin"])
    
    hdu0 = fits.PrimaryHDU()
    hdu1 = fits.BinTableHDU(tab1, name='TABLE1')
    hdu2 = fits.BinTableHDU(tab2, name='TABLE2')
    
    hdul = fits.HDUList([hdu0, hdu1, hdu2])
    hdul.writeto(os.path.join(path_spectrum,"{}_spherex.fits".format(obj_name) ), overwrite=True)

    ## PLOT FINAL FIGURE ####
    if MAKEFIGURE:
        print("Making figure . . .")
        fig = plt.figure(figsize=(9,4))
        ax1 = fig.add_subplot(1,1,1)
        
        xlims = [0.3-0.1,5+0.1]
        
        ## SPHEREx spectrum ==========
        ax1.scatter(LAM, FLUX , marker="o", color = "black", alpha=0.2, s=1, zorder=100, label="SPHEREx (native)")
        ax1.scatter(LAM, FLUX_SUB , marker="o", color = "blue", alpha=0.2, s=1, zorder=100, label="SPHEREx (foreground subtracted)")
        ax1.errorbar(LAM_bin , FLUX_bin, yerr=FLUXERR_bin, fmt="o",zorder=101,
                             markerfacecolor="black",markeredgecolor="white",alpha=1, markersize=4,
                     elinewidth=1, linewidth=1,capsize=2, ecolor="black", label="SPHEREx (native)" )
        ax1.errorbar(LAM_bin , FLUX_SUB_bin, yerr=FLUXERR_SUB_bin, fmt="o",zorder=101,
                             markerfacecolor="blue",markeredgecolor="white",alpha=1, markersize=4,
                     elinewidth=1, linewidth=1,capsize=2, ecolor="blue", label="SPHEREx (foreground subtracted)" )
        if propertys is not None:
            ax1.plot(wavelength/1e4, f_mJy , "-", color="red", alpha=0.7, linewidth=1, label="Model")
            #ax1.plot(wavelength/1e4, f_mJy_old , "-", color="green", alpha=0.7, linewidth=1, label="Model old")
            rp_mag = propertys[5]
            rp_Jy = 10**(-rp_mag*0.4)*3631*1e3
            ax1.errorbar(0.85 , rp_Jy, xerr=0.2, fmt="o",zorder=101,
                         markerfacecolor="purple",markeredgecolor="white",alpha=1, markersize=10,
                         elinewidth=1, linewidth=1,capsize=2, ecolor="purple", label=f"flux {rp_Jy}mJy" )
            ax1.set_title(f"{obj_name} Teff: {propertys[0]} logg:{propertys[1]} M/H:{propertys[2]}", fontsize=8)
        else:
            ax1.set_title(f"{obj_name}", fontsize=10)
            #ax1.set_ylim(0, rp_Jy*1.5)

        
        
        #ylims = [5-0.2,17]
        #ylims = [ np.round( np.min( sigma_clip(FLUX_bin, sigma=3, masked=False).data)*0.8*10)/10 ,
        #        np.round( np.max( sigma_clip(FLUX_bin, sigma=3, masked=False).data)*1.2*10)/10
        #        ]
        #ylims = [ np.round( np.min( sigma_clip(FLUX_bin, sigma=3, masked=False).data)*0.2*10)/10 ,
        #        np.round( np.max( sigma_clip(FLUX_bin, sigma=3, masked=False).data)*2*10)/10
        #        ]
        ylims = [ np.min(FLUX_bin) ,
                np.max(FLUX_bin)
                ]
        ylims[1] *= 1.5
        if ylims[0] < 0: ylims[0] *= 2
        if ylims[0] > 0: ylims[0] /= 2
        
        ## DESI spectrum:
        spec_desi_path = os.path.join(main_path , "spectra/{}".format(obj_name) , "{}_desi.fits".format(obj_name) )
        if os.path.exists(spec_desi_path):
            spectab_desi = Table.read(spec_desi_path , hdu=2)
        
            ax1.plot(spectab_desi["wavelength"]/1e4, spectab_desi["fnu_mjy"] , "-", color=def_cols[0], alpha=0.5, linewidth=0.5, label="DESI")
        else:
            print("no DESI spectrum available")
        
        ## indicate He feature ========
        he_lam = [1.05,1.15]
        ax1.axvspan(he_lam[0], he_lam[1], alpha=0.2, color='gray', linewidth=0)
        
        ## Emission lines =========
        if redshift <= 0:
            redshift = 0
        lines = [0.6563, 1.282, 1.6, 1.875, 3.3, 4.05]
        line_names = [r"H$\alpha$",r"Pa-$\beta$",r"$1.6\,{\rm \mu m}$",r"Pa-$\alpha$", "PAH\n({})".format(r"$3.3\,{\rm \mu m}$") , r"Br-$\alpha$"]
        h = 0.8
        for l,ln in zip(lines,line_names):
            if (l*(1+redshift) > xlims[0]) & (l*(1+redshift) < xlims[1]):
                ax1.axvline(l*(1+redshift), ymin=0.02, ymax=h, linewidth=0.5, linestyle="--", color="gray")
                ax1.text(l*(1+redshift) , (ylims[0]+np.diff(ylims)*(h*1.01)) , ln, va="bottom", ha="center",
                         fontsize=11, rotation=0, color="black")
        
        # Indicate redshift
        ax1.text(0.05,0.95 , r"z = {:2.3f}".format(redshift), va="top", ha="left", transform=ax1.transAxes,
                     fontsize=13, color="black")
        
        leg = ax1.legend(loc="upper left", fontsize=11,ncol=1,bbox_to_anchor=(1,1),
                         frameon=True, numpoints=1,scatterpoints=5, markerscale=1.2)
        leg.legend_handles[0].set_alpha(0.7)
        #leg.legend_handles[0]._sizes = [10]
        #leg.legend_handles[1]._sizes = [20]
        
        #ax1.set_title("euclid-{} | zspec = {:2.3f} | zphot = {:2.3f} | Hmag = {:2.2f}".format(obj_id,zspec,zphot,hmag) , fontsize=10)
        
        #ax1.set_yscale('log')
        ax1.set_xlim(xlims[0], xlims[1])
        #ax1.set_ylim(ylims[0], ylims[1])
        ax1.set_ylim(0, ylims[1])
        ax1.set_xlabel(r"Observed Wavelength ($\rm \mu m$)")
        ax1.set_ylabel(r"Flux ($\rm mJy$)")
        
        plt.savefig(os.path.join(main_path , "plots/{}.pdf".format(obj_name)) , bbox_inches="tight")
        plt.show()

        return(tab1 , tab2)
