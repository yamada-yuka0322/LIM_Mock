import numpy as np
import os

import pandas as pd

from scipy.interpolate import interp1d
from scipy.interpolate import interpn
from scipy.integrate import cumulative_trapezoid

from astropy.cosmology import Planck18 as cosmo
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic
from astropy.coordinates import FK5
from astropy.table import Table

import healpy as hp

import matplotlib.pyplot as plt

import zodipy

import datetime

import multiprocessing

#from DustEM_loader import *
from lim.create_JWST import frequency as Frequency

from functools import lru_cache

import h5py

Oiii = 5.997e5 #[GHz]
ha = 4.856e5 #[GHz]
pa = 1.599e5 #[GHz]
c = 2.9979e10 #[cm/s]
Jy = 1.0e-23
arcsec = 4.848136811094e-6
pc = 3085677857083127000
Rsun_cm = 6.957e10
pc_cm = 3.0856775814913673e18

def make_foreground(ra, dec, params, zodi_model = "dirbe"):
    flist, dflist = Frequency(params) # Hz, Hz
    star_path = '/mnt/data_cat3/yuka/data/LIM_mock/stars/output.fits'
    star, mask_ra, mask_dec = Star(ra, dec, params, flist, dflist, star_path)
    dust_path = '/mnt/data_cat3/yuka/repository/DustEM'
    dgl = DGL(ra, dec, params, flist, dust_path, makefig=True)
    mask, star_mask = Mask(ra, dec, params, mask_ra, mask_dec)
    zodiac, zodiac1 = Zodiac(coarse=32, center_ra=ra, center_dec=dec, params=params, flist=flist, zodi_model = zodi_model, day_bin=2, count=100, nprocesses=None)
    
    total = (star + dgl + (zodiac+zodiac1)/2.0) * (mask*star_mask)[:, :, np.newaxis]
    
    return star, dgl, [zodiac, zodiac1], mask, star_mask, total

#####################################################
def Star(center_ra, center_dec, params, flist, dflist, path):
    #いくつか出力
    Nx = int(3600*params.radius / params.resolution)
    Nz = len(flist)
    star_intensity = np.zeros((Nx, Nx, Nz), dtype=np.float32)
    
    ra_min = center_ra - params.radius
    dec_min = center_dec - params.radius
    
    hdu = fits.open(path)
    data = hdu[1].data
    J = data['J']
    logg = data['logg']
    ra = data['RAJ2000']
    ra = np.array(ra)
    dec = data['DECJ2000']
    dec = np.array(dec)
    
    distance = np.sqrt((ra - center_ra)**2+(dec - center_dec)**2)
    #area = (distance < 7.5)
    area = (np.abs(ra - center_ra) < params.radius/2.0)&(np.abs(dec - center_dec) < params.radius/2.0)
    bright = (J<22) & area
    mask_ra = ra[bright]
    mask_dec = dec[bright]
    
    stars = (logg<6) & area
    
    hdu0 = fits.PrimaryHDU()
    
    _path = '/mnt/data_cat3/yuka/data/grp/redcat/trds/grid/phoenix/'
    
    bright_stars = pd.DataFrame({})
    
    bright_star_ra = []
    bright_star_ira = []
    bright_star_dec = []
    bright_star_idec = []
    bright_star_path = []
    bright_star_logg = []
    
    ix = []
    iy = []
    spectrum = []
    T_eff = data[stars]['Teff']
    metallicity = data[stars]['[M/H]']
    _logg = logg[stars]
    
    
    for i in data[stars]:
        _ra = i['RAJ2000']
        _dec = i['DECJ2000']
        _J = i['J']
        
        i_ra = (_ra - ra_min) / params.resolution *3600
        i_ra = i_ra.astype(np.int32)
        iy.append(i_ra)
        i_dec = (_dec - dec_min) / params.resolution *3600
        i_dec = i_dec.astype(np.int32)
        ix.append(i_dec)
        
        if (i_ra<Nx)&(i_dec<Nx):
            filename = _path
            
            dist = i['dist']
            MH = i['[M/H]']
            filename = filename + metalicity(MH) + '/'
            
            Teff = i['Teff']
            
            filename = filename + Temperature(Teff, metalicity(MH))
            
            column = Logg(i['logg'], filename)
            
            radius = i['radius'] #stellar radius [R_sun]
            dist = i['dist'] #distance [kpc]
            
            hdu1 = fits.open(filename)
            
            data1 = hdu1[1].data
            wavelength = data1['WAVELENGTH'] #[Å]
            flux = data1[column] #surface flux [erg/s/cm^2/sr/Å]
            angle = (radius*6.957*10**10)/(dist*1000.0*pc)
            flux *= angle**2 * np.pi # total flux [erg/s/cm^2/Å]
            
            f_nu = compute_tophat_flux(wavelength, flux, flist) # flux f_nu [Jy]
            F = f_nu / (params.resolution * arcsec)**2 #[Jy/sr]
            spectrum.append(F)
            
            star_intensity[i_dec, i_ra] += F
            if(_J<16):
                bright_star_ra.append(_ra)
                bright_star_dec.append(_dec)
                bright_star_ira.append(i_ra)
                bright_star_idec.append(i_dec)
                bright_star_path.append(filename)
                bright_star_logg.append(column)
                
    spectrum = np.array(spectrum)
    
    outdir = '../output/catalog'
                
    save_meta_parquet(outdir, ra[stars], dec[stars], ix, iy, metallicity, T_eff, _logg, flist)

    # 3) HDF5 保存（連続光＋輝線＋合成）
    save_spectra_h5(outdir, flist, dflist, spectrum)
    
    bright_stars['ra'] = bright_star_ra
    bright_stars['i_ra'] = bright_star_ira
    bright_stars['dec'] = bright_star_dec
    bright_stars['i_dec'] = bright_star_idec
    bright_stars['path'] = bright_star_path
    bright_stars['logg'] = bright_star_logg
    
    table = Table.from_pandas(bright_stars)
    table.write('../output/bright_stars.fits', format='fits', overwrite=True)
    
    return star_intensity, mask_ra, mask_dec
            

            
def metalicity(MH):
    metalicities = np.array([-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.3, 0.5])
    name = ['m40', 'm35', 'm30', 'm25', 'm20', 'm15', 'm10', 'm05', 'm00', 'p03', 'p05']
    closest = np.argmin((metalicities - MH)**2)
    name = 'phoenix' + name[closest]
    return name

def Temperature(Teff, metalicity):
    _T = np.hstack((np.linspace(2000, 7000, 51), np.linspace(7200, 12000, 25), np.linspace(12500, 20000, 16), np.linspace(21000, 70000, 50)))
    closest = int(_T[np.argmin((_T - Teff)**2)])
    name = f'{metalicity}_{closest}.fits'
    return name


def Logg(logg, filename):
    hdu = fits.open(filename)
    data = hdu[1].data
    loggs = []
    name = data.columns.names
    available = []
    
    for i,g in enumerate(name):
        if np.max(data[g])>0:
            available.append(g)
            loggs.append(i*0.5)
    loggs = np.array(loggs)
    closest = np.argmin((loggs - logg)**2)
    hdu.close()
    return available[closest]

def compute_tophat_flux(wavelength, flux_lambda, frequency):
    """
    Parameters
    wavelength: wavelength of the SED template in Å
    flux_lambda:f_lambda [erg/s/cm^2/Å]
    frequency: frequency bin of the observation [Hz]
    """
    cspeed = c * 1e8 #speed of light [Å/s]
    
    freq_sed = cspeed / wavelength #frequency bin of SED templateHz
    f_nu_sed = flux_lambda * wavelength**2 / cspeed #[erg/s/cm^2/Hz]
    
    f_nu_interp = interp1d(freq_sed[::-1], f_nu_sed[::-1], bounds_error=False, fill_value=0.0)
    
    f_vals = f_nu_interp(frequency)
    f_Jy = f_vals / Jy
    
    return f_Jy
    #f_nu_list = []
    
    #for i in range(len(frequency) - 1):
        #nu1 = frequency[i]
        #nu2 = frequency[i+1]
        
        #nu_grid = np.linspace(nu1, nu2, 100)
        #f_vals = f_nu_interp(nu_grid)
        
        #f_bin_avg = np.trapz(f_vals, nu_grid) / (nu2 - nu1) #[erg/s/cm^2/Hz]
        #f_Jy = f_bin_avg /Jy
        #f_nu_list.append(f_Jy)
        
    #return np.array(f_nu_list)

def _bandavg_Fnu_cgs_from_surface_Flambda(wavelength_A, Fsurf_lambda_per_sr, flist_Hz, dflist_Hz):
    """
    wavelength_A: (Nwave,) Å
    Fsurf_lambda_per_sr: (Nwave,)  PHOENIX surface flux [erg/s/cm^2/sr/Å]
    Returns:
      Fnu_base_surface: (Nz,)  stellar-surface hemispheric flux density [erg/s/cm^2/Hz]
        (surface per area; already integrated over hemisphere via *pi)
    """
    wl = np.asarray(wavelength_A, dtype=np.float64)
    Fsurf = np.asarray(Fsurf_lambda_per_sr, dtype=np.float64)

    # hemispheric surface flux per Å
    F_lambda_surface = Fsurf * np.pi  # erg/s/cm^2/Å

    # nu grid (Hz)
    nu = c / (wl * 1e-8) #Hz
    order = np.argsort(nu)  # ascending nu
    nu = nu[order]
    wl = wl[order]
    F_lambda_surface = F_lambda_surface[order]

    ca = c * 1e8
    Fnu_surface = F_lambda_surface * wl**2 / ca  # erg/s/cm^2/Hz

    # cumulative integral over nu
    I = cumulative_trapezoid(Fnu_surface, nu, initial=0.0)

    nu1 = flist_Hz - 0.5 * dflist_Hz
    nu2 = flist_Hz + 0.5 * dflist_Hz

    I1 = np.interp(nu1, nu, I, left=I[0], right=I[-1])
    I2 = np.interp(nu2, nu, I, left=I[0], right=I[-1])

    dnu = np.clip(nu2 - nu1, 1e-60, None)
    return (I2 - I1) / dnu  # erg/s/cm^2/Hz  (at stellar surface)


@lru_cache(maxsize=1024)
def phoenix_base_Fnu_surface(filename, column, flist_tuple, dflist_tuple):
    """
    Cache per (filename, column): returns base surface Fnu (Nz,)
    """
    flist = np.asarray(flist_tuple, dtype=np.float64)
    dflist = np.asarray(dflist_tuple, dtype=np.float64)

    with fits.open(filename, memmap=True) as h:
        d = h[1].data
        wl = np.asarray(d["WAVELENGTH"], dtype=np.float64)   # Å
        Fsurf = np.asarray(d[column], dtype=np.float64)      # erg/s/cm^2/sr/Å

    return _bandavg_Fnu_cgs_from_surface_Flambda(wl, Fsurf, flist, dflist)  # (Nz,)


def Star_save_perstar_spectra(
    center_ra, center_dec, params, flist, dflist, path,
    phoenix_root="/mnt/data_cat3/yuka/data/grp/redcat/trds/grid/phoenix/",
    outdir="/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/stars",
    faint_star_J=22.0,
    bright_star_J=16.0,
):
    """
    Save per-star spectra (Nstar, Nz) in erg/s/cm^2/Hz.

    Outputs:
      - {outdir}/stars_meta.parquet
      - {outdir}/stars_spectra.h5 (via your save_spectra_h5; dataset name fnu_cont_jy is just a name)
      - {outdir}/bright_stars.fits
    Returns:
      spectra: (Nstar, Nz) float32  [erg/s/cm^2/Hz]
      meta_df: pandas DataFrame
      mask_ra, mask_dec: (for mask candidates)
    """

    os.makedirs(outdir, exist_ok=True)

    flist = np.asarray(flist, dtype=np.float64)
    dflist = np.asarray(dflist, dtype=np.float64)
    Nz = len(flist)

    ra_min = center_ra - params.radius/2.0
    dec_min = center_dec - params.radius/2.0
    Nx = int(3600 * params.radius / params.resolution)

    # --- read catalog once ---
    with fits.open(path, memmap=True) as h:
        data = h[1].data
        J = np.asarray(data["J"])
        logg_all = np.asarray(data["logg"])
        ra = np.asarray(data["RAJ2000"], dtype=np.float64)
        dec = np.asarray(data["DECJ2000"], dtype=np.float64)

        area = (np.abs(ra - center_ra) < params.radius/2.0) & (np.abs(dec - center_dec) < params.radius/2.0)

        bright_for_mask = (J < faint_star_J) & area
        mask_ra = ra[bright_for_mask]
        mask_dec = dec[bright_for_mask]

        stars = (logg_all < 6) & area

        # pull fields for stars
        ra_s = ra[stars]
        dec_s = dec[stars]
        J_s = J[stars]
        Teff_s = np.asarray(data["Teff"][stars])
        MH_s = np.asarray(data["[M/H]"][stars])
        logg_s = np.asarray(data["logg"][stars])
        radius_s = np.asarray(data["radius"][stars])  # R_sun
        dist_s = np.asarray(data["dist"][stars])      # kpc

    # --- pixel indices (vectorized) ---
    iy = ((ra_s - ra_min) / params.resolution * 3600.0).astype(np.int32)   # ra-index
    ix = ((dec_s - dec_min) / params.resolution * 3600.0).astype(np.int32) # dec-index
    inside = (iy >= 0) & (iy < Nx) & (ix >= 0) & (ix < Nx)

    ra_s = ra_s[inside]; dec_s = dec_s[inside]; J_s = J_s[inside]
    Teff_s = Teff_s[inside]; MH_s = MH_s[inside]; logg_s = logg_s[inside]
    radius_s = radius_s[inside]; dist_s = dist_s[inside]
    ix = ix[inside]; iy = iy[inside]

    Nstar = len(ra_s)
    if Nstar == 0:
        # save empty outputs if you want
        empty_meta = pd.DataFrame()
        return np.zeros((0, Nz), dtype=np.float32), empty_meta, mask_ra, mask_dec

    # --- build template keys (filename, column) for each star ---
    filenames = np.empty(Nstar, dtype=object)
    columns = np.empty(Nstar, dtype=object)

    # (ここは metalicity/Temperature/Logg の実装次第で I/O ではなく文字列生成なので軽い)
    for i in range(Nstar):
        mh = MH_s[i]
        teff = Teff_s[i]
        lg = logg_s[i]

        subdir = metalicity(mh)  # <- your function
        folder = os.path.join(phoenix_root, subdir)
        fname = os.path.join(folder, Temperature(teff, subdir))  # <- your function
        col = Logg(lg, fname)  # <- your function

        filenames[i] = fname
        columns[i] = col

    # --- geometric scaling (R/D)^2 ---
    R_cm = radius_s.astype(np.float64) * Rsun_cm
    D_cm = dist_s.astype(np.float64) * 1000.0 * pc_cm
    scale = (R_cm / D_cm)**2  # dimensionless

    # --- cache keys for lru_cache (hashable) ---
    flist_t = tuple(flist.tolist())
    dflist_t = tuple(dflist.tolist())

    # --- group by template key to minimize PHOENIX reads ---
    # pandas.factorize is fast for object pairs
    keys = pd.Series(list(zip(filenames, columns)))
    gid, uniq = pd.factorize(keys, sort=False)

    spectra = np.empty((Nstar, Nz), dtype=np.float32)

    for g, (fname, col) in enumerate(uniq):
        sel = (gid == g)
        if not np.any(sel):
            continue
        try:
            base = phoenix_base_Fnu_surface(fname, col, flist_t, dflist_t)  # (Nz,) surface Fnu
        except Exception as e:
            print(f"[WARN] PHOENIX load failed: {fname} col={col} err={e}")
            # fill zeros for this group
            spectra[sel, :] = 0.0
            continue

        # observed per-star Fnu: base * (R/D)^2
        spectra[sel, :] = (scale[sel, None] * base[None, :]).astype(np.float32)

    # --- save meta ---
    meta_df = pd.DataFrame({
        "star_id": np.arange(Nstar, dtype=np.int64),
        "ra": ra_s.astype(np.float64),
        "dec": dec_s.astype(np.float64),
        "J": J_s.astype(np.float32),
        "Teff": Teff_s.astype(np.float32),
        "MH": MH_s.astype(np.float32),
        "logg": logg_s.astype(np.float32),
        "radius_Rsun": radius_s.astype(np.float32),
        "dist_kpc": dist_s.astype(np.float32),
        "ix": ix.astype(np.int32),
        "iy": iy.astype(np.int32),
        "template_path": filenames.astype(str),
        "template_col": columns.astype(str),
        "Nz": int(Nz),
        "nu_min_Hz": float(flist.min()),
        "nu_max_Hz": float(flist.max()),
    })
    meta_path = os.path.join(outdir, "stars_meta_all.parquet")
    meta_df.to_parquet(meta_path, index=False)

    # --- save spectra (erg/s/cm^2/Hz) ---
    # NOTE: your save_spectra_h5 names the dataset "fnu_cont_jy" but it will store whatever you pass.
    # If you want, you can later rename dataset key; for now, keep compatibility.
    save_spectra_h5(outdir, flist, dflist, spectrum=spectra)

    # --- save bright star list (optional) ---
    #bright = (J_s < bright_star_J)
    #if np.any(bright):
        #bright_df = meta_df.loc[bright, ["ra", "dec", "ix", "iy", "template_path", "template_col", "J"]].copy()
        #Table.from_pandas(bright_df).write(os.path.join(outdir, "bright_stars.fits"),
                                           #format="fits", overwrite=True)

    return spectra, meta_df, mask_ra, mask_dec

############################################################################################################
@lru_cache(maxsize=8)
def _load_SED(path):
    # SED.res: skiprows=8 は固定仕様
    SED = np.loadtxt(path, skiprows=8)
    wave_um = SED[:, 0].astype(np.float64)
    nu      = c / (wave_um * 1e-4)
    nuI_nu  = SED[:, -1].astype(np.float64)            # 4π ν Iν / NH [erg/s/H]
    I_nu_over_NH = nuI_nu / nu / (4*np.pi)             # Iν / NH [erg/s/Hz/sr/H]
    return wave_um, I_nu_over_NH

@lru_cache(maxsize=8)
def _load_ISRF(path):
    arr = np.loadtxt(path, skiprows=7)
    wave_um = arr[:, 0].astype(np.float64)
    fourpi_I = arr[:, 1].astype(np.float64) / (4*np.pi)           # Iν [erg/cm^2/sr/s/Hz]
    # return np.interp 関数化版（scipy不要の場合）
    # def f(x):
    #     return np.interp(x, wave_um, fourpi_I, left=0.0, right=0.0)
    # return f
    return interp1d(wave_um, fourpi_I, kind='linear',
                    bounds_error=False, fill_value=0.0)

@lru_cache(maxsize=8)
def _load_albedo(path):
    # 任意の実装に合わせてください。ここでは (λ[um], albedo) を返す想定。
    lam, alb = get_albedo(path)  # ユーザの関数
    lam = np.asarray(lam, dtype=np.float64)
    alb = np.asarray(alb, dtype=np.float64)
    return lam, alb

@lru_cache(maxsize=8)
def _load_g(path):
    lam, g = get_g_factor(path)  # ユーザの関数
    lam = np.asarray(lam, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64)
    return lam, g

# ------------------------------
#  等角変換: ICRS (RA,Dec) -> Galactic(l,b) を healpy で高速変換
# ------------------------------
_rot_C2G = hp.Rotator(coord=['C', 'G'])  # Equatorial(J2000) → Galactic

def _icrs_to_gal_l_b(ra_deg, dec_deg):
    # healpy の入力は (theta[colat], phi[lon]) [rad]
    theta = np.deg2rad(90.0 - dec_deg)
    phi   = np.deg2rad(ra_deg)
    # (theta, phi) を銀河座標へ回転
    theta_g, phi_g = _rot_C2G(theta, phi)
    l_deg = np.rad2deg(phi_g)
    b_deg = 90.0 - np.rad2deg(theta_g)
    return l_deg, b_deg

def interp_rows_common_x(y_rows, x_src, x_tgt, *, left=0.0, right=0.0):
    """
    y_rows: (Nrow, M)   各行が x_src 上の系列
    x_src:  (M,)        単調増加
    x_tgt:  (K,)        欲しい座標
    戻り値: (Nrow, K)
    """
    y = np.asarray(y_rows, dtype=np.float64)
    xs = np.asarray(x_src,  dtype=np.float64)
    xt = np.asarray(x_tgt,  dtype=np.float64)

    # x_src が昇順でなければ反転
    if xs[0] > xs[-1]:
        xs = xs[::-1]
        y  = y[:, ::-1]

    M = xs.size
    K = xt.size

    # xt を xs で挟むインデックス（右側 i1、左側 i0）
    i1 = np.searchsorted(xs, xt, side='left')
    i0 = i1 - 1

    # 範囲内フラグ（外側は left/right を使う）
    valid = (i1 > 0) & (i1 < M)

    # 端はクリップしておく（値は後で置換）
    i0c = np.clip(i0, 0, M-1)
    i1c = np.clip(i1, 0, M-1)

    x0 = xs[i0c]            # (K,)
    x1 = xs[i1c]            # (K,)
    dx = x1 - x0
    dx[dx == 0.0] = 1.0     # 0割回避（同一点）

    t = (xt - x0) / dx      # (K,)

    # y0, y1 を (Nrow, K) で取得
    y0 = y[:, i0c]          # (Nrow, K)
    y1 = y[:, i1c]          # (Nrow, K)

    out = y0 + (y1 - y0) * t  # (Nrow, K)

    # 外側を left/right で埋める
    if left is not None or right is not None:
        out = out.copy()
        if np.any(~valid):
            maskL = (i1 == 0)
            maskR = (i1 == M)
            if left is not None and np.any(maskL):
                out[:, maskL] = left
            if right is not None and np.any(maskR):
                out[:, maskR] = right
    return out


# ------------------------------
#  メイン: 高速化 DGL
# ------------------------------
def DGL(center_ra, center_dec, params, flist, base_path, makefig=False):
    """
    戻り値:
      dict {
        "Lambda_um": (Nz,),
        "cube": (Nx, Nx, Nz),         # erg/s/cm^2/Hz
        "mask_in": (Nx, Nx),          # True: 円領域内（NaN埋めと同じ位置）
        "ipix": (Npix_in_disc,),      # 使った HEALPix ピクセル
        "pix_grid": (Nx, Nx),         # 各格子点の HEALPix 番号
      }
    """
    # ------------------ 観測の波長レンジ ------------------
    flist = np.asarray(flist, dtype=np.float64)   # Hz, 長さ Nz
    Nz = flist.size
    Lambda = (c / flist) * 1e4                # um
    # 両端±(1/40)帯域分を窓に
    lambda_min = c / flist[-1] * 1e4 - (c / flist[-1]) * 1e4 / 40.0
    lambda_max = c / flist[ 0] * 1e4 + (c / flist[ 0]) * 1e4 / 40.0

    # ------------------ DustEM 出力の読み込み（キャッシュ） ------------------
    wave_SED, I_nu_over_NH = _load_SED(f"{base_path}/out/SED.RES")   # (Nlam_SED,)
    f_ISRF                  = _load_ISRF(f"{base_path}/data/ISRF.DAT") # 関数
    lam_alb, alb_full       = _load_albedo(base_path)                 # (Nlam_alb,)
    lam_g  , g_full         = _load_g(base_path)                      # (Nlam_g,)

    # 波長窓でフィルタ
    sel = (wave_SED > lambda_min) & (wave_SED < lambda_max)
    wl_sel    = wave_SED[sel]                       # (Nlam,)
    I_nu_sel  = I_nu_over_NH[sel]                   # Iν / NH [erg/s/Hz/sr/H]
    ISRF_sel  = f_ISRF(wl_sel)                      # Iν [erg/cm^2/sr/s/Hz]

    # albedo, g を wl_sel に線形補間（高速な np.interp）
    albedo_sel = np.interp(wl_sel, lam_alb, alb_full, left=alb_full[0], right=alb_full[-1])
    g_sel      = np.interp(wl_sel, lam_g,   g_full,   left=g_full[0],   right=g_full[-1])

    # ------------------ NH（円領域内） ------------------
    ebv_map = hp.read_map("/mnt/data_cat3/yuka/data/LIM_mock/csfd_ebv.fits")
    I100 = ebv_map/0.0184               # [MJy/sr] at 100 μm
    nside = hp.npix2nside(I100.size)

    # --- ICRS (RA, Dec) -> Galactic (l, b) ---
    c_icrs = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg, frame='icrs')
    l_deg = c_icrs.galactic.l.value   # degrees
    b_deg = c_icrs.galactic.b.value   # degrees

    #HEALPix の (theta, phi) に変換（theta=co-latitude=90-b, phi=l）
    theta0 = np.deg2rad(90.0 - b_deg)
    phi0   = np.deg2rad(l_deg)
    vec0   = hp.ang2vec(theta0, phi0)

    # 円内ピクセル
    ipix = hp.query_disc(nside, vec0, np.deg2rad(params.radius)*1.3,
                         inclusive=False, fact=4)
    I100_in = I100[ipix].astype(np.float64)          # (Npix,)
    ra_all, dec_all = hp.pix2ang(nside, ipix, lonlat=True)
    all_pix = SkyCoord(ra=ra_all * u.deg, dec=dec_all * u.deg, frame='icrs')
    beta = all_pix.galactic.b.value   # degrees

    
    # MJy/sr → nW m^-2 sr^-1 Hz^-1
    I100_in *= 1e-11
    # ν(100 μm)
    nu100 = c / 100*1e4
    # N_H へ（νIν[100μm] = 18.6e-20 * N_H）
    Nh = (I100_in * nu100) / 18.6*1e20               # (Npix,) [/cm^2]
    Nh[Nh<0] = 0.0
    plt.figure()
    plt.hist(np.log10(Nh), bins=50)
    plt.figure()
    plt.scatter(beta, np.log10(Nh))

    # ------------------ tau(λ) ------------------
    # get_tau は (λ[um], τ[Npix, Nlam]) を返す設計（上の実装と整合）
    wave_tau, tau = get_tau(base_path, Nh)          # wave_tau: (Nlam_tau,), tau: (Npix, Nlam_tau)
    # wl_sel へ補間（ベクトル化）
    # tau: (Npix, Nlam_sel) に
    # np.interp は axis を持たないので、転置→各列独立補間→転置が最速
    tau_pixlam = interp_rows_common_x(tau, wave_tau, wl_sel, left=0.0, right=0.0)  # (Npix, Nlam_sel)

    # ------------------ 散乱 + 熱放射（ベクトル化） ------------------
    alb = albedo_sel[None, :]                         # (1, Nlam)
    gg  = g_sel[None, :]                              # (1, Nlam)
    IS  = ISRF_sel[None, :]                           # (1, Nlam)

    beta  = np.sqrt((1.0 - alb*gg) / (1.0 - alb))     # (1, Nlam)
    alpha = np.sqrt((1.0 - alb) * (1.0 - alb*gg))     # (1, Nlam)

    # exp は (Npix, Nlam) でまとめて
    exp_p = np.exp(alpha * tau_pixlam)
    exp_m = np.exp(-alpha * tau_pixlam)
    denom = (beta + 1.0) * exp_p + (beta - 1.0) * exp_m
    I_scat = 2.0 * beta * IS / denom  - IS*np.exp(-tau_pixlam)                 # (Npix, Nlam)　Iν [erg/cm^2/sr/s/Hz]
    I_scat_1 = alb/(1 - alb) * IS * (1 - np.exp(-(1 - alb)*tau_pixlam))
    I_scat_2 = IS * np.exp(-(1 - alb)*tau_pixlam)*(1 - np.exp(-alb*tau_pixlam))

    I_therm = I_nu_sel[None, :] * Nh[:, None]         # (Npix, Nlam) Iν [erg/cm^2/sr/s/Hz]

    I_tot = I_scat + I_therm                          # (Npix, Nlam) Iν [erg/cm^2/sr/s/Hz]

    # 観測波長グリッド Lambda[um] へ補間（ベクトル化）
    # 転置→各列補間→転置（np.interp は axis を持たないため）
    I_dgl = interp_rows_common_x(np.log10(I_tot), wl_sel, Lambda, left=np.nan, right=np.nan)    # (Npix, Nz)
    I_dgl = 10**I_dgl #[erg/cm^2/sr/s/Hz]

    if makefig:
        # 1) 事前計算
        freq = c / (Lambda * 1e-4)                   # (Nz,)  Hz
        re_100 = interp1d(wave_SED, I_nu_over_NH, bounds_error=False, fill_value=np.nan)
        const_Jy = re_100(100.0) / Jy / 1e6          # [MJy/sr per H]

        # 2) 有効ピクセル抽出（Nh>0, 有限値）
        valid = np.isfinite(Nh) & (Nh > 0)
        I_dgl_valid = I_dgl[valid, :]                # (Npix_valid, Nz)
        Nh_valid = Nh[valid]                          # (Npix_valid,)

        # 3) 全ピクセルの比 R(λ) を一括計算
        #    νIν [erg s^-1 cm^-2 sr^-1] -> [nW m^-2 sr^-1] の 1e6 係数は元コード踏襲
        re_scale_all = I_dgl_valid * 1e6 * freq[None, :]           # (Npix_valid, Nz)
        Jy_100_all   = const_Jy * Nh_valid                          # (Npix_valid,)
        ratio = re_scale_all / Jy_100_all[:, None]                  # (Npix_valid, Nz)

        # 4) パーセンタイル集計（列=λ方向）
        p16, p50, p84 = np.nanpercentile(ratio, [16, 50, 84], axis=0)

        # 5) 図示：中央値 + シェーディング(16–84%)
        plt.figure()
        plt.plot(Lambda, p50, color='k', lw=1.8, label='median')
        plt.fill_between(Lambda, p16, p84, alpha=0.25, label='16–84%')

        # （オプション）範囲を広く見たい場合は 5–95% を追加
        # p5, p95 = np.nanpercentile(ratio, [5, 95], axis=0)
        # plt.fill_between(Lambda, p5, p95, alpha=0.15, label='5–95%')

        plt.xscale('log'); plt.yscale('log')
        plt.xlabel(r'wavelength [$\mu$m]')
        plt.ylabel(r'$\nu_i b_i$ [nW m$^{-2}$ sr$^{-1}$ / MJy sr$^{-1}$]')
        plt.legend()
        plt.tight_layout()
        plt.show()

    # ------------------ 3D キューブへ再配置 ------------------
    Nx = int(2.0 * 3600.0 * params.radius / params.resolution)  # 角秒解像度 → ピクセル数

    ramin  = center_ra  - params.radius
    ramax  = center_ra  + params.radius
    decmin = center_dec - params.radius
    decmax = center_dec + params.radius

    # グリッド（meshgrid はデフォルトで 'xy'）
    ra_lin  = np.linspace(ramin,  ramax,  Nx, dtype=np.float64)
    print(ra_lin)
    dec_lin = np.linspace(decmin, decmax, Nx, dtype=np.float64)
    print(dec_lin)
    ra_grid, dec_grid = np.meshgrid(ra_lin, dec_lin, indexing='xy')  # (Nx, Nx)

    # ICRS → 銀河座標（高速）
    l_deg, b_deg = _icrs_to_gal_l_b(ra_grid.ravel(), dec_grid.ravel())
    pix = hp.ang2pix(nside, l_deg, b_deg, lonlat=True)             # (Nx*Nx,)

    # ipix をソートして searchsorted（ハッシュ辞書より速い）
    order = np.argsort(ipix)
    ipix_sorted = ipix[order]
    idx = np.searchsorted(ipix_sorted, pix)

    # 修正：範囲内だけ比較 → 有効マスクを作る
    valid = (idx < ipix_sorted.size)

    # 比較用のブール配列を用意（範囲外は False のまま）
    match = np.zeros_like(valid, dtype=bool)
    # ※ idx[valid] だけを使って ipix_sorted[...] を参照する
    match[valid] = (ipix_sorted[idx[valid]] == pix[valid])

    # 最終的に「ipix 内に存在する」位置
    valid &= match

    # rows（I_dgl の行番号）を埋める
    rows = np.full(pix.size, -1, dtype=int)
    rows[valid] = order[idx[valid]]

    # 出力配列（外側は NaN）
    cube_flat = np.full((pix.size, Nz), 0.0, dtype=np.float64)
    cube_flat[valid, :] = I_dgl[rows[valid], :] / Jy
    cube = cube_flat.reshape(Nx, Nx, Nz)
    mask_in = valid.reshape(Nx, Nx)

    return cube


    

def old_DGL(center_ra, center_dec, params, flist):
    ebv_map = hp.read_map("/mnt/data_cat3/yuka/data/LIM_mock/csfd_ebv.fits")

    I100 = ebv_map/0.0184 #[MJy/sr]

    nside = hp.get_nside(ebv_map)

    center_frequency = (flist[1::]+flist[:-1:])/2.0
    lambdas = c * 1e4 / (center_frequency) #um
    
    Nx = int(2.0*3600*params.radius / params.resolution)
    
    ramin = center_ra - params.radius
    decmin = center_dec - params.radius
    ramax = center_ra + params.radius
    decmax = center_dec + params.radius
    
    # ベクトルでRA/Dec配列を生成
    ra_grid, dec_grid = np.meshgrid(
        np.linspace(ramin, ramax, Nx),
        np.linspace(decmin, decmax, Nx)
    )
    ra_flat = ra_grid.ravel()
    dec_flat = dec_grid.ravel()

    coords = SkyCoord(ra=ra_flat * u.deg, dec=dec_flat * u.deg, frame='icrs')
    gal_coords = coords.galactic
    l = gal_coords.l.deg
    b = gal_coords.b.deg
    pix = hp.ang2pix(nside, l, b, lonlat=True)

    I100_vals = I100[pix]  # [MJy/sr]
    DGL = np.loadtxt('/mnt/data_cat3/yuka/repository/LIM_mock/DGL.txt', skiprows=1)
    wavelength = DGL[:,0] #um
    mu_bi = DGL[:,1] #[nW m-2 sr-1 / MJy sr-1]

    interp = interp1d(wavelength, mu_bi)
    bi = interp(lambdas)/center_frequency #[nW m-2 sr-1 Hz-1/ MJy sr-1]

    # 波長方向 (Nz,) を (1, 1, Nz) に reshape
    bi = bi[None, None, :]  # shape: (1, 1, Nz)

    # E(B-V) 空間マップ (Nx*Nx,) → (Nx, Nx, 1) に reshape
    I100_vals = I100_vals.reshape(Nx, Nx, 1)  # shape: (Nx, Nx, 1)

    #map
    map = bi * I100_vals * 1e-9  # W m-2 sr-1 Hz-1
    map *= 1e26 #Jy/sr
    return map

######################################################################################################

def interpolate_zodi(result_coarse, Nx):
    Nz = result_coarse.shape[-1]
    x = np.linspace(0, 1, result_coarse.shape[0])
    y = np.linspace(0, 1, result_coarse.shape[1])
    z = np.arange(Nz)
    points = (x, y, z)

    xi = np.linspace(0, 1, Nx)
    yi = np.linspace(0, 1, Nx)
    zi = z
    grid = np.meshgrid(xi, yi, zi, indexing="ij")
    coords = np.stack([g.ravel() for g in grid], axis=-1)

    interp = interpn(points, result_coarse, coords, method="linear", bounds_error=False, fill_value=0)
    result_fine = interp.reshape((Nx, Nx, Nz))
    return result_fine

def evaluate_zodiac_for_day(args):
    day_index, day_string, wavelength, ra_flat, dec_flat, Nx, zodi_model = args
    result = np.zeros((Nx, Nx, len(wavelength)))

    for j, l in enumerate(wavelength):
        model = zodipy.Model(l * u.micron, name=zodi_model)
        skycoords = SkyCoord(ra=ra_flat * u.deg, dec=dec_flat * u.deg,
                             frame="icrs", obstime=day_string)
        intensity = model.evaluate(skycoords).value  # [MJy/sr]
        result[:, :, j] = intensity.reshape((Nx, Nx)) * 1e6  # [Jy/sr]

    return result


def Zodiac(coarse, center_ra, center_dec, params, flist, zodi_model = "dirbe", day_bin=2, count=100, nprocesses=None):
    #ある日時でのimageを出力
    #center_frequency = (flist[1::]+flist[:-1:])/2.0
    #wavelength = c / center_frequency * 1e4  # [um]
    wavelength = c / flist * 1e4
 
    Nx = int(2.0*3600*params.radius / params.resolution/coarse)
    
    ramin = center_ra - params.radius
    decmin = center_dec - params.radius
    ramax = center_ra + params.radius
    decmax = center_dec + params.radius
    
    ra_grid, dec_grid = np.meshgrid(
        np.linspace(ramin, ramax, Nx),
        np.linspace(decmin, decmax, Nx)
    )
    ra_flat = ra_grid.ravel()
    dec_flat = dec_grid.ravel()

    init_day = datetime.date(2025, 3, 22)
    args_list = []
    args_list_1 = []
    for i in range(count):
        day = init_day + datetime.timedelta(days=day_bin * i)
        day_string = day.strftime('%Y-%m-%d')
        args_list.append((i, day_string, wavelength, ra_flat, dec_flat, Nx, zodi_model))

    init_day = day
    for i in range(count):
        day = init_day + datetime.timedelta(days=day_bin * i)
        day_string = day.strftime('%Y-%m-%d')
        args_list_1.append((i, day_string, wavelength, ra_flat, dec_flat, Nx, zodi_model))

    # 並列化実行
    with multiprocessing.Pool(10) as pool:
        results = pool.map(evaluate_zodiac_for_day, args_list)

    with multiprocessing.Pool(10) as pool:
        results1 = pool.map(evaluate_zodiac_for_day, args_list_1)

    # 結果をまとめて4次元配列に
    full_spectra = np.mean(results, axis=0)  # shape:  (Nx/32, Nx/32, Nz)
    full_spectra_1 = np.mean(results1, axis=0)  # shape:  (Nx/32, Nx/32, Nz)
    
    full_Nx = int(2.0*3600*params.radius / params.resolution)
    zodi_final = interpolate_zodi(full_spectra, full_Nx)
    zodi_final_1 = interpolate_zodi(full_spectra_1, full_Nx)
    return zodi_final, zodi_final_1

#################################################################################################

def Mask(center_ra, center_dec, params, mask_ra, mask_dec):
    Nx = int(2.0*3600*params.radius / params.resolution)

    ramin = center_ra - params.radius
    decmin = center_dec - params.radius
    
    #観測領域外を0に
    ix = np.arange(Nx)
    iy = np.arange(Nx)
    xx, yy = np.meshgrid(ix, iy)
    center_i = 3600 * params.radius / params.resolution
    radius_i = 2 * 3600 * params.radius / params.resolution
    xx = xx - center_i
    yy = yy - center_i
    mask = ((xx**2)+(yy**2))<radius_i**2
    
    star_mask = np.ones([Nx, Nx])
    # Mag ≤ 19の星のピクセルもマスクする
    for ra_star, dec_star in zip(mask_ra, mask_dec):
        x = int((ra_star - ramin) * 3600.0 / params.resolution)
        y = int((dec_star - decmin) * 3600.0 / params.resolution)
        if 0 <= x < Nx and 0 <= y < Nx:
            # 半径Nピクセル（例: 2）以内をマスク
            radius = 2
            x_min, x_max = max(0, x - radius), min(Nx, x + radius + 1)
            y_min, y_max = max(0, y - radius), min(Nx, y + radius + 1)
            star_mask[y_min:y_max, x_min:x_max] = 0.0
    
    return mask, star_mask

def save_meta_parquet(outdir, ra, dec, ix, iy, metal, T_eff, logg, flist):
    os.makedirs(outdir, exist_ok=True)

    # ここで全てネイティブ化（dtype も明示）
    ra_n  = to_native(ra, dtype=np.float32)
    dec_n  = to_native(dec,  dtype=np.float32)
    ix_n    = to_native(ix,    dtype=np.int64)
    iy_n    = to_native(iy,    dtype=np.int64)
    Z_n     = to_native(metal, dtype=np.float32)   # ★ ここが >f4 で落ちていた
    Teff_n  = to_native(T_eff, dtype=np.float32)
    logg_n  = to_native(logg,  dtype=np.float32)
    nu_min  = float(np.min(flist))
    nu_max  = float(np.max(flist))

    df = pd.DataFrame({
        "star_id": np.arange(ix_n.size, dtype=np.int64),
        "ra": ra_n,
        "dec": dec_n,
        "x": ix_n,
        "y": iy_n,
        "metallicity": Z_n,
        "Teff": Teff_n,
        "logg": logg_n,
        "nu_min_Hz": nu_min,
        "nu_max_Hz": nu_max,
    })

    # 念のため DataFrame 側でも統一
    df = df.astype({
        "star_id": "int64",
        "ra": "float32",
        "dec": "float32",
        "x": "int64",
        "y": "int64",
        "metallicity": "float32",
        "Teff": "float32",
        "logg": "float32",
        "nu_min_Hz": "float64",
        "nu_max_Hz": "float64",
    })

    df.to_parquet(os.path.join(outdir, "stars_meta_all.parquet"), index=False)


def save_spectra_h5(outdir, flist, dflist, spectrum):
    path = os.path.join(outdir, "stellar_spectra_all.h5")
    with h5py.File(path, "w") as h5:
        h5.create_dataset("flist_Hz",  data=flist,  compression="gzip")
        h5.create_dataset("dflist_Hz", data=dflist, compression="gzip")
        h5.create_dataset(
            "spectrum", data=spectrum.astype(np.float32),
            compression="gzip", compression_opts=4,
            chunks=(min(1024, spectrum.shape[0]), min(2048, spectrum.shape[1]))
        )
        
def to_native(a, dtype=None):
    """配列 a をネイティブエンディアン & 連続メモリに揃える"""
    x = np.asarray(a, dtype=dtype if dtype is not None else a.dtype)
    # '|' はバイト順序なし（例：bool, bytes）、'=' はネイティブ
    if hasattr(x, 'dtype') and x.dtype.byteorder not in ('=', '|'):
        x = x.byteswap().newbyteorder()      # 実データの並びを入れ替えてネイティブへ
    return np.ascontiguousarray(x)           # ついでに C 連続に
        