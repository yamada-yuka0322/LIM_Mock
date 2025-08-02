import numpy as np
import os

import pandas as pd

from scipy.interpolate import interp1d
from scipy.interpolate import interpn

from astropy.cosmology import Planck18 as cosmo
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic
from astropy.coordinates import FK5
from astropy.io import fits
from astropy.table import Table

import healpy as hp

import matplotlib.pyplot as plt

import zodipy

import datetime

import multiprocessing

Oiii = 5.997e5 #[GHz]
ha = 4.856e5 #[GHz]
pa = 1.599e5 #[GHz]
c = 2.9979e10 #[cm/s]
Jy = 1.0e-23
arcsec = 4.848136811094e-6
pc = 3085677857083127000

class Params(object):
    def __init__(self):
        self.frequency = None
        self.R = None
        self.resolution = None
        self.radius = None
    def Default(self):
        self.frequency = [120000, 400000] #[GHz, GHz]
        self.R = 41
        self.resolution = 6.5 #[arcsec]
        self.radius = 7.5 #[deg]

def frequency(params):
    flist = []
    dflist = []
    
    fmin, fmax = params.frequency #[GHz, GHz]
    fnow = fmin * 1e9 #[Hz]
    while fnow <= fmax * 1e9:
        flist.append(fnow)
        df = fnow/params.R
        dflist.append(df)
        
        fnow += df
        
    flist = np.array(flist, dtype=np.float32)
    dflist = np.array(dflist, dtype=np.float32)
    return flist, dflist

def make_mock(path, params):
    flist, dflist = frequency(params) # Hz, Hz
    
    df = pd.read_csv(path)
    lumi_dist = cosmo.luminosity_distance(np.array(df['truez'])).to(u.cm)
    lumi_dist = lumi_dist.value
    
    ra = np.array(df['ra']) * 3600.0 #[arcsec]
    dec = np.array(df['dec']) * 3600.0 #[arcsec]
    
    ix = np.floor((ra + params.radius * 3600) / params.resolution).astype(np.int64)
    iy = np.floor((dec + params.radius * 3600) / params.resolution).astype(np.int64)
    print(ix)
    
    ##initialize
    Nx = int(2.0*3600*params.radius / params.resolution)
    Nz = len(flist) - 1
    
    npix = np.array([Nx, Nx, Nz])
    ########################### Ha intensity
    Ha_intensity = np.zeros([Nx, Nx, Nz], dtype=np.float32)
    
    freq_obs = ha/(1 + np.array(df['obsz']))*1e9 #[Hz]
    iz = np.searchsorted(flist, freq_obs, side='right') - 1
    
    indices = np.array([ix, iy, iz]).T
    valid_mask = np.all((indices >=0)&(indices < npix), axis=1)
    
    indices_valid = indices[valid_mask]
    lumi_dist_valid = lumi_dist[valid_mask]
    lumi_valid = np.array(df['Ha'])[valid_mask] #[erg/s]
    flux_valid = lumi_valid / (4. * np.pi * lumi_dist_valid**2) #[erg/s/cm^2]
    intensity_valid = flux_valid / dflist[indices_valid[:,2]] /Jy / (params.resolution * arcsec)**2 #[Jy/sr]
    
    np.add.at(Ha_intensity, (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]), intensity_valid)
    ########################### Pa intensity
    Pa_intensity = np.zeros([Nx, Nx, Nz], dtype=np.float32)
    
    freq_obs = pa/(1 + np.array(df['obsz']))*1e9 #[Hz]
    iz = np.searchsorted(flist, freq_obs, side='right') - 1
    
    indices = np.array([ix, iy, iz]).T
    valid_mask = np.all((indices >=0)&(indices < npix), axis=1)
    
    indices_valid = indices[valid_mask]
    lumi_dist_valid = lumi_dist[valid_mask]
    lumi_valid = np.array(df['Pa'])[valid_mask] #[erg/s]
    flux_valid = lumi_valid / (4. * np.pi * lumi_dist_valid**2) #[erg/s/cm^2]
    intensity_valid = flux_valid / dflist[indices_valid[:,2]] /Jy / (params.resolution * arcsec)**2 #[Jy/sr]
    
    np.add.at(Pa_intensity, (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]), intensity_valid)
    
    ########################### [OIII] intensity
    OIII_intensity = np.zeros([Nx, Nx, Nz], dtype=np.float32)
    
    freq_obs = Oiii/(1 + np.array(df['obsz'])) * 1e9 #[Hz]
    iz = np.searchsorted(flist, freq_obs, side='right') - 1
    
    indices = np.array([ix, iy, iz]).T
    valid_mask = np.all((indices >=0)&(indices < npix), axis=1)
    
    indices_valid = indices[valid_mask]
    lumi_dist_valid = lumi_dist[valid_mask]
    lumi_valid = np.array(df['OIII'])[valid_mask] #[erg/s]
    flux_valid = lumi_valid / (4. * np.pi * lumi_dist_valid**2) #[erg/s/cm^2]
    intensity_valid = flux_valid / dflist[indices_valid[:,2]] /Jy / (params.resolution * arcsec)**2 #[Jy/sr]
    
    np.add.at(OIII_intensity, (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]), intensity_valid)
    
    total_intensity = Ha_intensity + OIII_intensity + Pa_intensity
    
    return total_intensity, Ha_intensity, OIII_intensity, Pa_intensity

####################################################################################################
def make_foreground(ra, dec, params):
    flist, dflist = frequency(params) # Hz, Hz
    star_path = '/mnt/data_cat3/yuka/data/LIM_mock/stars/output.fits'
    star, mask_ra, mask_dec = Star(ra, dec, params, flist, star_path)
    dgl = DGL(ra, dec, params, flist)
    mask, star_mask = Mask(ra, dec, params, mask_ra, mask_dec)
    zodiac = Zodiac(32, ra, dec, params, flist)
    
    total = (star + dgl + zodiac) * (mask*star_mask)[:, :, np.newaxis]
    
    return star, dgl, zodiac, mask, star_mask, total

#####################################################
def Star(center_ra, center_dec, params, flist, path):
    Nx = int(2.0*3600*params.radius / params.resolution)
    Nz = len(flist) - 1
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
    area = (distance < 7.5)
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
    
    for i in data[stars]:
        _ra = i['RAJ2000']
        _dec = i['DECJ2000']
        _J = i['J']
        
        i_ra = (_ra - ra_min) / params.resolution *3600
        i_ra = i_ra.astype(np.int32)
        i_dec = (_dec - dec_min) / params.resolution *3600
        i_dec = i_dec.astype(np.int32)
        
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
            
            star_intensity[i_dec, i_ra] += F
            if(_J<16):
                bright_star_ra.append(_ra)
                bright_star_dec.append(_dec)
                bright_star_ira.append(i_ra)
                bright_star_idec.append(i_dec)
                bright_star_path.append(filename)
                bright_star_logg.append(column)
    
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
    closest = int(_T[np.argmin(_T - Teff)**2])
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
    
    f_nu_list = []
    
    for i in range(len(frequency) - 1):
        nu1 = frequency[i]
        nu2 = frequency[i+1]
        
        nu_grid = np.linspace(nu1, nu2, 100)
        f_vals = f_nu_interp(nu_grid)
        
        f_bin_avg = np.trapz(f_vals, nu_grid) / (nu2 - nu1) #[erg/s/cm^2/Hz]
        f_Jy = f_bin_avg /Jy
        f_nu_list.append(f_Jy)
        
    return np.array(f_nu_list)

############################################################################################################

def DGL(center_ra, center_dec, params, flist):
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
    day_index, day_string, wavelength, ra_flat, dec_flat, Nx = args
    result = np.zeros((Nx, Nx, len(wavelength)))

    for j, l in enumerate(wavelength):
        model = zodipy.Model(l * u.micron)
        skycoords = SkyCoord(ra=ra_flat * u.deg, dec=dec_flat * u.deg,
                             frame="icrs", obstime=day_string)
        intensity = model.evaluate(skycoords).value  # [MJy/sr]
        result[:, :, j] = intensity.reshape((Nx, Nx)) * 1e6  # [Jy/sr]

    return result


def Zodiac(coarse, center_ra, center_dec, params, flist, day_bin=2, count=100, nprocesses=None):

    center_frequency = (flist[1::]+flist[:-1:])/2.0
    wavelength = c / center_frequency * 1e4  # [um]
 
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
    for i in range(count):
        day = init_day + datetime.timedelta(days=day_bin * i)
        day_string = day.strftime('%Y-%m-%d')
        args_list.append((i, day_string, wavelength, ra_flat, dec_flat, Nx))

    # 並列化実行
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        results = pool.map(evaluate_zodiac_for_day, args_list)

    # 結果をまとめて4次元配列に
    full_spectra = np.mean(results, axis=0)  # shape:  (Nx/32, Nx/32, Nz)
    
    full_Nx = int(2.0*3600*params.radius / params.resolution)
    zodi_final = interpolate_zodi(full_spectra, full_Nx)
    return zodi_final

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
    radius_i = 3600 * params.radius / params.resolution
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
        