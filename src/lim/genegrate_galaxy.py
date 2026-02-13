import astropy.io.fits as fits
import sys
import numpy as np
from astropy.cosmology import Planck18 as cosmo
from astropy import units as u
from astropy.cosmology import z_at_value

from multiprocessing  import Pool

import matplotlib.pyplot as plt

#from speculator import *

from functools import partial

import os
os.environ['SPS_HOME'] = '/mnt/data_cat3/yuka/repository/fsps'

import fsps
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid

import time

from lim.create_JWST import frequency
from lim.AbundanceMatch import GetEL

import h5py
import pandas as pd

log_Lsun = 33.58297 # log10 of solar luminosity in erg/s/Hz

# グローバル変数（ワーカーごとに中身が別になる）
prospector_alpha_UV = None
prospector_alpha_opticalNIR = None
prospector_alpha_IR = None

sp = None
wavelength = None
filters = None
#log_Lsun = np.log10(L_sun.cgs.value)

Oiii = 5.997e5 #[GHz]
ha = 4.856e5 #[GHz]
pa = 1.599e5 #[GHz]
c = 2.9979e10 #[cm/s]

def init_worker_fsps():
    """各プロセス起動時に一度だけ呼ばれる。"""
    global sp
    sp = fsps.StellarPopulation(
        zcontinuous=1,
        add_neb_emission=False,
        add_neb_continuum=True,
        add_dust_emission=True,
        redshift_colors=1,
        sfh=3,
        imf_type=1,
        dust_type=0,
    )
    
def Get_Sides(name):
    file = '/mnt/data_cat3/yuka/data/SIDES/' + f'pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits'
    hdu = fits.open(file)
    data = hdu[1].data
    ra = data['ra']
    dec = data['dec']
    redshift = data['redshift']
    Mhalo = np.log10(data['Mhalo']) #[Msun]
    Mstar = data['Mstar'] #[Msun]
    
    selection = redshift < 4.0

    return ra[selection], dec[selection], redshift[selection], Mhalo[selection], Mstar[selection]

def get_indices_with_fallback(key, bin_to_indices, max_zbin, max_mbin, max_radius=3):
    zb, mb = key
    # まずはそのまま
    idxs = bin_to_indices.get(key, None)
    if idxs is not None and len(idxs) > 0:
        return idxs

    # 周囲をだんだん広げながら探す
    for r in range(1, max_radius+1):
        cand_keys = []
        for dz in range(-r, r+1):
            for dm in range(-r, r+1):
                nz, nm = zb + dz, mb + dm
                if 0 <= nz < max_zbin and 0 <= nm < max_mbin:
                    cand_keys.append((nz, nm))

        for k in cand_keys:
            idxs = bin_to_indices.get(k, None)
            if idxs is not None and len(idxs) > 0:
                return idxs

    # どうしても見つからなければ None
    return None


def sample_params(redshift, logM, samples, bin_to_indices, z_edges, logM_edges):

    z_bins = np.digitize(redshift, z_edges) - 1
    logM_bins = np.digitize(logM, logM_edges) - 1

    z_bins = np.clip(z_bins, 0, len(z_edges)-1)
    logM_bins = np.clip(logM_bins, 0, len(logM_edges)-1)

    B = len(redshift)
    max_zbin = len(z_edges) - 1
    max_mbin = len(logM_edges) - 1

    sample_list = []
    for i in range(B):
        key = (z_bins[i], logM_bins[i])
        #idxs = bin_to_indices.get(key, None)
        idxs = get_indices_with_fallback(key, bin_to_indices, max_zbin, max_mbin)

        if idxs is None or len(idxs) == 0:
            print(f"No data for key = {key}")
            sys.exit(1)
        else:
            j = np.random.randint(len(idxs))
            chosen_idx = idxs[j]
            sample_list.append( samples[chosen_idx] )

    params = np.stack(sample_list, axis=0) # (B, num_params)

    return params

def sample_SEDs(redshift, logM, SED, SFR, bin_to_indices, z_edges, logM_edges, f115w, f150w, f277w, f444w):

    z_bins = np.digitize(redshift, z_edges) - 1
    logM_bins = np.digitize(logM, logM_edges) - 1

    z_bins = np.clip(z_bins, 0, len(z_edges)-1)
    logM_bins = np.clip(logM_bins, 0, len(logM_edges)-1)

    B = len(redshift)
    max_zbin = len(z_edges) - 1
    max_mbin = len(logM_edges) - 1

    sample_list = []
    sfr_list = []
    f115w_list = []
    f150w_list = []
    f277w_list = []
    f444w_list = []
    idx_list = []
    
    for i in range(B):
        key = (z_bins[i], logM_bins[i])
        #idxs = bin_to_indices.get(key, None)
        idxs = get_indices_with_fallback(key, bin_to_indices, max_zbin, max_mbin)

        if idxs is None or len(idxs) == 0:
            print(f"No data for key = {key}")
            sys.exit(1)
        else:
            j = np.random.randint(len(idxs))
            chosen_idx = idxs[j]
            sample_list.append( SED[chosen_idx] )
            sfr_list.append( SFR[chosen_idx] )
            f115w_list.append(f115w[chosen_idx])
            f150w_list.append(f150w[chosen_idx])
            f277w_list.append(f277w[chosen_idx])
            f444w_list.append(f444w[chosen_idx])
            idx_list.append(chosen_idx)

    spectrum = np.stack(sample_list, axis=0) # (B, num_params)
    sfr = np.stack(sfr_list, axis=0)
    F115W = np.array(f115w_list)
    F150W = np.array(f150w_list)
    F277W = np.array(f277w_list)
    F444W = np.array(f444w_list)
    idx = np.array(idx_list)

    return spectrum, sfr, F115W, F150W, F277W, F444W, idx

def sample_ids(redshift, logM, bin_to_indices, z_edges, logM_edges):
    z_bins = np.digitize(redshift, z_edges) - 1
    logM_bins = np.digitize(logM, logM_edges) - 1

    z_bins = np.clip(z_bins, 0, len(z_edges)-1)
    logM_bins = np.clip(logM_bins, 0, len(logM_edges)-1)

    B = len(redshift)
    max_zbin = len(z_edges) - 1
    max_mbin = len(logM_edges) - 1

    idx_list = []
    
    for i in range(B):
        key = (z_bins[i], logM_bins[i])
        #idxs = bin_to_indices.get(key, None)
        idxs = get_indices_with_fallback(key, bin_to_indices, max_zbin, max_mbin)

        if idxs is None or len(idxs) == 0:
            print(f"No data for key = {key}")
            sys.exit(1)
        else:
            j = np.random.randint(len(idxs))
            chosen_idx = idxs[j]
            idx_list.append(chosen_idx)

    idx = np.array(idx_list)

    return idx

def get_true_spectrum(args):
    global sp
    theta_samples, logM_formed, mw_age, M_frac_bins, t_edge_bins = args
    N_pred, zmeta_pred, sfr1_pred, sfr2_pred, sfr3_pred, sfr4_pred, sfr5_pred, sfr6_pred,  tau2_pred, dust_n_pred, tau1_pred, fAGN_pred, tauAGN_pred, gasZ_pred, gasU_pred, redshift_pred = theta_samples
    SFH_time = t_edge_bins[-1] - t_edge_bins
    SFR = 10**M_frac_bins

    tage = mw_age
    SFH_time = SFH_time[::-1]
    SFR = SFR[::-1]

    sp.params['zred']      = redshift_pred
    sp.params['logzsol']   = zmeta_pred
    sp.params['gas_logz']  = gasZ_pred
    sp.params['gas_logu']  = gasU_pred
    sp.params['dust1']     = tau2_pred * tau1_pred
    sp.params['dust2']     = tau2_pred
    sp.params['dust_index']= dust_n_pred
    sp.params['fagn']      = 10**fAGN_pred
    sp.params['agn_tau']   = 10**tauAGN_pred
    sp.params['tage']      = tage
    age_seq = (SFH_time[1::]+SFH_time[:-1:])*0.5
    sp.set_tabular_sfh(age_seq, SFR)

    wave, spec = sp.get_spectrum(tage=tage, peraa=False)  # spec in Lsun/Hz
    log_spec = np.log10(spec) + log_Lsun  # in erg/s/Hz
    f115w, f150w, f277w, f444w = sp.get_mags(tage=tage, redshift=redshift_pred, bands=['jwst_f115w', 'jwst_f150w', 'jwst_f277w', 'jwst_f444w'])

    return wave, log_spec, f115w, f150w, f277w, f444w

def _bandavg_group_top_hat(wavelength, log_Lnu_SED, nu0_obs, dnu_obs, z):
        """
        テンプレ (nu_grid, Lnu_grid) を、各銀河の z に対して
        観測チャンネル (nu0_obs, dnu_obs) を top-hat として帯域平均 Lν を返す。
        返り値 shape = (N_gal, N_chan)
        wavelength: Å
        Lnu_SED: erg/s/Hz
        nu0_obs: Hz
        dnu_obs: Hz
        z: redshift
        -----------
        Lnu_bandavg: erg/s/cm2/Hz
        -----------
        """
        # 1) 累積積分 I(ν) を作成（テンプレ1本につき一回）
        # cumulative_trapezoid(x,y): ∫ y dx を返す（len=N-1）。先頭0を付けて長さを揃える
        nu_grid = c / (wavelength * 1e-8)   # Hz
        
        lumi_dist = cosmo.luminosity_distance(np.array(z)).to(u.cm)
        lumi_dist = lumi_dist.value
        log_geom = 2* np.log10(lumi_dist)[:, None] + np.log10(4.0 * np.pi)
        log_Lnu_SED = log_Lnu_SED - log_geom   # erg/s/cm2/Hz
        Lnu_SED = 10**log_Lnu_SED          # erg/s/cm2/Hz
        I = cumulative_trapezoid(Lnu_SED, nu_grid, initial=0.0)  # shape = (N_grid,)

        # 2) 各銀河の rest-frame のバンド端を作る
        #    ν_rest = (1+z) * ν_obs
        nu1 = (1.0 + z)[:, None] * (nu0_obs[None, :] - 0.5 * dnu_obs[None, :])   # (N_gal, N_chan)
        nu2 = (1.0 + z)[:, None] * (nu0_obs[None, :] + 0.5 * dnu_obs[None, :])

        # 3) I(ν) を ν1, ν2 に補間（範囲外は端値）
        I1 = np.interp(nu1.ravel(), nu_grid, I, left=I[0], right=I[-1]).reshape(nu1.shape)
        I2 = np.interp(nu2.ravel(), nu_grid, I, left=I[0], right=I[-1]).reshape(nu2.shape)

        # 4) 帯域平均 Lν
        dnu_rest = np.clip(nu2 - nu1, 1e-30, None)
        Lnu_bandavg = (I2 - I1) / dnu_rest     # (N_gal, N_chan)
        return Lnu_bandavg
    
def bandavg_group_top_hat(
    wavelength_A,
    log_Lnu,          # (Nobj, Nwave) or (Nwave,)  [erg/s/Hz]
    nu0_obs_Hz,       # (Nchan,) [Hz]
    dnu_obs_Hz,       # (Nchan,) [Hz]
    z,                # (Nobj,) or scalar
    chunk_size=2048,  # adjust for memory
    dtype=np.float64
):
    """
    Top-hat band-averaged observed-frame flux density Fnu_obs [erg/s/cm^2/Hz].

    Inputs
    ------
    wavelength_A : (Nwave,) Å (assumed monotonic, usually increasing)
    log_Lnu      : log10(Lnu) [erg/s/Hz], shape (Nobj, Nwave) or (Nwave,)
    nu0_obs_Hz   : channel centers [Hz], shape (Nchan,)
    dnu_obs_Hz   : channel widths  [Hz], shape (Nchan,)
    z            : redshift per object, shape (Nobj,) or scalar

    Returns
    -------
    Fnu_bandavg : (Nobj, Nchan) [erg/s/cm^2/Hz]
                  (or (Nchan,) if input log_Lnu was 1D and z scalar)
    """

    # --- sanitize shapes ---
    wl = np.asarray(wavelength_A, dtype=dtype)
    nu0 = np.asarray(nu0_obs_Hz, dtype=dtype)
    dnu = np.asarray(dnu_obs_Hz, dtype=dtype)

    logL = np.asarray(log_Lnu, dtype=dtype)
    one_d = (logL.ndim == 1)
    if one_d:
        logL = logL[None, :]  # (1, Nwave)

    z = np.asarray(z, dtype=dtype)
    if z.ndim == 0:
        z = np.full((logL.shape[0],), float(z), dtype=dtype)
    elif z.shape[0] != logL.shape[0]:
        raise ValueError(f"z length {z.shape[0]} must match log_Lnu first dim {logL.shape[0]}")

    Nobj, Nwave = logL.shape
    Nchan = nu0.shape[0]

    # --- build frequency grid, ensure increasing order ---
    nu_grid = (c / (wl * 1e-8)).astype(dtype)  # Hz
    order = np.argsort(nu_grid)                # ensure ascending
    nu_grid = nu_grid[order]
    logL = logL[:, order]

    # Precompute distances and geometric dimming
    dL_cm = cosmo.luminosity_distance(z).to(u.cm).value.astype(dtype)  # (Nobj,)
    #geom = (4.0 * np.pi * dL_cm**2).astype(dtype)                      # (Nobj,)
    geom = np.log10(4.0 * np.pi) + 2.0 * np.log10(dL_cm)

    # observed top-hat edges (obs frame)
    nu1_obs = (nu0 - 0.5*dnu)[None, :]  # (1, Nchan)
    nu2_obs = (nu0 + 0.5*dnu)[None, :]

    # rest-frame edges for each object
    # nu_rest = (1+z)*nu_obs
    zp1 = (1.0 + z)[:, None]  # (Nobj,1)
    nu1 = zp1 * nu1_obs       # (Nobj, Nchan)
    nu2 = zp1 * nu2_obs       # (Nobj, Nchan)

    # rest bandwidth
    dnu_rest = (nu2 - nu1)
    dnu_rest = np.clip(dnu_rest, 1e-60, None)

    # Output
    out = np.empty((Nobj, Nchan), dtype=dtype)

    # Helper: vectorized linear interp of y(x) with y shape (B, Nwave), xq shape (B, Nchan)
    def interp_2d(x_grid, y_grid, xq):
        """
        x_grid: (Nwave,) ascending
        y_grid: (B, Nwave)
        xq    : (B, Nchan)
        returns yq: (B, Nchan)
        """
        Nw = x_grid.shape[0]
        # clip query to grid range (constant extrapolation)
        xq_clip = np.clip(xq, x_grid[0], x_grid[-1])

        hi = np.searchsorted(x_grid, xq_clip, side='right')
        hi = np.clip(hi, 1, Nw-1)
        lo = hi - 1

        x_lo = x_grid[lo]                 # (B, Nchan) via broadcasting? nope, lo is (B,Nchan)
        x_hi = x_grid[hi]
        # gather y at lo/hi
        y_lo = np.take_along_axis(y_grid, lo, axis=1)
        y_hi = np.take_along_axis(y_grid, hi, axis=1)

        t = (xq_clip - x_lo) / (x_hi - x_lo)
        return y_lo + (y_hi - y_lo) * t

    # --- chunk loop to control memory ---
    for s in range(0, Nobj, chunk_size):
        e = min(Nobj, s + chunk_size)

        # Lnu [erg/s/Hz]
        Lnu = np.power(10.0, logL[s:e, :] - geom[s:e, None], dtype=dtype)

        # cumulative integral I(nu) = ∫ Lnu dnu along nu_grid
        # shape: (B, Nwave)
        I = cumulative_trapezoid(Lnu, nu_grid, axis=1, initial=0.0)

        # interpolate I at rest-frame edges
        I1 = interp_2d(nu_grid, I, nu1[s:e, :])
        I2 = interp_2d(nu_grid, I, nu2[s:e, :])

        # band-averaged rest-frame <Lnu>_rest
        Lnu_avg_rest = (I2 - I1) / dnu_rest[s:e, :]

        # Observed-frame Fnu:
        # Fnu_obs = (1+z)^(-1) * <Lnu>_rest / (4π dL^2)
        out[s:e, :] = (Lnu_avg_rest) * (1.0 + z[s:e, None])

    if one_d:
        # if original was 1D and z was scalar, return (Nchan,)
        if out.shape[0] == 1:
            return out[0]
    return out

def plot_LF(mag, name):
    bins_m = np.arange(15, 31, 0.5)
    hist_m, _ = np.histogram(mag, bins=bins_m)
    dNdm = hist_m / np.diff(bins_m) * 10

    m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
    err_dNdm = np.sqrt(hist_m) / np.diff(bins_m)

    data = np.loadtxt(f'../{name}.txt', delimiter=',')

    plt.figure()
    plt.errorbar(m_cent, dNdm, yerr=err_dNdm, fmt='o', ms=3, label='mock (cont+lines)')
    plt.scatter(data[:,0], data[:,1], label='observation', color='red', s=10)
    plt.yscale('log'); plt.xlabel(f'm_AB ({name})'); plt.ylabel('dN/dm [deg$^{-2}$ mag$^{-1}$]')
    plt.legend(); plt.tight_layout()
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_SIDES_{name}.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_SIDES_{name}.pdf")
    
def generate_catalog(names, params):
    init_worker_fsps()
    
    flist, dflist = frequency(params)
    
    with Pool(processes=20) as pool:
        results = pool.map(Get_Sides, names)
    ra, dec, photoz, Mhalo, Mstar = zip(*results)
    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    photoz = np.concatenate(photoz)
    Mhalo = np.concatenate(Mhalo)
    Mstar = np.concatenate(Mstar)
    
    f_sample = "/mnt/data_cat3/yuka/repository/LIM_mock/theta_samples_100000_seed1_mag.npz"
    data = np.load(f_sample, allow_pickle=True)
    samples = data["samples"]
    logM_samples = data["logM_samples"]
    logM_formed = data['logM_formed']
    tage = data['tage']
    #sfr = data['SFR']
    log10SFR_bins = data['log10SFR_bins']
    t_edge_bins = data['t_edge_bins']
    redshift = samples[:,-1]

    print("Load {}".format(f_sample))

    worker_iter = zip(samples, logM_formed, tage, log10SFR_bins, t_edge_bins)
    with Pool(processes=20) as pool:
        results = pool.map(get_true_spectrum, worker_iter)
        #results = pool.map(get_true_spectrum, samples)
    waves, log_abs_spectrum, f115w, f150w, f277w, f444w = zip(*results)
    log_abs_spectrum = np.array(log_abs_spectrum) # (N_samples, N_lambda) rest frame luminosity[erg/s/Hz]
    f115w = np.array(f115w)
    f150w = np.array(f150w)
    f277w = np.array(f277w)
    f444w = np.array(f444w)
    
    global wavelength
    wavelength = np.array(waves[0])
    
    # create (z, logM) bins
    Nbin = 20
    z_edges = np.linspace(0, 4, Nbin)
    logM_edges = np.linspace(7, 11, Nbin)

    z_idx = np.digitize(redshift, z_edges) - 1
    logM_idx = np.digitize(logM_samples, logM_edges) - 1

    z_idx = np.clip(z_idx, 0, Nbin - 1)
    logM_idx = np.clip(logM_idx, 0, Nbin - 1)

    from collections import defaultdict
    bin_to_indices = defaultdict(list)

    for i in range(len(samples)):
        key = (z_idx[i], logM_idx[i])
        bin_to_indices[key].append(i)

    for key in bin_to_indices:
        bin_to_indices[key] = np.array(bin_to_indices[key])
    
    z_array = photoz
    M_array = np.log10(Mstar)
    SFR = 10**log10SFR_bins[:,0]
    
    obs_spectrum = bandavg_group_top_hat(wavelength, log_abs_spectrum, flist, dflist, redshift)
    nu_band = c / (np.array([1.15, 1.50, 2.77, 4.44])/1e4) #Hz
    plt.figure()
    colorlist = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd']
    for i in range(5):
        index = i*1000
        band_flux = abmag_to_fnu(np.array([f115w[index], f150w[index], f277w[index], f444w[index]]))
        spectrum = obs_spectrum[index]
        plt.plot(flist/1e9, spectrum, color=colorlist[i])
        plt.scatter(nu_band/1e9, band_flux, color=colorlist[i])
    plt.xlim(60000, 400000)
    plt.xscale('log')
    plt.yscale('log')
    log_spectrum, sfr, f115w_pred, f150w_pred, f277w_pred, f444w_pred, idx = sample_SEDs(z_array, M_array, obs_spectrum, SFR, bin_to_indices, z_edges, logM_edges, f115w, f150w, f277w, f444w)
    
    plt.figure()
    for i in range(5):
        index = i * 100
        sample_index = idx[index]
        nu_band = c / (np.array([1.15, 1.50, 2.77, 4.44])/1e4) #Hz
        band_flux = abmag_to_fnu(np.array([f115w_pred[index], f150w_pred[index], f277w_pred[index], f444w_pred[index]]))
        spectrum = log_spectrum[index]
        
        obs_band_flux = abmag_to_fnu(np.array([f115w[sample_index], f150w[sample_index], f277w[sample_index], f444w[sample_index]]))
        obs_spec = obs_spectrum[sample_index]
        #plt.plot(flist/1e9, spectrum, color=colorlist[i])
        plt.plot(flist/1e9, obs_spec, color=colorlist[i], ls="--")
        #plt.scatter(nu_band/1e9, band_flux, color=colorlist[i])
        plt.scatter(nu_band/1e9, obs_band_flux, marker="^",color=colorlist[i])
    plt.xlim(60000, 400000)
    plt.xscale('log')
    plt.yscale('log')
    
    Ha = GetEL(sfr, z_array, area=1.0, name='Halpha') #erg/s
    OIII = GetEL(sfr, z_array, area=1.0, name='OIII') #erg/s
    log_Ha = np.log10(Ha)
    log_Pa = log_Ha - 0.60
    Pa = 10**log_Pa #erg/s
    
    geom = 4.0 * np.pi * (cosmo.luminosity_distance(np.array(z_array)).to(u.cm).value)**2  # cm2
    
    fedges = np.empty(len(flist) + 1, dtype=float)
    fedges[1:-1] = 0.5*(flist[:-1] + flist[1:])
    fedges[0]    = flist[0]  - 0.5*(dflist[0])
    fedges[-1]   = flist[-1] + 0.5*(dflist[-1])
    
    def bin_index_from_freq_edges(fedges, nu_obs):
        iz = np.searchsorted(fedges, nu_obs, side='right') - 1
        return np.clip(iz, 0, len(fedges)-2)
    
    freq_obs = ha/(1 + np.array(z_array))*1e9 #[Hz]
    iz_Ha = bin_index_from_freq_edges(fedges, freq_obs)
    Ha_intensity = Ha / geom / dflist[iz_Ha]  # erg/s/cm2/Hz
    
    freq_obs = pa/(1 + np.array(z_array))*1e9 #[Hz]
    iz_Pa = bin_index_from_freq_edges(fedges, freq_obs)
    Pa_intensity = Pa / geom / dflist[iz_Pa]  # erg/s/cm2/Hz
    
    freq_obs = Oiii/(1 + np.array(z_array)) * 1e9 #[Hz]
    iz_OIII = bin_index_from_freq_edges(fedges, freq_obs)
    OIII_intensity = OIII / geom / dflist[iz_OIII]  # erg/s/cm2/Hz
    
    outdir = '../output/catalog/spectra/SIDES'
    
    save_meta_parquet(
        outdir, Mstar, sfr, ra, dec, photoz,
        iz_Ha, iz_Pa, iz_OIII,
        Ha, Pa, OIII,
        Ha_intensity, Pa_intensity, OIII_intensity,
        f115w_pred, f150w_pred, f277w_pred, f444w_pred,
        flist, dflist
        )

    # 2) 輝線スペクトル 2D を作成（点置き or LSF分配）
    fnu_lines_jy = build_line_spectra_jy(
        Ngal=log_spectrum.shape[0], flen=log_spectrum.shape[1],
        iz_Ha=iz_Ha, I_Ha=Ha_intensity,
        iz_Pa=iz_Pa, I_Pa=Pa_intensity,
        iz_OIII=iz_OIII, I_OIII=OIII_intensity,
        gaussian_sigma_bins=None        # 例：None=点置き、0.7等でガウシアン分配
        )

    # 3) HDF5 保存（連続光＋輝線＋合成）
    fnu = log_spectrum  # erg/s/cm2/Hz
    save_spectra_h5(outdir, flist, dflist, fnu_cont_jy=fnu, fnu_lines_jy=fnu_lines_jy)
    
def abmag_to_fnu(mag):
    """
    AB magnitude を erg/s/cm^2/Hz の F_nu に変換する
    """
    f_nu = 10.0**(-0.4*(mag + 48.6))
    return f_nu

def generate_mocks(names, params):
    init_worker_fsps()
    
    flist, dflist = frequency(params)
    Nz = len(flist)
    
    Nx = params.radius * 3600.0/params.resolution
    
    f_sample = "/mnt/data_cat3/yuka/repository/LIM_mock/theta_samples_100000_seed1_mag.npz"
    data = np.load(f_sample, allow_pickle=True)
    samples = data["samples"]
    logM_samples = data["logM_samples"]
    logM_formed = data['logM_formed']
    tage = data['tage']
    #sfr = data['SFR']
    log10SFR_bins = data['log10SFR_bins']
    t_edge_bins = data['t_edge_bins']
    redshift = samples[:,-1]
    SFR = 10**log10SFR_bins[:,0]

    print("Load {}".format(f_sample))
    
    # create (z, logM) bins
    Nbin = 20
    z_edges = np.linspace(0, 4, Nbin)
    logM_edges = np.linspace(7, 11, Nbin)

    z_idx = np.digitize(redshift, z_edges) - 1
    logM_idx = np.digitize(logM_samples, logM_edges) - 1

    z_idx = np.clip(z_idx, 0, Nbin - 1)
    logM_idx = np.clip(logM_idx, 0, Nbin - 1)

    from collections import defaultdict
    bin_to_indices = defaultdict(list)

    for i in range(len(samples)):
        key = (z_idx[i], logM_idx[i])
        bin_to_indices[key].append(i)

    for key in bin_to_indices:
        bin_to_indices[key] = np.array(bin_to_indices[key])
    
    with Pool(processes=20) as pool:
        results = pool.map(Get_Sides, names)
    ra, dec, photoz, Mhalo, Mstar = zip(*results)
    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    photoz = np.concatenate(photoz)
    Mhalo = np.concatenate(Mhalo)
    Mstar = np.concatenate(Mstar)
    
    ra_min = np.min(ra)
    ra_max = ra_min + params.radius
    dec_min = np.min(dec)
    dec_max = dec_min + params.radius
    
    selection = (ra>ra_min)&(ra<ra_max)&(dec>dec_min)&(dec<dec_max)
    ra = ra[selection]
    dec = dec[selection]
    photoz = photoz[selection]
    Mhalo = Mhalo[selection]
    Mstar = Mstar[selection]
    logM = np.log10(Mstar)
        
    ix = ((ra - ra_min)*3600.0/params.resolution).astype(np.int32)
    iy = ((dec - dec_min)*3600.0/params.resolution).astype(np.int32)
    
    idx = sample_ids(photoz, logM, bin_to_indices, z_edges, logM_edges)
    
    sfr = SFR[idx]
    area = params.radius**2
    Ha = GetEL(sfr, photoz, area=area, name='Halpha') #erg/s
    OIII = GetEL(sfr, photoz, area=area, name='OIII') #erg/s
    log_Ha = np.log10(Ha)
    log_Pa = log_Ha - 0.60
    Pa = 10**log_Pa #erg/s
    
    geom = 4.0 * np.pi * (cosmo.luminosity_distance(np.array(photoz)).to(u.cm).value)**2  # cm2
    
    fedges = np.empty(len(flist) + 1, dtype=float)
    fedges[1:-1] = 0.5*(flist[:-1] + flist[1:])
    fedges[0]    = flist[0]  - 0.5*(dflist[0])
    fedges[-1]   = flist[-1] + 0.5*(dflist[-1])
    
    def bin_index_from_freq_edges(fedges, nu_obs):
        iz = np.searchsorted(fedges, nu_obs, side='right') - 1
        return np.clip(iz, 0, len(fedges)-2)
    
    freq_obs = ha/(1 + np.array(photoz))*1e9 #[Hz]
    iz_Ha = bin_index_from_freq_edges(fedges, freq_obs)
    Ha_intensity = Ha / geom / dflist[iz_Ha]  # erg/s/cm2/Hz
    
    freq_obs = pa/(1 + np.array(photoz))*1e9 #[Hz]
    iz_Pa = bin_index_from_freq_edges(fedges, freq_obs)
    Pa_intensity = Pa / geom / dflist[iz_Pa]  # erg/s/cm2/Hz
    
    freq_obs = Oiii/(1 + np.array(photoz)) * 1e9 #[Hz]
    iz_OIII = bin_index_from_freq_edges(fedges, freq_obs)
    OIII_intensity = OIII / geom / dflist[iz_OIII]  # erg/s/cm2/Hz
    
    outdir = '../output/catalog/spectra/SIDES'
            
    save_meta_parquet1(outdir, Mstar, sfr, ra, dec, ix, iy, photoz,
                      iz_Ha, iz_Pa, iz_OIII,
                      Ha, Pa, OIII,
                      Ha_intensity, Pa_intensity, OIII_intensity,
                      idx)
    
    
    
    
    
    

def main():
    init_worker_fsps()
    ra, dec, obs_redshift, Mhalo, mstar = Get_Sides()
    f_sample = "../../theta_samples_100000.npz"
    data = np.load(f_sample, allow_pickle=True)
    samples = data["samples"]
    logM_samples = data["logM_samples"]
    redshift = samples[:,-1]

    print("Load {}".format(f_sample))

    # create (z, logM) bins
    Nbin = 20
    z_edges = np.linspace(0, 4, Nbin)
    logM_edges = np.linspace(7, 11, Nbin)

    z_idx = np.digitize(redshift, z_edges) - 1
    logM_idx = np.digitize(logM_samples, logM_edges) - 1

    z_idx = np.clip(z_idx, 0, Nbin - 1)
    logM_idx = np.clip(logM_idx, 0, Nbin - 1)

    from collections import defaultdict
    bin_to_indices = defaultdict(list)

    for i in range(len(samples)):
        key = (z_idx[i], logM_idx[i])
        bin_to_indices[key].append(i)

    for key in bin_to_indices:
        bin_to_indices[key] = np.array(bin_to_indices[key])

    z_array = obs_redshift[::10]
    M_array = np.log10(mstar)[::10]
    
    print(f'number of galaxies {len(z_array)}')
    params = sample_params(z_array, M_array, samples, bin_to_indices, z_edges, logM_edges)
    
    start = time.process_time()
    
    worker_iter = zip(params, z_array, M_array)
    with Pool(processes=20) as pool:
        #results = pool.map(emulate_spectrum, worker_iter)
        results = pool.map(get_true_spectrum, worker_iter)
    waves, log_abs_spectrum= zip(*results)
    log_abs_spectrum = np.array(log_abs_spectrum) # (N_samples, N_lambda) [erg/s/Hz]
    global wavelength
    wavelength = np.array(waves[0])
    
    end = time.process_time()
    print(f'time: {end-start}')
    
    filter_files = {}
    filter_names = ['F115W', 'F150W', 'F277W', 'F444W']
    for fname in filter_names:
        if 'Spitzer' in fname:
            filter_files[fname] = f"/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125{fname[-3:].lower()}trans_full.txt"
        else:
            filter_files[fname] = f"/mnt/data_cat3/yuka/data/nircam_throughputs/mean_throughputs/{fname}_May2024_mean_system_throughput.txt"

    global filters
    filters = {}
    for fname, path in filter_files.items():
        nu_flt, Tnu_flt = load_filter_throughput_txt(fname, path)
        denom_T = np.trapz(Tnu_flt, nu_flt)
        filters[fname] = {
            "nu": nu_flt,
            "T":  Tnu_flt,
            "denom_T": max(denom_T, 1e-300),
        }

    worker_iter = zip(log_abs_spectrum, z_array)
    with Pool(processes=20) as pool:
        results = pool.map(get_magnitudes_from_spectrum, worker_iter)
    f115w_pred, f150w_pred, f277w_pred, f444w_pred = zip(*results)
    f115w_pred = np.array(f115w_pred)
    f150w_pred = np.array(f150w_pred)
    f277w_pred = np.array(f277w_pred)
    f444w_pred = np.array(f444w_pred)
    
    plot_LF(f115w_pred, 'F115W')
    plot_LF(f150w_pred, 'F150W')
    plot_LF(f277w_pred, 'F277W')
    plot_LF(f444w_pred, 'F444W')
    
def save_meta_parquet(outdir, Mstar, sfr, ra, dec, redshift,
                      iz_Ha, iz_Pa, iz_OIII,
                      L_Ha, L_Pa, L_OIII, 
                      F_Ha, F_Pa, F_OIII,
                      f115w, f150w, f277w, f444w):
    os.makedirs(outdir, exist_ok=True)

    df = pd.DataFrame({
        "gal_id"       : np.arange(len(Mstar), dtype=np.int64),
        "ra"            : np.asarray(ra),
        "dec"            : np.asarray(dec),
        "z_obs"        : np.asarray(redshift),
        "Mstar"        : np.asarray(Mstar),
        "SFR"          : np.asarray(sfr),
        # line info per galaxy（落ちたbinと面輝度）
        "Ha_iz"        : np.asarray(iz_Ha,   dtype=np.int32),
        "Ha_L"   : np.asarray(L_Ha,    dtype=np.float32),
        "Ha_F"   : np.asarray(F_Ha,    dtype=np.float32),
        "Pa_iz"        : np.asarray(iz_Pa,   dtype=np.int32),
        "Pa_L"   : np.asarray(L_Pa,    dtype=np.float32),
        "Pa_F"   : np.asarray(F_Pa,    dtype=np.float32),
        "OIII_iz"      : np.asarray(iz_OIII, dtype=np.int32),
        "OIII_L" : np.asarray(L_OIII,  dtype=np.float32),
        "OIII_F" : np.asarray(F_OIII,  dtype=np.float32),
        # JWST magnitude
        "f115w"  : np.asarray(f115w,  dtype=np.float32),
        "f150w"  : np.asarray(f150w,  dtype=np.float32),
        "f277w"  : np.asarray(f277w,  dtype=np.float32),
        "f444w"  : np.asarray(f444w,  dtype=np.float32),
    })
    path = outdir + "_galaxies_all_meta.parquet"
    #df.to_parquet(os.path.join(outdir, "galaxies_meta_TNG.parquet"), index=False)
    df.to_parquet(path, index=False)
    
def save_meta_parquet1(outdir, Mstar, sfr, ra, dec, ix, iy, redshift,
                      iz_Ha, iz_Pa, iz_OIII,
                      L_Ha, L_Pa, L_OIII, 
                      F_Ha, F_Pa, F_OIII,
                      ID):
    os.makedirs(outdir, exist_ok=True)

    df = pd.DataFrame({
        "gal_id"       : np.asarray(ID, dtype=np.int64),
        "ra"            : np.asarray(ra),
        "dec"            : np.asarray(dec),
        "ix"            : np.asarray(ix, dtype=np.int32),
        "iy"            : np.asarray(iy, dtype=np.int32),
        "z_obs"        : np.asarray(redshift),
        "Mstar"        : np.asarray(Mstar),
        "SFR"          : np.asarray(sfr),
        # line info per galaxy（落ちたbinと面輝度）
        "Ha_iz"        : np.asarray(iz_Ha,   dtype=np.int32),
        "Ha_L"   : np.asarray(L_Ha,    dtype=np.float32),
        "Ha_F"   : np.asarray(F_Ha,    dtype=np.float32),
        "Pa_iz"        : np.asarray(iz_Pa,   dtype=np.int32),
        "Pa_L"   : np.asarray(L_Pa,    dtype=np.float32),
        "Pa_F"   : np.asarray(F_Pa,    dtype=np.float32),
        "OIII_iz"      : np.asarray(iz_OIII, dtype=np.int32),
        "OIII_L" : np.asarray(L_OIII,  dtype=np.float32),
        "OIII_F" : np.asarray(F_OIII,  dtype=np.float32),
    })
    path = outdir + "_galaxies_all_meta.parquet"
    #df.to_parquet(os.path.join(outdir, "galaxies_meta_TNG.parquet"), index=False)
    df.to_parquet(path, index=False)
    
    
def build_line_spectra_jy(Ngal, flen, iz_Ha, I_Ha, iz_Pa, I_Pa, iz_OIII, I_OIII,
                          gaussian_sigma_bins=None):
    """
    各銀河の輝線スペクトル（[Jy/sr]）を (Ngal, Nz) に展開。
    - gaussian_sigma_bins を与えると、各線をそのbin中心にガウシアン分配（チャネル単位のσ）。
      None なら 1チャネルに点置き。
    """
    lines = np.zeros((Ngal, flen), dtype=np.float32)

    if gaussian_sigma_bins is None:
        # 点置き（bin外/NaNを弾く）
        for iz, I in [(iz_Ha, I_Ha), (iz_Pa, I_Pa), (iz_OIII, I_OIII)]:
            m = np.isfinite(iz) & np.isfinite(I) & (iz >= 0) & (iz < flen)
            np.add.at(lines, (np.where(m)[0], iz[m].astype(int)), I[m].astype(np.float32))
        return lines

    # ガウシアン分配
    win = int(np.ceil(gaussian_sigma_bins*5))  # ±2.5σくらいまで
    idxs = np.arange(flen)
    for iz, I in [(iz_Ha, I_Ha), (iz_Pa, I_Pa), (iz_OIII, I_OIII)]:
        m = np.isfinite(iz) & np.isfinite(I) & (iz >= 0) & (iz < flen)
        g_ids = np.where(m)[0]
        for g in g_ids:
            c = int(iz[g])
            lo = max(0, c - win)
            hi = min(flen - 1, c + win)
            w = np.exp(-0.5*((idxs[lo:hi+1]-c)/gaussian_sigma_bins)**2)
            w /= (w.sum() or 1.0)
            lines[g, lo:hi+1] += (I[g] * w).astype(np.float32)
    return lines

def save_spectra_h5(outdir, flist, dflist, fnu_cont_jy, fnu_lines_jy=None):
    path = outdir + "_spectra_1deg.h5"
    #path = os.path.join(outdir, "spectra_TNG.h5")
    with h5py.File(path, "w") as h5:
        h5.create_dataset("flist_Hz",  data=flist,  compression="gzip")
        h5.create_dataset("dflist_Hz", data=dflist, compression="gzip")
        h5.create_dataset(
            "fnu_cont", data=fnu_cont_jy.astype(np.float32),
            compression="gzip", compression_opts=4,
            chunks=(min(1024, fnu_cont_jy.shape[0]), min(2048, fnu_cont_jy.shape[1]))
        )
        if fnu_lines_jy is not None:
            h5.create_dataset(
                "fnu_lines", data=fnu_lines_jy.astype(np.float32),
                compression="gzip", compression_opts=4,
                chunks=(min(1024, fnu_lines_jy.shape[0]), min(2048, fnu_lines_jy.shape[1]))
            )
            fnu_total = fnu_cont_jy + fnu_lines_jy
            h5.create_dataset(
                "fnu_total", data=fnu_total.astype(np.float32),
                compression="gzip", compression_opts=4,
                chunks=(min(1024, fnu_total.shape[0]), min(2048, fnu_total.shape[1]))
            )
    
if __name__ == "__main__":
    main()