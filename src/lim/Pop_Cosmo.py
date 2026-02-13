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

import time

from lim.AbundanceMatch import GetStellarMass
from lim.AbundanceMatch import GetEL, GetEL_nb

from lim.utils import wtheta_landy_szalay

from scipy.stats import binned_statistic

log_Lsun = 33.58297 # log10 of solar luminosity in erg/s/Hz

# グローバル変数（ワーカーごとに中身が別になる）
prospector_alpha_UV = None
prospector_alpha_opticalNIR = None
prospector_alpha_IR = None

sp = None
wavelength = None
filters = None
#log_Lsun = np.log10(L_sun.cgs.value)

def Get_Sides(name):
    file = '/mnt/data_cat3/yuka/data/SIDES/' + f'pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits'
    if os.path.exists(file) == False:
        print(f'No file: pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits')
        sys.exit(1)
    else:
        print(f'Load file: pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits')
    hdu = fits.open(file)
    data = hdu[1].data

    redshift = data['redshift']
    Mstar = data['Mstar'] #[Msun]
    SFR = data['SFR'] #[Msun/yr]
    Mhalo = data['Mhalo'] #[Msun]
    ra = data['ra']
    dec = data['dec']
    hdu.close()
    
    #Mstar = GetStellarMass(np.log10(Mhalo), redshift)

    #return redshift[redshift<5.85], SFR[redshift<5.85], Mhalo[redshift<5.85]
    return ra[(redshift > 0)&(redshift < 6.5)], dec[(redshift > 0)&(redshift < 6.5)], redshift[(redshift > 0)&(redshift < 6.5)], Mstar[(redshift > 0)&(redshift < 6.5)], SFR[(redshift > 0)&(redshift < 6.5)], Mhalo[(redshift > 0)&(redshift < 6.5)]

def init_worker_fsps():
    """各プロセス起動時に一度だけ呼ばれる。"""
    global sp
    sp = fsps.StellarPopulation(
        zcontinuous=1,
        add_neb_emission=True,
        add_neb_continuum=True,
        add_dust_emission=True,
        redshift_colors=1,
        sfh=3,
        imf_type=1,
        dust_type=0,
    )

    
    

#def init_worker():
    #"""各ワーカーで一回だけ呼ばれる初期化関数"""
    #global prospector_alpha_UV, prospector_alpha_opticalNIR, prospector_alpha_IR

    #prospector_alpha_UV = Speculator(
        #restore=True,
        #restore_filename='/mnt/data_cat3/yuka/repository/speculator/trained_models/prospector_alpha/100_400/model'
    #)
    #prospector_alpha_opticalNIR = Speculator(
        #restore=True,
        #restore_filename='/mnt/data_cat3/yuka/repository/speculator/trained_models/prospector_alpha/400_1100/model'
    #)
    #prospector_alpha_IR = Speculator(
        #restore=True,
        #restore_filename='/mnt/data_cat3/yuka/repository/speculator/trained_models/prospector_alpha/1100_30000/model'
    #)
    
def emulate_spectrum(args):
    """pool.map に投げる関数（トップレベルに定義すること）"""
    global prospector_alpha_UV, prospector_alpha_opticalNIR, prospector_alpha_IR

    # ここで prospector_* を使ってスペクトル or マグを計算
    # 例えば:
    # spectrum = calc_spectrum(prospector_alpha_UV,
    #                          prospector_alpha_opticalNIR,
    #                          prospector_alpha_IR,
    #                          params, z, M)
    # return spectrum
    return calc_spectrum(
        prospector_alpha_UV,
        prospector_alpha_opticalNIR,
        prospector_alpha_IR,
        args
    )

def load_data():
    path = '/mnt/data_cat3/yuka/data/Cosmos_Web/COSMOSWeb_mastercatalog_v1_photom_primary.fits'
    hdu = fits.open(path)
    data = hdu[1].data

    f115w = data['mag_model_f115w']
    f150w = data['mag_model_f150w']
    f277w = data['mag_model_f277w']
    f444w = data['mag_model_f444w']

    path = '/mnt/data_cat3/yuka/data/Cosmos_Web/COSMOSWeb_mastercatalog_v1_cigale.fits'
    hdu = fits.open(path)
    data = hdu[1].data
    print(repr(hdu[1].header))

    total_SFH = data['sfh_integrated']
    tage = data['age_form']
    SFH1 = data['sfh_sfr_bin1']
    SFH1_time = data['sfh_time_bin1'] #look back time in Myr
    SFH2 = data['sfh_sfr_bin2']
    SFH2_time = data['sfh_time_bin2']
    SFH3 = data['sfh_sfr_bin3']
    SFH3_time = data['sfh_time_bin3']
    SFH4 = data['sfh_sfr_bin4']
    SFH4_time = data['sfh_time_bin4']
    SFH5 = data['sfh_sfr_bin5']
    SFH5_time = data['sfh_time_bin5']
    SFH6 = data['sfh_sfr_bin6']
    SFH6_time = data['sfh_time_bin6']
    SFH7 = data['sfh_sfr_bin7']
    SFH7_time = data['sfh_time_bin7']
    SFH8 = data['sfh_sfr_bin8']
    SFH8_time = data['sfh_time_bin8']
    SFH9 = data['sfh_sfr_bin9']
    SFH9_time = data['sfh_time_bin9']

    zmeta = data['metallicity']

    Mstar = data['mass']
    Mstar_err = data['mass_err']

    #Mstar *= 10
    hdu.close()

    path = '/mnt/data_cat3/yuka/data/Cosmos_Web/COSMOSWeb_mastercatalog_v1_lephare.fits'
    hdu = fits.open(path)
    data = hdu[1].data
    print(repr(hdu[1].header))

    photoz = data['zfinal']
    Type = data['type']
    tage1 = data['age_med']

    hdu.close()

    selection = (Type==0)&(photoz>0)&(photoz<5.85)&(Mstar>1e7)
    return f115w[selection], f150w[selection], f277w[selection], f444w[selection], total_SFH[selection], photoz[selection], Type[selection], Mstar[selection]

def load_filter_throughput_txt(fname, path, col_wave=0, col_thr=1):
    """path: λ vs T(λ)。λ単位は 'um' or 'A'。返り値は (nu, Tnu)（ν昇順）"""
    if 'Spitzer' in fname:
        arr = np.loadtxt(path)
        lam = arr[:, col_wave].astype(float)  # 波長
        Tlam = arr[:, col_thr].astype(float)  # 透過
        lam_A = lam * 1e4
    else:
        arr = np.loadtxt(path, skiprows=1)
        lam = arr[:, col_wave].astype(float)  # 波長
        Tlam = arr[:, col_thr].astype(float)  # 透過
        lam_A = lam * 1e4  # Å

    c_A_per_s = 2.99792458e18
    nu = c_A_per_s / lam_A
    order = np.argsort(nu)
    return nu[order], Tlam[order]

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

def sample_SEDs(redshift, logM, SED, bin_to_indices, z_edges, logM_edges, f115w_pred, f150w_pred, f277w_pred, f444w_pred):

    z_bins = np.digitize(redshift, z_edges) - 1
    logM_bins = np.digitize(logM, logM_edges) - 1

    z_bins = np.clip(z_bins, 0, len(z_edges)-1)
    logM_bins = np.clip(logM_bins, 0, len(logM_edges)-1)

    B = len(redshift)
    max_zbin = len(z_edges) - 1
    max_mbin = len(logM_edges) - 1

    sample_list = []
    f115w = []
    f150w = []
    f277w = []
    f444w = []
    
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
            
            f115w.append(f115w_pred[chosen_idx])
            f150w.append(f150w_pred[chosen_idx])
            f277w.append(f277w_pred[chosen_idx])
            f444w.append(f444w_pred[chosen_idx])

    spectrum = np.stack(sample_list, axis=0) # (B, num_params)
    f115w = np.stack(f115w, axis=0)
    f150w = np.stack(f150w, axis=0)
    f277w = np.stack(f277w, axis=0)
    f444w = np.stack(f444w, axis=0)

    return spectrum, f115w, f150w, f277w, f444w

def sample_mags(redshift, logM, bin_to_indices, z_edges, logM_edges, f115w_pred, f150w_pred, f277w_pred, f444w_pred):

    z_bins = np.digitize(redshift, z_edges) - 1
    logM_bins = np.digitize(logM, logM_edges) - 1

    z_bins = np.clip(z_bins, 0, len(z_edges)-1)
    logM_bins = np.clip(logM_bins, 0, len(logM_edges)-1)

    B = len(redshift)
    max_zbin = len(z_edges) - 1
    max_mbin = len(logM_edges) - 1

    f115w = []
    f150w = []
    f277w = []
    f444w = []
    
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
            
            f115w.append(f115w_pred[chosen_idx])
            f150w.append(f150w_pred[chosen_idx])
            f277w.append(f277w_pred[chosen_idx])
            f444w.append(f444w_pred[chosen_idx])

    f115w = np.stack(f115w, axis=0)
    f150w = np.stack(f150w, axis=0)
    f277w = np.stack(f277w, axis=0)
    f444w = np.stack(f444w, axis=0)

    return f115w, f150w, f277w, f444w


def sample_ELs(redshift, logM, SFR, bin_to_indices, z_edges, logM_edges):

    z_bins = np.digitize(redshift, z_edges) - 1
    logM_bins = np.digitize(logM, logM_edges) - 1

    z_bins = np.clip(z_bins, 0, len(z_edges)-1)
    logM_bins = np.clip(logM_bins, 0, len(logM_edges)-1)

    B = len(redshift)
    max_zbin = len(z_edges) - 1
    max_mbin = len(logM_edges) - 1

    sfr_list = []
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
            sfr_list.append( SFR[chosen_idx] )

    sfr_all = np.array(sfr_list)

    return sfr_all

def distance_modulus(z,cosmo=cosmo):
    """
    PyTorch routine for computing distance modulus.

    Parameters
    ----------
    z : torch.Tensor
        Redshifts to evaluate for.
    H0 : float, optional
        Hubble constant in km/s/Mpc. Defaults to Planck 2018.
    OmegaM : float, optional
        Matter density. Defaults to Planck 2018.

    Returns
    -------
    mu : torch.Tensor
        Distance modulus in magnitudes.

    See Also
    --------
    `comoving_distance` : Approximate distance integral.
    """
    DL = cosmo.luminosity_distance(z).to(u.Mpc).value
    return 5.0*np.log10(DL) + 25.0

def _calculate_SFH_time(index, redshift):
    t_univ = cosmo.age(redshift).to(u.Myr).value # in Gyr
    t = 100.0 + (index - 3.0)/4.0 * (0.85 * t_univ - 100.0)
    return t # in Gyr

def _get_SFHs(ratios, SFH_time, Mstar):
    dt = np.diff(SFH_time)   # duration of each bin [Gyr]
    M_form = 10**Mstar    # total formed stellar mass [Msun]

    # compute multiplicative cumulative ratios
    mult = np.ones(len(dt))
    for i in range(1, len(dt)):
        mult[i] = mult[i-1] * ratios[i-1]

    # solve for SFR_1
    SFR1 = M_form / np.sum(mult * dt)

    # compute all SFRs
    SFR = SFR1 * mult
    SFH = SFR / (dt*1e9)
    return SFR, SFH

def _t_age(Mstar, SFR, SFH_time):
    tage = 0.0
    for i in range(6):
        tage += 1/2.0 * (SFH_time[i+1]**2 - SFH_time[i]**2) * (SFR[i])/Mstar
    return tage

def redshift_from_cosmic_time(t_gyr):
    # t_gyr: 宇宙の始まりからの時間 [Gyr]
    t = t_gyr * u.Gyr
    z = z_at_value(cosmo.age, t)  # cosmo.age(z) = 宇宙年齢
    z = z.value
    return z

def calc_spectrum(prospector_UV, prospector_opt, prospector_IR, param_all):
    param, redshift_true, Mstar_true = param_all
    N_pred, zmeta_pred, sfr1_pred, sfr2_pred, sfr3_pred, sfr4_pred, sfr5_pred, sfr6_pred,  tau2_pred, dust_n_pred, tau1_pred, fAGN_pred, tauAGN_pred, gasZ_pred, gasU_pred, redshift_pred = param
    ln_zmeta = zmeta_pred
    
    ln_sfr1 = sfr1_pred
    ln_sfr2 = sfr2_pred
    ln_sfr3 = sfr3_pred
    ln_sfr4 = sfr4_pred
    ln_sfr5 = sfr5_pred
    ln_sfr6 = sfr6_pred
    
    tau2 = np.sqrt(tau2_pred)
    
    ln_gasZ = gasZ_pred
    
    ratios = 10 ** np.array([sfr1_pred, sfr2_pred, sfr3_pred, sfr4_pred, sfr5_pred, sfr6_pred])
    SFH_time = np.array([0,30])
    SFH_time = np.append(SFH_time, calculate_SFH_time(np.array([3,4,5,6,7]),redshift_pred))
    SFH_time = np.array(SFH_time)/1e3  # in Gyr

    SFR, SFH = get_SFHs(ratios, SFH_time, Mstar_true)
    tage = t_age(10**Mstar_true, SFR, SFH_time)

    redshift_eff = redshift_from_cosmic_time(tage)
    print(f"effective redshift : {redshift_eff}")


    

    theta = np.array([ln_zmeta, # log metallicity
                      ln_sfr1, ln_sfr2, ln_sfr3, ln_sfr4, ln_sfr5, ln_sfr6, # log SFH ratios
                      tau2, # square root of diffuse dust optical depth
                      dust_n_pred, # dust attenuation index
                      tau1_pred, # ratio of birth cloud vs diffuse optical depths   
                      fAGN_pred, # fraction of bolometric luminosity from AGN
                      tauAGN_pred, # optical depth of AGN dust taurus
                      ln_gasZ, # log gas phase metallicity
                      redshift_eff # redshift (used to set lookback time only - see notes above)
                      ])

    UV_spectrum = prospector_UV.log_spectrum_(theta) + Mstar_true + log_Lsun # log spectrum [erg/s/Hz]
    opticalNIR_spectrum = prospector_opt.log_spectrum_(theta) + Mstar_true  + log_Lsun # log spectrum [erg/s/Hz]
    IR_spectrum = prospector_IR.log_spectrum_(theta) + Mstar_true + log_Lsun # log spectrum [erg/s/Hz]

    spectrum = np.concatenate([UV_spectrum, opticalNIR_spectrum, IR_spectrum], axis=0)

    return spectrum

def get_magnitudes_from_spectrum(params):
    """
    wavelength: 1D [Angstrom], rest-frame
    log_abs_spectrum: 2D (Nobj, Nwave), log10(L_lambda) [erg/s/Angstrom] の想定
    redshift: 1D (Nobj,)
    """

    global wavelength
    log_abs_spectrum, redshift = params
    #wavelength = np.asarray(wavelength)               # (Nwave,)
    log_abs_spectrum = np.asarray(log_abs_spectrum)   # (Nobj, Nwave)
    redshift = np.asarray(redshift)                   # (Nobj,)

    lumi_dist = cosmo.luminosity_distance(redshift).to(u.cm).value  # (Nobj,)
    four_pi = 4.0 * np.pi

    # --- rest-frame frequency grid ---
    c_A_per_s = 2.99792458e18  # [Angstrom/s]
    lam_A = wavelength         # [Angstrom]
    nu = c_A_per_s / lam_A     # [Hz]

    # 周波数が昇順になるようにソート（np.interp 用）
    if not np.all(np.diff(nu) > 0):
        order = np.argsort(nu)
        nu = nu[order]
        lam_A = lam_A[order]
        log_abs_spectrum = log_abs_spectrum[order]

    # --- L_lambda -> L_nu への変換を log でやる ---
    # log10 L_nu = log10 L_lambda + 2 log10(lambda) - log10(c)
    log10_lambda = np.log10(lam_A)           # (Nwave,)
    log10_c      = np.log10(c_A_per_s)
    log_L_nu_all = log_abs_spectrum
    #log_L_nu_all = log_abs_spectrum + 1.0 * log10_lambda[None, :] + 4.0 - log10_c

    # 出力用
    f115w_pred = f150w_pred = f277w_pred = f444w_pred = None
    
    z = redshift
    D_L = lumi_dist

    # この銀河の log L_nu(ν_rest) スペクトル
    log_L_nu = log_L_nu_all  # (Nwave,)

    # rest -> observed frame
    nu_obs = nu / (1.0 + z)     # (Nwave,)
    # F_nu(obs) = L_nu(rest) / [4π D_L^2 (1+z)]
    log_F_nu = (
        log_L_nu
        - 2.0 * np.log10(D_L)
        - np.log10(four_pi)
        - np.log10(1.0 + z)
    )
    F_nu_obs = 10.0**log_F_nu   # (Nwave,)

    for fname, filter_data in filters.items():
        nu_flt = filter_data["nu"]           # (Nflt,)  [Hz], 昇順想定
        Tnu_flt = filter_data["T"]           # (Nflt,)
        denom_T = filter_data["denom_T"]     # ∫ T(ν) dν のはず

        # フィルタ周波数で補間（nu_obs, F_nu_obs は昇順前提）
        f_nu_interp = np.interp(nu_flt, nu_obs, F_nu_obs, left=0.0, right=0.0)

        # フィルタ平均
        numerator = np.trapz(f_nu_interp * Tnu_flt, nu_flt)
        f_nu_flt = numerator / denom_T

        mag_AB = -2.5 * np.log10(f_nu_flt) - 48.6

        if fname == 'F115W':
            f115w_pred = mag_AB
        elif fname == 'F150W':
            f150w_pred = mag_AB
        elif fname == 'F277W':
            f277w_pred = mag_AB
        elif fname == 'F444W':
            f444w_pred = mag_AB

    return f115w_pred, f150w_pred, f277w_pred, f444w_pred

def get_magnitudes_from_spectrum_vec(wavelength, log_abs_spectrum, redshift, filters):
    """
    wavelength: 1D array [Angstrom], rest-frame
    log_abs_spectrum: 2D array (Nobj, Nwave), log10(L_nu) [erg/s/Hz]
    redshift: 1D array (Nobj,)
    filters: dict[fname] -> {"nu": nu_flt, "T": Tnu_flt, "denom_T": denom_T}
    """
    wavelength       = np.asarray(wavelength)         # (Nwave,)
    log_abs_spectrum = np.asarray(log_abs_spectrum)   # (Nobj, Nwave)
    redshift         = np.asarray(redshift)           # (Nobj,)

    Nobj, Nwave = log_abs_spectrum.shape

    # ---- luminosity distance (vectorized) ----
    D_L = cosmo.luminosity_distance(redshift).to(u.cm).value   # (Nobj,)
    four_pi = 4.0 * np.pi

    # ループ内で何度も log を取らないように、先にまとめて計算
    log10_4pi   = np.log10(four_pi)
    log10_DL2   = 2.0 * np.log10(D_L)           # (Nobj,)
    log10_1pz   = np.log10(1.0 + redshift)      # (Nobj,)

    # ---- rest-frame frequency grid (shared) ----
    c_A_per_s = 2.99792458e18
    nu_rest = c_A_per_s / wavelength    # Hz, (Nwave,)

    # 昇順でない場合はソート
    if not np.all(np.diff(nu_rest) > 0):
        order = np.argsort(nu_rest)
        nu_rest = nu_rest[order]
        log_abs_spectrum = log_abs_spectrum[:, order]

    # ---- log10 L_nu(rest) を補間する関数 ----
    interp_L_abs = interp1d(
        nu_rest,
        log_abs_spectrum,   # (Nobj, Nwave)
        kind="linear",
        axis=-1,            # 最後の次元が ν
        bounds_error=False,
        fill_value=0.0,
        assume_sorted=True,
    )

    mags_dict = {}

    for fname, fdata in filters.items():
        nu_flt  = np.asarray(fdata["nu"])   # (Nflt,)
        Tnu_flt = np.asarray(fdata["T"])    # (Nflt,)
        denom_T = fdata["denom_T"]          # scalar

        # ν_rest = (1+z)*ν_obs  → shape (Nobj, Nflt)
        nu_rest_target = (1.0 + redshift)[:, None] * nu_flt[None, :]

        # log10 L_nu(rest) at these ν_rest: (Nobj, Nflt)
        log_L_abs_nu_at = interp_L_abs(nu_rest_target)

        # log10 F_nu(obs) = log10 L_nu - log10(4π) - 2log10(D_L) - log10(1+z)
        # shape は (Nobj, Nflt) に揃う
        log_F_nu_at = (
            log_L_abs_nu_at
            - log10_4pi
            - log10_DL2[:, None]
            - log10_1pz[:, None]
        )

        # ここで初めて 10** するので、値のスケールはだいぶ抑えられる
        F_nu_at = np.power(10.0, log_F_nu_at)    # (Nobj, Nflt)

        # フィルター平均: ∫ F_ν T(ν) dν / ∫ T(ν) dν
        # F_nu_at * Tnu_flt[None, :] → (Nobj, Nflt)
        numerator = np.trapz(F_nu_at * Tnu_flt[None, :], nu_flt, axis=-1)  # (Nobj,)
        f_nu_flt  = numerator / denom_T

        mag_AB = -2.5 * np.log10(f_nu_flt) - 48.6
        mags_dict[fname] = mag_AB

    return (
        mags_dict.get("F115W", None),
        mags_dict.get("F150W", None),
        mags_dict.get("F277W", None),
        mags_dict.get("F444W", None),
    )
    
def plot_Nz_for_magnitude_bins(
    redshift, 
    mag, 
    mag_bins=[(18,19),(19,20),(20,21),(21,22),(22,23),(23,24)],
    label_prefix="JWST_F115W"
):

    redshift = np.asarray(redshift).ravel()
    mag      = np.asarray(mag).ravel()
    plt.figure(figsize=(8,6))
    bins = np.linspace(0.0, 3.0, 31)

    for (m1, m2) in mag_bins:
        sel = (mag >= m1) & (mag < m2)
        if sel.sum() == 0:
            continue

        z_sel = redshift[sel]

        plt.hist(z_sel, bins=bins, histtype='step',
                 label=f"{label_prefix}: {m1}–{m2} (N={sel.sum()})",
                 )

    plt.xlabel("Redshift")
    plt.ylabel("Normalized counts")
    plt.legend()
    plt.yscale('log')
    plt.ylim(1.0, 1e3)
    plt.xlim(0.0, 3.0)
    plt.title(f"N(z) per magnitude bin — {label_prefix}")
    plt.grid()
    plt.show()

    plt.savefig(f"/mnt/data_cat3/yuka/output/N_z_{label_prefix}.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/N_z_{label_prefix}.pdf")

def distance_modulus(z,cosmo=cosmo):
    """
    PyTorch routine for computing distance modulus.

    Parameters
    ----------
    z : torch.Tensor
        Redshifts to evaluate for.
    H0 : float, optional
        Hubble constant in km/s/Mpc. Defaults to Planck 2018.
    OmegaM : float, optional
        Matter density. Defaults to Planck 2018.

    Returns
    -------
    mu : torch.Tensor
        Distance modulus in magnitudes.

    See Also
    --------
    `comoving_distance` : Approximate distance integral.
    """
    DL = cosmo.luminosity_distance(z).to(u.Mpc).value
    return 5.0*np.log10(DL) + 25.0 

def get_true_spectrum(args):
    global sp
    theta_samples, logM_formed, mw_age, M_frac_bins, t_edge_bins = args
    N_pred, zmeta_pred, sfr1_pred, sfr2_pred, sfr3_pred, sfr4_pred, sfr5_pred, sfr6_pred,  tau2_pred, dust_n_pred, tau1_pred, fAGN_pred, tauAGN_pred, gasZ_pred, gasU_pred, redshift_pred = theta_samples

    t_diff = np.diff(t_edge_bins) * 1e9 # time within bin in yr
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

    remaining_ratio = sp.stellar_mass
    if (remaining_ratio<=0):
        print(f'remaining ratio is 0')
        #return wave, np.zeros(len(wave)), False, 0, 0, 0, 0
        return False, 0, 0, 0, 0
    else:
        log_spec -= np.log10(remaining_ratio) #erg/s/Hz per solar mass (remaining)
        #ratio = 10**(Mstar_true - Mstar_pred)
        f115w, f150w, f277w, f444w = sp.get_mags(tage=tage, redshift=redshift_pred, bands=['jwst_f115w', 'jwst_f150w', 'jwst_f277w', 'jwst_f444w'])
        #return wave, log_spec, True, f115w, f150w, f277w, f444w
        return True, f115w, f150w, f277w, f444w

def main():
    init_worker_fsps()
    f115w, f150w, f277w, f444w, total_SFH, photoz, Type, Mstar = load_data()
    # Further processing...
    # load pre-created data 
    f_sample = "../../theta_samples_100000_seed1_mag.npz"
    data = np.load(f_sample, allow_pickle=True)
    samples = data["samples"] #100000このsampleに対するSED parameter
    logM_samples = data["logM_samples"] #remaining stellar mass (log M_sun)
    logM_formed = data['logM_formed']
    tage = data['tage']
    Mfrac_bins = data['log10SFR_bins']
    t_edge_bins = data['t_edge_bins']
    redshift = samples[:,-1]

    print("Load {}".format(f_sample))

    start = time.process_time()
    
    worker_iter = zip(samples, logM_formed, tage, Mfrac_bins, t_edge_bins)
    with Pool(processes=20) as pool:
        results = pool.map(get_true_spectrum, worker_iter)
        #results = pool.map(get_true_spectrum, samples)
    waves, log_abs_spectrum, flag, f115w_pred, f150w_pred, f277w_pred, f444w_pred = zip(*results)
    flag = np.array(flag)
    log_abs_spectrum = np.array(log_abs_spectrum) # (N_samples, N_lambda) [erg/s/Hz]
    global wavelength
    wavelength = np.array(waves[0])
    
    f115w_pred = np.array(f115w_pred)
    f150w_pred = np.array(f150w_pred)
    f277w_pred = np.array(f277w_pred)
    f444w_pred = np.array(f444w_pred)
    
    end = time.process_time()
    print(f'time: {end-start}')

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
    #M_array = Get_Mstar(redshift, logM_samples, photoz, np.log10(Mstar))
    
    print(f'number of galaxies {len(z_array)}')
    log_spectrum, f115w_pred, f150w_pred, f277w_pred, f444w_pred = sample_SEDs(z_array, M_array, log_abs_spectrum, bin_to_indices, z_edges, logM_edges, f115w_pred, f150w_pred, f277w_pred, f444w_pred)
    log_spectrum += M_array[:, None]
    
    
    #filter_files = {}
    #filter_names = ['F115W', 'F150W', 'F277W', 'F444W']
    #for fname in filter_names:
        #if 'Spitzer' in fname:
            #filter_files[fname] = f"/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125{fname[-3:].lower()}trans_full.txt"
        #else:
            #filter_files[fname] = f"/mnt/data_cat3/yuka/data/nircam_throughputs/mean_throughputs/{fname}_May2024_mean_system_throughput.txt"

    #global filters
    #filters = {}
    #for fname, path in filter_files.items():
        #nu_flt, Tnu_flt = load_filter_throughput_txt(fname, path)
        #denom_T = np.trapz(Tnu_flt, nu_flt)
        #filters[fname] = {
            #"nu": nu_flt,
            #"T":  Tnu_flt,
            #"denom_T": max(denom_T, 1e-300),
        #}

    #worker_iter = zip(log_spectrum, z_array)
    #with Pool(processes=20) as pool:
        #results = pool.map(get_magnitudes_from_spectrum, worker_iter)
    #f115w_pred, f150w_pred, f277w_pred, f444w_pred = zip(*results)
    #f115w_pred = np.array(f115w_pred)
    #f150w_pred = np.array(f150w_pred)
    #f277w_pred = np.array(f277w_pred)
    #f444w_pred = np.array(f444w_pred)
    
    #f115w_pred, f150w_pred, f277w_pred, f444w_pred = get_magnitudes_from_spectrum(wavelength, log_abs_spectrum, z_array, filters)

    plot_Nz_for_magnitude_bins(z_array, f115w, label_prefix=f"F115W distribution")
    plot_Nz_for_magnitude_bins(z_array, f115w_pred, label_prefix=f"F115W distribution predicted")

    plt.figure()
    bins_m = np.arange(15, 31, 0.5)
    hist_mock, _ = np.histogram(f115w_pred, bins=bins_m)
    hist_JWST, _ = np.histogram(f115w, bins=bins_m)
    m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
    plt.scatter(m_cent, hist_mock, color='blue', label='predicted')
    plt.scatter(m_cent, hist_JWST, color='red', label='true')
    plt.legend()
    plt.yscale('log')
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_F115W.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_F115W.pdf")
    
    plt.figure()
    plt.scatter(f115w, f115w_pred, s=1, 
               c=z_array,   # 色：赤方偏移
               cmap='viridis',            # 好きなカラーマップに変更OK
               alpha=0.7)
    plt.colorbar(label='redshift (z)')
    x = np.linspace(15, 40, 50)
    y = x
    plt.plot(x, y, color='red', linestyle='dashed')
    plt.xlim(15, 40)
    plt.ylim(15, 40)
    plt.xlabel("true mag")
    plt.ylabel("predicted mag")
    plt.savefig(f"/mnt/data_cat3/yuka/output/magnitude_F115W.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/magnitude_F115W.pdf")
            

    plot_Nz_for_magnitude_bins(z_array, f150w, label_prefix=f"F150W distribution")
    plot_Nz_for_magnitude_bins(z_array, f150w_pred, label_prefix=f"F150W distribution predicted")

    plt.figure()
    bins_m = np.arange(15, 31, 0.5)
    hist_mock, _ = np.histogram(f150w_pred, bins=bins_m)
    hist_JWST, _ = np.histogram(f150w, bins=bins_m)
    m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
    plt.scatter(m_cent, hist_mock, color='blue', label='predicted')
    plt.scatter(m_cent, hist_JWST, color='red', label='true')
    plt.yscale('log')
    plt.legend()
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_F150W.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_F150W.pdf")
    
    plt.figure()
    plt.scatter(f150w, f150w_pred, s=1, 
               c=z_array,   # 色：赤方偏移
               cmap='viridis',            # 好きなカラーマップに変更OK
               alpha=0.7)
    plt.colorbar(label='redshift (z)')
    x = np.linspace(15, 40, 50)
    y = x
    plt.plot(x, y, color='red', linestyle='dashed')
    plt.xlim(15, 40)
    plt.ylim(15, 40)
    plt.xlabel("true mag")
    plt.ylabel("predicted mag")
    plt.savefig(f"/mnt/data_cat3/yuka/output/magnitude_F150W.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/magnitude_F150W.pdf")

    plot_Nz_for_magnitude_bins(z_array, f277w, label_prefix=f"F277W distribution")
    plot_Nz_for_magnitude_bins(z_array, f277w_pred, label_prefix=f"F277W distribution predicted")

    plt.figure()
    bins_m = np.arange(15, 31, 0.5)
    hist_mock, _ = np.histogram(f277w_pred, bins=bins_m)
    hist_JWST, _ = np.histogram(f277w, bins=bins_m)
    m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
    plt.scatter(m_cent, hist_mock, color='blue', label='predicted')
    plt.scatter(m_cent, hist_JWST, color='red', label='true')
    plt.yscale('log')
    plt.legend()
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_F277W.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_F277W.pdf")
    
    plt.figure()
    plt.scatter(f277w, f277w_pred, s=1, 
               c=z_array,   # 色：赤方偏移
               cmap='viridis',            # 好きなカラーマップに変更OK
               alpha=0.7)
    plt.colorbar(label='redshift (z)')
    x = np.linspace(15, 40, 50)
    y = x
    plt.plot(x, y, color='red', linestyle='dashed')
    plt.xlim(15, 40)
    plt.ylim(15, 40)
    plt.xlabel("true mag")
    plt.ylabel("predicted mag")
    plt.savefig(f"/mnt/data_cat3/yuka/output/magnitude_F277W.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/magnitude_F277W.pdf")

    plot_Nz_for_magnitude_bins(z_array, f444w, label_prefix=f"F444W distribution")
    plot_Nz_for_magnitude_bins(z_array, f444w_pred, label_prefix=f"F444W distribution predicted")

    plt.figure()
    bins_m = np.arange(15, 31, 0.5)
    hist_mock, _ = np.histogram(f444w_pred, bins=bins_m)
    hist_JWST, _ = np.histogram(f444w, bins=bins_m)
    m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
    plt.scatter(m_cent, hist_mock, color='blue', label='predicted')
    plt.scatter(m_cent, hist_JWST, color='red', label='true')
    plt.yscale('log')
    plt.legend()
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_F444W.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/LF_F444W.pdf")
    
    plt.figure()
    plt.scatter(f444w, f444w_pred, s=1, 
               c=z_array,   # 色：赤方偏移
               cmap='viridis',            # 好きなカラーマップに変更OK
               alpha=0.7)
    plt.colorbar(label='redshift (z)')
    x = np.linspace(15, 40, 50)
    y = x
    plt.plot(x, y, color='red', linestyle='dashed')
    plt.xlim(15, 40)
    plt.ylim(15, 40)
    plt.xlabel("true mag")
    plt.ylabel("predicted mag")
    plt.savefig(f"/mnt/data_cat3/yuka/output/magnitude_F444W.png")
    plt.savefig(f"/mnt/data_cat3/yuka/output/magnitude_F444W.pdf")

def match_distribution_monotone(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    n = len(b)

    # a の分布から n 個欲しいので、a が n より短い場合は「復元抽出」で増やす
    if len(a) < n:
        idx = np.random.randint(0, len(a), size=n)
        a_use = a[idx]
    else:
        # 長いなら n 個だけ使う（ランダムに選ぶと分布を保ちやすい）
        idx = np.random.choice(len(a), size=n, replace=False)
        a_use = a[idx]

    a_sorted = np.sort(a_use)
    b_order = np.argsort(b)

    a_prime = np.empty(n, dtype=a_sorted.dtype)
    a_prime[b_order] = a_sorted
    return a_prime


def Get_Mstar(sample_redshift, sample_Mstar, redshift, Mstar):
    redshift_bin = [0.4, 0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8]
    Mstar_ = Mstar
    for i, z in enumerate(redshift_bin):
        if i==0:
            zmin = 0.0
            zmax = redshift_bin[i]
        else:
            zmin = redshift_bin[i-1]
            zmax = redshift_bin[i]
        target_sel = (redshift>zmin) & (redshift<zmax)
        sample_sel = (sample_redshift>zmin) & (sample_redshift<zmax)
        Mstar_[target_sel] = match_distribution_monotone(sample_Mstar[sample_sel], Mstar[target_sel])
    return Mstar_

def test_mag(names):
    init_worker_fsps()
    
    with Pool(processes=20) as pool:
        results = pool.map(Get_Sides, names)
    #file = '/mnt/data_cat3/yuka/data/SIDES/' + f'pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits'
    #photoz, Mstar, SFR, Mhalo = Get_Sides(file)
    ra, dec, photoz_list, Mstar_list, sfr_list, mhalo_list = zip(*results)

    photoz = np.concatenate(photoz_list)
    SFR    = np.concatenate(sfr_list)
    Mhalo  = np.concatenate(mhalo_list)
    Mstar  = np.concatenate(Mstar_list)
    
    #Mstar = GetStellarMass(np.log10(Mhalo), photoz)

    f_sample = "/mnt/data_cat3/yuka/repository/LIM_mock/theta_samples_100000_seed1_mag.npz"
    data = np.load(f_sample, allow_pickle=True)
    samples = data["samples"] #100000このsampleに対するSED parameter
    logM_samples = data["logM_samples"] #remaining stellar mass (log M_sun)
    logM_formed = data['logM_formed']
    tage = data['tage']
    #sfr = data['SFR']
    log10SFR_bins = data['log10SFR_bins']
    t_edge_bins = data['t_edge_bins']
    redshift = samples[:,-1]

    print("Load {}".format(f_sample))

    start = time.process_time()
    
    #worker_iter = zip(params, z_array, total_array)
    worker_iter = zip(samples, logM_formed, tage, log10SFR_bins, t_edge_bins)
    with Pool(processes=20) as pool:
        results = pool.map(get_true_spectrum, worker_iter)
        #results = pool.map(get_true_spectrum, samples)
    #waves, log_abs_spectrum, flag, f115w_pred, f150w_pred, f277w_pred, f444w_pred = zip(*results)
    flag, f115w_pred, f150w_pred, f277w_pred, f444w_pred = zip(*results)
    flag = np.array(flag)
    #log_abs_spectrum = np.array(log_abs_spectrum) # (N_samples, N_lambda) [erg/s/Hz]
    #global wavelength
    #wavelength = np.array(waves[0])
    
    f115w_pred = np.array(f115w_pred)
    f150w_pred = np.array(f150w_pred)
    f277w_pred = np.array(f277w_pred)
    f444w_pred = np.array(f444w_pred)
    
    end = time.process_time()
    print(f'time: {end-start}')

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

    #M_array = Get_Mstar(redshift, logM_samples, photoz, np.log10(Mstar))
    M_array = np.log10(Mstar)
    z_array = photoz
    #M_array = np.log10(Mstar)
    
    print(f'number of galaxies {len(z_array)}')
    #log_spectrum, f115w_pred, f150w_pred, f277w_pred, f444w_pred = sample_SEDs(z_array, M_array, log_abs_spectrum, bin_to_indices, z_edges, logM_edges, f115w_pred, f150w_pred, f277w_pred, f444w_pred)
    f115w_pred, f150w_pred, f277w_pred, f444w_pred = sample_mags(z_array, M_array, bin_to_indices, z_edges, logM_edges, f115w_pred, f150w_pred, f277w_pred, f444w_pred)
    #log_spectrum += M_array[:, None]
    
    
    #filter_files = {}
    filter_names = ['F115W', 'F150W', 'F277W', 'F444W']
    #for fname in filter_names:
        #if 'Spitzer' in fname:
            #filter_files[fname] = f"/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125{fname[-3:].lower()}trans_full.txt"
        #else:
            #filter_files[fname] = f"/mnt/data_cat3/yuka/data/nircam_throughputs/mean_throughputs/{fname}_May2024_mean_system_throughput.txt"

    #global filters
    #filters = {}
    #for fname, path in filter_files.items():
        #nu_flt, Tnu_flt = load_filter_throughput_txt(fname, path)
        #denom_T = np.trapz(Tnu_flt, nu_flt)
        #filters[fname] = {
            #"nu": nu_flt,
            #"T":  Tnu_flt,
            #"denom_T": max(denom_T, 1e-300),
        #}

    #worker_iter = zip(log_spectrum, z_array)
    #with Pool(processes=20) as pool:
        #results = pool.map(get_magnitudes_from_spectrum, worker_iter)
    #f115w_pred, f150w_pred, f277w_pred, f444w_pred = zip(*results)
    #f115w_pred = np.array(f115w_pred)
    #f150w_pred = np.array(f150w_pred)
    #f277w_pred = np.array(f277w_pred)
    #f444w_pred = np.array(f444w_pred)

    def Get_LF(filter_name):
        if filter_name=='F115W':
            m_ab = f115w_pred
        elif filter_name=='F150W':
            m_ab = f150w_pred
        elif filter_name=='F227W':
            m_ab = f277w_pred
        else:
            m_ab =f444w_pred
        
        bins_m = np.arange(15, 28, 0.5)
        hist_m, _ = np.histogram(m_ab, bins=bins_m)
        dNdm = hist_m / np.diff(bins_m) / 117.0 # per deg^2

        m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
        err_dNdm = np.sqrt(hist_m) / np.diff(bins_m)

        data = np.loadtxt(f'../{filter_name}.txt', delimiter=',')
        plt.figure()
        plt.errorbar(m_cent, dNdm, yerr=err_dNdm, fmt='o', ms=3, label='This work')
        plt.plot(data[:,0], data[:,1], label='observation', color='red')
        plt.yscale('log'); plt.xlabel(f'AB magnitude ({filter_name})'); plt.ylabel('dN/dm [deg$^{-2}$ mag$^{-1}$]')
        plt.legend(); plt.tight_layout()
        plt.savefig(f"/mnt/data_cat3/yuka/output/LF_{filter_name}_SIDES.png")
        plt.savefig(f"/mnt/data_cat3/yuka/output/LF_{filter_name}_SIDES.pdf")
        
        return m_ab
    
    results = [Get_LF(f) for f in filter_names]
    
def Volume(zmin, zmax):
    solid_angle = (1.0 * u.deg**2).to(u.sr)
    # 赤方偏移範囲
    z1 = zmin
    z2 = zmax

    # comoving volumeの差（立体角で割ってない、全立体角の体積）
    vol1 = cosmo.comoving_volume(z1)
    vol2 = cosmo.comoving_volume(z2)

    # 差分を計算（全空間での体積の差）
    volume_total = (vol2 - vol1)

    # 求めたい体積 = 全体積 × (対象領域の立体角 / 4π)
    fraction = (solid_angle / (4 * np.pi * u.sr))
    volume_partial = volume_total * fraction
    return volume_partial.value
    
def test_EL(name):
    file = '/mnt/data_cat3/yuka/data/SIDES/' + f'pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits'
    photoz, Mstar, cat_SFR, Mhalo = Get_Sides(name)

    f_sample = "/mnt/data_cat3/yuka/repository/LIM_mock/theta_samples_100000_seed1_mag.npz"
    data = np.load(f_sample, allow_pickle=True)
    samples = data["samples"] #100000このsampleに対するSED parameter
    logM_samples = data["logM_samples"] #remaining stellar mass (log M_sun)
    #SFR = data['SFR']
    Ha = data["Ha"] #remaining stellar mass (log M_sun)
    Oiii = data['Oiii']
    Pa = data['Pa']
    log10SFR_bins =  data['log10SFR_bins']
    t_edge_bins = data['t_edge_bins']
    redshift = samples[:,-1]
    
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

    #M_array = Get_Mstar(redshift, logM_samples, photoz, np.log10(Mstar))
    z_array = photoz
    M_array = np.log10(Mstar)
    
    SFR = 10**log10SFR_bins[:,0] #mean SFR for the last 30Myr
    print(f'number of galaxies {len(z_array)}')
    sfr = sample_ELs(z_array, M_array, SFR, bin_to_indices, z_edges, logM_edges)
    Ha = GetEL(sfr, z_array)
    log_Ha = np.log10(Ha)
    
    def Get_LF(redshift, Ha, zmin, zmax):
        log_ha = Ha
        volume = Volume(zmin, zmax)
        Lbin = np.linspace(39.5, 44.5, 40)
        hist, bin = np.histogram(log_ha[(redshift>zmin) & (redshift<zmax)], bins=Lbin)
        density = hist/volume/(bin[1:] - bin[:-1])
        bin_center = (bin[1:] + bin[:-1])/2.0
        plt.figure()
        plt.scatter(bin_center, np.log10(density))
        plt.xlim(40.5, 43.5)
        plt.ylim(-5.5, -0.5)
        
    def Get_SFR(redshift, sfr, cat_SFR, zmin, zmax):
        volume = Volume(zmin, zmax)
        Lbin = 10**np.linspace(-4.0, 3.0, 40)
        hist, bin = np.histogram(sfr[(redshift>zmin) & (redshift<zmax)], bins=Lbin)
        density = hist/volume/(bin[1:] - bin[:-1])
        bin_center = (bin[1:] + bin[:-1])/2.0
        plt.figure()
        plt.scatter(bin_center, np.log10(density), color='red')
        
        hist, bin = np.histogram(cat_SFR[(redshift>zmin) & (redshift<zmax)], bins=Lbin)
        density = hist/volume/(bin[1:] - bin[:-1])
        bin_center = (bin[1:] + bin[:-1])/2.0
        plt.scatter(bin_center, np.log10(density), color='green')
        plt.xlim(1e-4, 1e3)
        plt.ylim(-6.0,2.0)
        plt.xscale('log')
        
    Get_LF(z_array, log_Ha, 0.35, 0.45)
    Get_LF(z_array, log_Ha, 0.75, 0.95)
    Get_LF(z_array, log_Ha, 1.35, 1.55)
    Get_LF(z_array, log_Ha, 1.7, 2.8)
    
    Get_SFR(z_array, sfr, cat_SFR, 0.35, 0.45)
    Get_SFR(z_array, sfr, cat_SFR, 0.75, 0.95)
    Get_SFR(z_array, sfr, cat_SFR, 1.35, 1.55)
    Get_SFR(z_array, sfr, cat_SFR, 1.7, 2.8)
    
    #plt.figure()
    #plt.scatter(Mhalo[(photoz>1.8)&(photoz<2.6)], sfr[(photoz>1.8)&(photoz<2.6)], s=1)
    #plt.xlim(1e8, 2e13)
    #plt.ylim(-8, 3)
    #plt.xscale("log")
    #plt.xlabel(r'M [M$_\odot$]')
    #plt.ylabel(r'log SFR [M$_\odot$yr$^$-1$]')
    
    #plt.figure()
    #plt.scatter(Mhalo[(photoz>1.8)&(photoz<2.6)], np.log10(cat_SFR)[(photoz>1.8)&(photoz<2.6)], s=1)
    #plt.xlim(1e8, 2e13)
    #plt.ylim(-8, 3)
    #plt.xscale("log")
    #plt.xlabel(r'M [M$_\odot$]')
    #plt.ylabel(r'log SFR [M$_\odot$yr$^$-1$]')
    
def test_clustering(names, redshift_name='084'):
    with Pool(processes=20) as pool:
        results = pool.map(Get_Sides, names)
    #file = '/mnt/data_cat3/yuka/data/SIDES/' + f'pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits'
    #photoz, Mstar, SFR, Mhalo = Get_Sides(file)
    ra, dec, photoz, Mstar, cat_SFR, Mhalo = zip(*results)
    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    photoz = np.concatenate(photoz)
    Mstar = np.concatenate(Mstar)
    cat_SFR = np.concatenate(cat_SFR)
    Mhalo = np.concatenate(Mhalo)
    
    selection = (photoz < 1.2)
    ra = ra[selection]
    dec = dec[selection]
    photoz = photoz[selection]
    Mstar = Mstar[selection]
    cat_SFR = cat_SFR[selection]
    Mhalo= Mhalo[selection]
    print(f'number of total galaxies: {len(photoz)}')
    #file = '/mnt/data_cat3/yuka/data/SIDES/' + f'pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits'
    #ra, dec, photoz, Mstar, cat_SFR, Mhalo = Get_Sides(name)

    f_sample = "/mnt/data_cat3/yuka/repository/LIM_mock/theta_samples_100000_seed1_mag.npz"
    data = np.load(f_sample, allow_pickle=True)
    samples = data["samples"] #100000このsampleに対するSED parameter
    logM_samples = data["logM_samples"] #remaining stellar mass (log M_sun)
    #SFR = data['SFR']
    #Ha = data["Ha"] #remaining stellar mass (log M_sun)
    #Oiii = data['Oiii']
    #Pa = data['Pa']
    log10SFR_bins =  data['log10SFR_bins']
    t_edge_bins = data['t_edge_bins']
    redshift = samples[:,-1]
    
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

    #M_array = Get_Mstar(redshift, logM_samples, photoz, np.log10(Mstar))
    z_array = photoz
    M_array = np.log10(Mstar)
    
    mask = (photoz>0.83)&(photoz<0.85)
    print("after z selection:", np.sum(mask))
    
    #SFR = 10**(log10SFR_bins[:,0])*0.3 + 10**(log10SFR_bins[:,1])*0.7 #mean SFR for the last 100Myr
    #SFR = cat_SFR
    print(f'number of galaxies {len(z_array)}')
    #sfr = sample_ELs(z_array, M_array,SFR, bin_to_indices, z_edges, logM_edges)
    sfr = cat_SFR
    #####################################################################
    M = Mhalo[mask] / 0.7
    S = sfr[mask] / 0.7
    bins = np.logspace(10, 13, 20)
    mean_sfr, edges, _ = binned_statistic(M, S, statistic='mean', bins=bins)
    bin_centers = 0.5 * (edges[1:] + edges[:-1])

    plt.figure()
    plt.scatter(M, S, s=0.01, alpha=0.1)
    plt.plot(bin_centers, mean_sfr, color='red', lw=2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e9, 1e13)
    plt.ylim(1e-2, 1e3)
    plt.show()
    
    M = Mhalo[mask] / 0.7
    S = cat_SFR[mask] / 0.7
    bins = np.logspace(10, 13, 20)
    mean_sfr, edges, _ = binned_statistic(M, S, statistic='mean', bins=bins)
    bin_centers = 0.5 * (edges[1:] + edges[:-1])

    plt.figure()
    plt.scatter(M, S, s=0.01, alpha=0.1)
    plt.plot(bin_centers, mean_sfr, color='green', lw=2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e9, 1e13)
    plt.ylim(1e-2, 1e3)
    plt.show()
    ####################################################################
    Oiii = GetEL_nb(sfr[mask], z_array[mask], 0.83, 0.83, area=81.0, name='OIII')
    Oiii_flux = Oiii / (4.0 * np.pi * cosmo.luminosity_distance(z_array[mask]).to(u.cm)**2).value #erg/s/cm2

    selection = np.loadtxt(f'/mnt/data_cat3/yuka/repository/LIM_mock/Hizelts/Hizelts_z{redshift_name}.txt', delimiter=',')
    #f = interp1d(selection[:,0], selection[:,1], bounds_error=False, fill_value=(selection[:,1][0], selection[:,1][-1]))
    f = interp1d(selection[:,0], selection[:,1], bounds_error=False, fill_value=(selection[:,1][0], 0.0))
    fraction = f(np.log10(Oiii_flux))
    
    #plt.figure()
    #plt.hist(np.log10(Oiii_flux[mask]), bins=50)
    
    random = np.random.random(len(Oiii_flux)) - fraction
    selected_Oiii = Oiii_flux[random<0]
    print(f'number of selected galaxies: {len(selected_Oiii)}')
    selected_ra = ra[mask][random<0]
    selected_ra -= np.min(selected_ra)
    
    selected_dec = dec[mask][random<0]
    selected_dec -= np.min(selected_dec)
    
    ra_rand = np.random.uniform(0, np.max(selected_ra), 50000)
    dec_rand = np.random.uniform(0, np.max(selected_dec), 50000)
    
    plt.figure(figsize=(6,6))
    plt.scatter(selected_ra, selected_dec, s=1)
    plt.xlabel("RA [deg]")
    plt.ylabel("DEC [deg]")
    
    bins_arcsec = 10**np.linspace(0, 3.5, 20)/3600.0
    
    theta_mid, w = wtheta_landy_szalay(selected_ra, selected_dec, ra_rand, dec_rand, bins_arcsec, nthreads=8)
    plt.figure()
    plt.scatter(theta_mid*3600.0, w, marker='o', color='blue', label='This work')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(2, 3200)
    plt.ylim(1e-4, 20)
    plt.xlabel("theta [arcsec]")
    plt.ylabel("w(theta)")
    
    data = np.loadtxt(f'/mnt/data_cat3/yuka/repository/LIM_mock/Hizelts/z{redshift_name}.txt', delimiter=',')
    plt.scatter(data[:,0], 10**data[:,1], color='red', label='observation')
    plt.legend(); plt.tight_layout()
    #plt.savefig(f"/mnt/data_cat3/yuka/output/clustering_SIDES.png")
    #plt.savefig(f"/mnt/data_cat3/yuka/output/clustering_SIDES.pdf")
    
    
    
    
    
    

if __name__ == "__main__":
    main()