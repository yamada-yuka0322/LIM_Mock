import numpy as np
import matplotlib.pyplot as plt
import torch
import pickle

from astropy.cosmology import Planck18

from pop_cosmos.constants import COSMOS_FILTERS_LATEX
from pop_cosmos.catalogue import CatalogueGenerator
from pop_cosmos.emlines import EmLineEmulator
from pop_cosmos.utils import compute_derived_quantities, compute_mass_remaining

import os
os.environ['SPS_HOME'] = '/mnt/data_cat3/yuka/repository/fsps'

import fsps

from astropy.cosmology import Planck18 as cosmo
from astropy import units as u

from multiprocessing  import Pool


seed = 1
np.random.seed(seed)

Nsample = 100000
#Nsample = 1000

sp = None
sp_cont = None

log_Lsun = 33.58297 # log10 of solar luminosity in erg/s/Hz
c = 2.9979e10 #[cm/s]

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
    
    global sp_cont
    sp_cont = fsps.StellarPopulation(
        zcontinuous=1,
        add_neb_emission=False,
        add_neb_continuum=True,
        add_dust_emission=True,
        redshift_colors=1,
        sfh=3,
        imf_type=1,
        dust_type=0,
    )

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

def get_emission_line(wave, spec, spec_cont):
    """
    wave: wavelength in Å
    spec: spectrum with nebular emission L_sun/Hz
    spec-cont: spectrum without nebular emission L_sun/Hz
    """
    ha = 6565
    Oiii = 5008.5
    Pa = 18748.6
    
    line = spec - spec_cont
    
    nu = c * 1e8 / wave #frequency Hz
    #Haは+-20Åで計算
    ha_range = (wave > ha - 10) & (wave < ha + 10)
    ha_nu = nu[ha_range]
    ha_line = line[ha_range]
    
    order = np.argsort(ha_nu)
    ha_nu = ha_nu[order]
    ha_line = ha_line[order]
    I_ha = np.trapz(ha_line, ha_nu)
    
    Oiii_range = (wave > Oiii - 10) & (wave < Oiii + 10)
    Oiii_nu = nu[Oiii_range]
    Oiii_line = line[Oiii_range]
    
    order = np.argsort(Oiii_nu)
    Oiii_nu = Oiii_nu[order]
    Oiii_line = Oiii_line[order]
    I_Oiii = np.trapz(Oiii_line, Oiii_nu)
    
    Pa_range = (wave > Pa - 750) & (wave < Pa + 750)
    Pa_nu = nu[Pa_range]
    Pa_line = line[Pa_range]
    
    order = np.argsort(Pa_nu)
    Pa_nu = Pa_nu[order]
    Pa_line = Pa_line[order]
    I_Pa = np.trapz(Pa_line, Pa_nu)
    return I_ha, I_Oiii, I_Pa
    

def generate_spectrum(args):
    theta_samples, logM_formed, mw_age, M_frac_bins, t_edge_bins = args
    N_pred, zmeta_pred, sfr1_pred, sfr2_pred, sfr3_pred, sfr4_pred, sfr5_pred, sfr6_pred,  tau2_pred, dust_n_pred, tau1_pred, fAGN_pred, tauAGN_pred, gasZ_pred, gasU_pred, redshift_pred = theta_samples

    Mstar_pred = logM_formed
    
    t_diff = np.diff(t_edge_bins) * 1e9 # time within bin in yr
    SFH_time = t_edge_bins[-1] - t_edge_bins
    SFR = 10**M_frac_bins #SFR M_sun/yr

    tage = mw_age
    SFH_time = SFH_time[::-1]
    SFR = SFR[::-1]
    age_seq = (SFH_time[1::]+SFH_time[:-1:])*0.5

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
    sp.set_tabular_sfh(age_seq, SFR)
    
    sp_cont.params['zred']      = redshift_pred
    sp_cont.params['logzsol']   = zmeta_pred
    sp_cont.params['gas_logz']  = gasZ_pred
    sp_cont.params['gas_logu']  = gasU_pred
    sp_cont.params['dust1']     = tau2_pred * tau1_pred
    sp_cont.params['dust2']     = tau2_pred
    sp_cont.params['dust_index']= dust_n_pred
    sp_cont.params['fagn']      = 10**fAGN_pred
    sp_cont.params['agn_tau']   = 10**tauAGN_pred
    sp_cont.params['tage']      = tage
    sp_cont.set_tabular_sfh(age_seq, SFR)

    remaining_ratio = sp.stellar_mass
    
    if remaining_ratio <= 0:
        print(f'remaining ratio 0 for redshift {redshift_pred} and stellar mass {Mstar_pred}')
        return 0, False, 0, np.zeros(5994), np.zeros(5994), np.zeros(5994),0
    else:
        log_Mstar = np.log10(remaining_ratio)
        mag = sp.get_mags(tage=tage, redshift=redshift_pred, bands=['jwst_f277w'])
        mag = mag[0]
        wave, spec = sp.get_spectrum(tage=tage, peraa=False)
        wave, spec_cont = sp_cont.get_spectrum(tage=tage, peraa=False)
        
        log_lbol = sp.log_lbol

        #log_spec = np.log10(spec) - np.log10(remaining_ratio)  # Lsun/Hz/M_sun SED per solar mass (for remaining stellar mass)
        return log_Mstar, True, mag, wave, spec, spec_cont, log_lbol



def main():
    init_worker_fsps()
    #Load model trained in Thorp et al. (2025) on the CPU
    catalogue_generator = torch.load("/home/yuka/anaconda3/envs/LIM/pop-cosmos/trained_models/catalogueModelT25.pt", weights_only=False)

    #Generate base samples for the diffusion model
    base_noise, base_sigma, base_phi = catalogue_generator.generate_base_samples(Nsample)

    #Generate a representative catalogue from the model
    noisy_fluxes, noisy_magnitudes, noisy_asinh_magnitudes, flux_sigmas, theta_samples, model_fluxes = catalogue_generator(base_noise, base_sigma, base_phi)

    # Selection
    lims = torch.inf * torch.ones(26)
    lims[-2] = 25.0
    lims[-1] = 25.0

    # Move to numpy
    #noisy_magnitudes = noisy_magnitudes.detach().numpy()
    logM_formed, mw_age, log10SFR, log10sSFR, M_frac_bins, log10SFR_bins, t_edge_bins = compute_derived_quantities(theta_samples, use_astropy=True, return_bins=True)
    logM_formed = logM_formed.detach().numpy()
    #logSFR = log10SFR.detach().numpy()
    mw_age = mw_age.detach().numpy()
    log10SFR_bins = log10SFR_bins.detach().numpy() #log SFR M_sun/yr
    t_edge_bins = t_edge_bins.detach().numpy() #Edges of SFH bins in lookback time, with units of Gyr.
    theta_samples = theta_samples.detach().numpy() # (Nsample, num_params)

    worker_iter = zip(theta_samples, logM_formed, mw_age, log10SFR_bins, t_edge_bins)

    with Pool(processes=20) as pool:
        results = pool.map(generate_spectrum, worker_iter)

    logM_samples, flag, mag, wave, spec, spec_cont, log_lbol = zip(*results)
    
    logM_samples = np.array(logM_samples)
    flag = np.asarray(flag).astype(bool).reshape(-1)
    mag  = np.asarray(mag).reshape(-1)
    wave = np.asarray(wave)
    spec = np.asarray(spec)
    spec_cont = np.asarray(spec_cont)
    
    log_lbol = np.asarray(log_lbol)

    mask = flag

    np.savez(
        "../../theta_samples_{:d}_seed{:d}_mag.npz".format(Nsample, seed),
        samples=theta_samples[mask, :],
        logM_samples=logM_samples[mask],
        logM_formed=logM_formed[mask],
        tage=mw_age[mask],
        log10SFR_bins=log10SFR_bins[mask, ...],
        t_edge_bins=t_edge_bins[mask, ...],
        wavelength = wave,
        spectrum = spec[mask, :],
        continuum = spec_cont[mask, :],
        log_lbol = log_lbol[mask],
        )

if __name__ == "__main__":
    main()