import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

from astropy.cosmology import Planck18 as cosmo
from astropy import units as u

####### luminosity function
def Ha_LF(z, lumi):
    alpha = -1.6
    L = 0.45*z + 41.87
    Lstar = 10**L
    f = interp1d([0, 0.4, 0.84, 1.47, 2.23], [-3.15, -2.9, -2.47, -2.61, -2.78], bounds_error=False, fill_value = 'extrapolate')
    psi = f(z)
    Psi = np.log(10)*10**psi*(lumi/Lstar)**alpha * np.exp(-lumi/Lstar)*(lumi/Lstar)
    return Psi

def OIII_LF(z, lumi):
    alpha = -1.6
    f = interp1d([0.84, 1.42, 2.23, 3.24], [41.79, 42.06, 42.66, 42.83], bounds_error=False, fill_value = 'extrapolate')
    L = f(z)
    Lstar = 10**L
    f = interp1d([0.84, 1.42, 2.23, 3.24], [-2.55, -2.61, -3.03, -3.31], bounds_error=False, fill_value = 'extrapolate')
    psi = f(z)
    Psi = np.log(10)*10**psi*(lumi/Lstar)**alpha * np.exp(-lumi/Lstar)*(lumi/Lstar)
    return Psi


def Volume(zmin, zmax):
    theta = 7.5 * u.deg
    solid_angle = 2 * np.pi * (1 - np.cos(theta.to(u.rad))) * u.sr

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

def MassFunc(logM, redshift, zmin, zmax):
    volume = Volume(zmin, zmax)
    massbin = np.linspace(12, 15, 40)
    hist, bin = np.histogram(logM[(redshift>zmin) & (redshift<zmax)], bins=massbin)
    density = hist/volume/(bin[1:] - bin[:-1])
    mass = 10**((bin[1:]+bin[:-1])/2.0)
    return mass, density

def Cum_func(L, Psi):
    dlogL = np.log10(L[1]) - np.log10(L[0])
    
    invPsi = Psi[::-1]
    
    sum = np.cumsum(invPsi * dlogL)
    sum = sum[::-1]
    return sum

def AbundanceMatch(logM, redshift, zmin, zmax):
    z_mean = (zmin + zmax)/2.0
    
    M, dndlogM = MassFunc(logM, redshift, zmin, zmax)
    Phi_M = Cum_func(M, dndlogM)
    
    L = 10**np.linspace(40.5, 44.5, 200)
    LF = Ha_LF(z_mean, L)
    LF_OIII = OIII_LF(z_mean, L)
    Phi_L = Cum_func(L, LF)
    Phi_L_OIII = Cum_func(L, LF_OIII)
    
    interp_L_of_Phi =interp1d(Phi_L[::-1], L[::-1],bounds_error=False, fill_value='extrapolate')
    L_match = interp_L_of_Phi(Phi_M)
    interp_m_to_L = interp1d(M, L_match,bounds_error=False, fill_value='extrapolate')
    
    interp_L_OIII =interp1d(Phi_L_OIII[::-1], L[::-1],bounds_error=False, fill_value='extrapolate')
    L_match_OIII = interp_L_OIII(Phi_M)
    interp_m_to_L_OIII = interp1d(M, L_match_OIII,bounds_error=False, fill_value='extrapolate')
    
    mass = 10**(logM[(redshift>zmin) & (redshift<zmax)])
    luminosity = interp_m_to_L(mass)
    luminosity_OIII = interp_m_to_L_OIII(mass)
    
    scatter = 0.2
    logL_scattered = np.log10(luminosity) + np.random.normal(0, scatter , size=luminosity.shape)
    logL_OIII_scattered = np.log10(luminosity_OIII) + np.random.normal(0, scatter , size=luminosity.shape)
    luminosity_scattered = 10**logL_scattered
    luminosity_OIII_scattered = 10**logL_OIII_scattered
    
    luminosity_scattered = np.nan_to_num(luminosity_scattered)
    luminosity_OIII_scattered = np.nan_to_num(luminosity_OIII_scattered)
    return luminosity_scattered, luminosity_OIII_scattered

