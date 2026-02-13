import numpy as np
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

###################Get stellar mass from abundance matching#########################
def Volume(zmin, zmax, area_deg2=1.0):
    omega = (area_deg2 * u.deg**2).to(u.sr)

    # 赤方偏移範囲
    z1 = zmin
    z2 = zmax

    # comoving volumeの差（立体角で割ってない、全立体角の体積）
    vol1 = cosmo.comoving_volume(z1)
    vol2 = cosmo.comoving_volume(z2)

    # 差分を計算（全空間での体積の差）
    volume_total = (vol2 - vol1)

    # 求めたい体積 = 全体積 × (対象領域の立体角 / 4π)
    fraction = (omega / (4 * np.pi * u.sr))
    volume_partial = volume_total * fraction
    return volume_partial.value

def MassFunc(logM, redshift, zmin, zmax):
    volume = Volume(zmin, zmax, 117.0)
    massbin = np.linspace(10, 15, 40)
    hist, bin = np.histogram(logM[(redshift>zmin) & (redshift<zmax)], bins=massbin)
    density = hist/volume/(bin[1:] - bin[:-1])
    mass = 10**((bin[1:]+bin[:-1])/2.0)
    plt.figure()
    plt.plot(np.log10(mass), np.log10(density), color='red')
    plt.xlim(10,15)
    plt.ylim(-6, -1)
    plt.xlabel('log halo mass [Msun]')
    plt.ylabel('log dN/dlogMstar [Mpc^-3 dex^-1]')
    plt.grid()
    return mass, density

def dN_dlogMstar(Npts_mass, Mstargrid, i):
    # Load evolution of mass function
    (redshift, logM, logPhi1, alpha1, logPhi2, alpha2) = np.loadtxt('../params/SMF.par', unpack=True)

    Mstar_z = 10**logM[i]
    Phi1_z = 10**logPhi1[i]
    Phi2_z = 10**logPhi2[i]
    
    alpha1_z = alpha1[i]
    alpha2_z = alpha2[i]
    
    print(f'logM_*:{np.log10(Mstar_z)} log_phi_1:{np.log10(Phi1_z)} alpha_1:{alpha1_z} alpha_2:{alpha2_z} log_phi_2:{np.log10(Phi2_z)}')
    
    # Calculate the stellar mass function phi for the given z
    term2 = 0.0
    if not (logPhi2[i] == 0 and alpha2[i] == 0):
        term2 = Phi2_z*(Mstargrid/Mstar_z)**alpha2_z

    phi = np.exp(-Mstargrid/Mstar_z) * (Phi1_z*(Mstargrid/Mstar_z)**alpha1_z + term2)
    phi = phi / Mstar_z * Mstargrid * np.log(10)
    
    plt.figure()
    plt.plot(np.log10(Mstargrid), np.log10(phi))
    plt.xlim(7,12)
    plt.ylim(-6, -1)
    plt.xlabel('log Mstar [Msun]')
    plt.ylabel('log dN/dlogMstar [Mpc^-3 dex^-1]')
    plt.grid()
    
    return phi

def Cum_func(L, Psi):
    dlogL = np.log10(L[1]) - np.log10(L[0])
    
    invPsi = Psi[::-1]
    
    sum = np.cumsum(invPsi * dlogL)
    sum = sum[::-1]
    return sum

def AbundanceMatch(logM, redshift, zmin, zmax, i):
    z_mean = (zmin + zmax)/2.0
    
    M, dndlogM = MassFunc(logM, redshift, zmin, zmax)
    Phi_M = Cum_func(M, dndlogM) #zminからzmaxまでの累積数密度

    #get stellar mass function at z_mean
    Npts_mass = (13 - 6) * 50
    Mstargrid = 10.**(6 + (13 - 6) / Npts_mass * np.arange(0, Npts_mass + 1))
    dlogMstargrid = 1. / 50
    dn_dlogmstar = dN_dlogMstar(Npts_mass, Mstargrid, i)
    Phi_Mstar = Cum_func(Mstargrid, dn_dlogmstar)

    #abundance matching
    interp_Mstar_of_Phi =interp1d(Phi_Mstar[::-1], Mstargrid[::-1],bounds_error=False, fill_value='extrapolate')
    Mstar_match = interp_Mstar_of_Phi(Phi_M)
    interp_M_to_Mstar = interp1d(M, Mstar_match,bounds_error=False, fill_value='extrapolate')
    mass = 10**(logM[(redshift>zmin) & (redshift<zmax)])
    Mstar = interp_M_to_Mstar(mass)

    scatter = 0.2
    logMstar_scattered = np.log10(Mstar) + np.random.normal(0, scatter , size=Mstar.shape)
    Mstar_scattered = 10**logMstar_scattered
    return Mstar_scattered

def GetStellarMass(logM, redshift):
    redshift_bin = [0.2, 0.5, 0.8, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5, 6.5]
    Mstar = np.zeros_like(logM)

    for i, z in enumerate(redshift_bin):
        if i==0:
            zmin = 0.0
        else:
            zmin = redshift_bin[i-1]
        zmax = redshift_bin[i]
        mask = (redshift>zmin) & (redshift<zmax)
        Mstar[mask] = AbundanceMatch(logM, redshift, zmin, zmax, i)
    return Mstar

def dN_dlogSFR(logSFR, redshift, zmin, zmax, area):
    volume = Volume(zmin, zmax, area)
    SFRbin = np.linspace(-2.0, 2.5, 40)
    hist, bin = np.histogram(logSFR[(redshift>zmin) & (redshift<zmax)], bins = SFRbin)
    density = hist/volume/(bin[1:] - bin[:-1])
    mass = 10**((bin[1:]+bin[:-1])/2.0)
    return mass, density

def AbundanceMatch_EL(logSFR, redshift, zmin, zmax, i, area, name ='Halpha'):
    z_mean = (zmin + zmax)/2.0
    SFR = 10**logSFR
    
    SFR_bin, dndlogSFR = dN_dlogSFR(logSFR, redshift, zmin, zmax, area)
    Phi_SFR = Cum_func(SFR_bin, dndlogSFR) #zminからzmaxまでの累積数密度

    #get luminosity function at z_mean
    logL_bin = np.linspace(39.5, 44.5, 40)
    L_bin = 10**logL_bin
    (z_param, logPsi_star, logL_star, alpha, A) = np.loadtxt(f'../params/{name}.par', unpack=True)
    z_param_i = z_param[i]
    Psi = 10**logPsi_star[i]
    L_star = 10**logL_star[i]
    alpha_i = alpha[i]
    A_i = A[i]
    
    # Calculate the stellar mass function phi for the given z
    phi = np.exp(-L_bin/L_star) * (Psi*(L_bin/L_star)**alpha_i)
    phi = phi / L_star * L_bin * np.log(10)
    Phi_L = Cum_func(L_bin, phi)

    #abundance matching
    interp_Lstar_of_Phi =interp1d(Phi_L[::-1], L_bin[::-1],bounds_error=False, fill_value='extrapolate')
    L_SFR_match = interp_Lstar_of_Phi(Phi_SFR)
    interp_SFR_to_L = interp1d(SFR_bin, L_SFR_match, bounds_error=False, fill_value='extrapolate')
    L = interp_SFR_to_L(SFR)

    scatter = 0.3
    logL_scattered = np.log10(L) + np.random.normal(0, scatter , size=L.shape)
    L_scattered = 10**logL_scattered*10**(-0.4*A_i)  #convert to observed luminosity
    return L_scattered #erg/s

def GetEL(SFR, redshift, area, name='Halpha'):
    if name=='Halpha':
        #redshift_bin = [0.16, 0.32, 0.62, 1.15, 1.85, 3.34, 4.88, 5.73, 6.50]
        redshift_bin = [0.32, 0.62, 1.15, 1.85, 3.34]
    elif name=='OIII':
        redshift_bin = [1.14, 1.87, 2.73]
    L = np.zeros_like(SFR)
    
    for i, z in enumerate(redshift_bin):
        if i==0:
            zmin = 0.01
        else:
            zmin = redshift_bin[i-1]
        zmax = redshift_bin[i]
        mask = (redshift>zmin) & (redshift<zmax)
        
        L[mask] = AbundanceMatch_EL(np.log10(SFR[mask]), redshift[mask], zmin, zmax, i, area=area, name= name)
        #LOiii[mask] = AbundanceMatch_EL(np.log10(SFR[mask]), redshift, zmin, zmax, i, name='OIII')
    return L

def GetEL_nb(SFR, redshift, zmin, zmax, area, name='Halpha'):
    L = np.zeros_like(SFR)
    
    mask = (redshift>zmin) & (redshift<zmax)
        
    L[mask] = AbundanceMatch_EL(np.log10(SFR[mask]), redshift[mask], zmin, zmax, 0, area=area, name= name)
        #LOiii[mask] = AbundanceMatch_EL(np.log10(SFR[mask]), redshift, zmin, zmax, i, name='OIII')
    return L
        
    