import numpy as np
import os

import pandas as pd

from scipy.interpolate import interp1d
from scipy.interpolate import interpn
from scipy.special import erf
from scipy.integrate import cumulative_trapezoid 

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

import ReadPinocchio5 as rp

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

class Catalog(object):
    def __init__(self):
        self.x = None
        self.y = None
        self.redshift_obs = None
        self.redshift_real = None
        self.logm = None
        self.Mstar = None
        self.SFR = None
        self.JLumi = None # [erg/s/Hz] absolute luminosity
        self.quenched = None
        self.Starburst = None

    def loadcatalog(self, path):
        self.x, self.y, self.redshift_obs, self.redshift_real, self.logm = LoadLightCone(path)
        self.Mstar = GetStellarMass(self.logm, self.redshift_real)
        self.SFR, self.quenched, self.Starburst = get_SFR(self.redshift_real, np.log10(self.Mstar))
        self.JLumi = GetJband(self.redshift_real, self.Mstar, self.quenched) # [erg/s/Hz] absolute luminosity

        remove = (self.redshift_real > 5.85) | (self.Mstar <= 0) | (self.SFR < 0) | (self.JLumi <= 0)
        self.x = self.x[~remove]
        self.y = self.y[~remove]
        self.redshift_obs = self.redshift_obs[~remove]
        self.redshift_real = self.redshift_real[~remove]
        self.logm = self.logm[~remove]
        self.Mstar = self.Mstar[~remove]
        self.SFR = self.SFR[~remove]
        self.JLumi = self.JLumi[~remove]
        self.quenched = self.quenched[~remove]
        self.Starburst = self.Starburst[~remove]


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

def LoadLightCone(path):
    myplc = rp.plc(path)
           
    logm = np.log10( myplc.data["Mass"] ) # [Msun/h]
    theta = myplc.data["theta"] 
    phi = myplc.data["phi"]

    redshift_obs = myplc.data["obsz"]
    redshift_real = myplc.data["truez"]

    theta = ( 90. - theta )  # [deg]
    pos_x = theta * np.cos( phi * np.pi / 180. ) # [deg]
    pos_y = theta * np.sin( phi * np.pi / 180. ) # [deg]

    return pos_x[redshift_real<5.85], pos_y[redshift_real<5.85], redshift_obs[redshift_real<5.85], redshift_real[redshift_real<5.85], logm[redshift_real<5.85]

###################Get stellar mass from abundance matching#########################
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
    massbin = np.linspace(11, 15, 40)
    hist, bin = np.histogram(logM[(redshift>zmin) & (redshift<zmax)], bins=massbin)
    density = hist/volume/(bin[1:] - bin[:-1])
    mass = 10**((bin[1:]+bin[:-1])/2.0)
    return mass, density

def dN_dlogMstar(Npts_mass, Mstargrid, redshift):
    # Load evolution of mass function
    (zmean, logMknee, phiknee1, phiknee2, alpha1, alpha2) = np.loadtxt('../params/Bethermin.par', unpack=True)

    z=redshift
    # Interpolate for the given redshift z
    Mknee_z = 10.**np.interp(1. + z, 1. + zmean, logMknee)
    phiknee1_z = np.interp(1. + z, 1. + zmean, phiknee1)
    alpha1_z = np.interp(1. + z, 1. + zmean, alpha1)
    alpha2_z = np.interp(1. + z, 1. + zmean, alpha2)

    phiknee2_z = np.exp(np.interp(np.log(1. + z), np.log(1. + zmean), np.log(phiknee2)))
    
    print(f'logM_*:{np.log10(Mknee_z)} log_phi_1:{np.log10(phiknee1_z)} alpha_1:{alpha1_z} alpha_2:{alpha2_z} log_phi_2:{np.log10(phiknee2_z)}')
    
    # Calculate the stellar mass function phi for the given z
    phi = np.exp(-Mstargrid / Mknee_z) * (phiknee1_z * (Mstargrid / Mknee_z)**alpha1_z + phiknee2_z * (Mstargrid / Mknee_z)**alpha2_z) / Mknee_z* Mstargrid * np.log(10)
    
    return phi

def Cum_func(L, Psi):
    dlogL = np.log10(L[1]) - np.log10(L[0])
    
    invPsi = Psi[::-1]
    
    sum = np.cumsum(invPsi * dlogL)
    sum = sum[::-1]
    return sum

def AbundanceMatch(logM, redshift, zmin, zmax, faint=0.0):
    z_mean = (zmin + zmax)/2.0
    
    M, dndlogM = MassFunc(logM, redshift, zmin, zmax)
    Phi_M = Cum_func(M, dndlogM) #zminからzmaxまでの累積数密度

    #get stellar mass function at z_mean
    Npts_mass = (13 - 6) * 50
    Mstargrid = 10.**(6 + (13 - 6) / Npts_mass * np.arange(0, Npts_mass + 1))
    dlogMstargrid = 1. / 50
    dn_dlogmstar = dN_dlogMstar(Npts_mass, Mstargrid, z_mean)
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
    redshift_bin = [0.2, 0.44,  0.64, 0.84, 1.1, 1.45, 1.7, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0]
    Mstar = np.zeros_like(logM)

    for i, z in enumerate(redshift_bin):
        if i==0:
            continue
        zmin = redshift_bin[i-1]
        zmax = redshift_bin[i]
        mask = (redshift>zmin) & (redshift<zmax)
        Mstar[mask] = AbundanceMatch(logM, redshift, zmin, zmax)
    return Mstar


##################Get SFR from Bethermin et.al. #################################
def get_SFR(redshift, log_smass):
    full_path = '../params/SIDES.par'
    params = load_params(full_path)
    print('Generate the star-formation properties...')

    print('Draw quenched galaxies...')

    Ngal = len(log_smass)

    #Draw randomly which galaxies are quenched using the recipe from Bethermin+17
    Mtz = params['Mt0'] + params['alpha1'] * redshift + params['alpha2'] * redshift**2
    sigmaz =  params['sigma0'] +  params['beta1'] * redshift +  params['beta2'] * redshift**2
    qfrac0z = params['qfrac0'] * (1.+redshift)**params['gamma']
    
    Prob_SF = (1.-qfrac0z) * 0.5 * (1. - erf( ( log_smass - Mtz) /sigmaz ) )

    Xuni = np.random.rand(Ngal)

    qflag = Xuni > Prob_SF

    #Generate SFR for non-quenched objects

    print('Generate SFRs...')
    Starburst = np.zeros_like(log_smass, dtype=bool)
    index_SF = np.where(qflag == False)

    m_SF = np.array(log_smass[index_SF[0]]+ np.log10(params['Chab2Salp'] / 1.e9 ))
    z_SF = np.asarray(redshift)[index_SF[0]]
    r = np.log10(1.+z_SF)
    expr = np.maximum(m_SF - params['m1'] - params['a2'] * r, 0.)

    logSFRms_SF = m_SF - params['m0'] + params['a0'] * r - params['a1'] * expr**2 - np.log10(params['Chab2Salp'])

    logSFRms_SF += params['corr_zmean_lowzcorr'] * (params['zmax_lowzcorr'] - np.minimum(z_SF, params['zmax_lowzcorr'])) / (params['zmax_lowzcorr'] - params['zmean_lowzcorr'])

    Psb = params['Psb_hz'] + params['slope_Psb'] * (params['z_Psb_knee'] - np.minimum(z_SF, params['z_Psb_knee']))

    Xuni = np.random.rand(np.size(Psb))

    issb = (Xuni < Psb )
    Starburst[index_SF[0][issb]] = True

    SFR_SF = 10. ** ( logSFRms_SF + params['sigma_MS'] * np.random.randn(np.size(logSFRms_SF))
                      + params['logx0'] + issb * (params['logBsb'] - params['logx0']) )

    print('Deal with SFR drawn initially above the SFR limit...')

    too_high_SFRs = np.where( SFR_SF > params['SFR_max'])
    while np.size(too_high_SFRs) > 0:
        SFR_SF[too_high_SFRs[0]] = 10. ** ( logSFRms_SF[too_high_SFRs] + params['sigma_MS'] *
                                         np.random.randn(np.size(too_high_SFRs))
                                         + params['logx0'] + issb[too_high_SFRs] * (params['logBsb'] - params['logx0']) )
        
        too_high_SFRs = np.where( SFR_SF > params['SFR_max'])

    print('Store the results...')
    SFR = np.zeros(Ngal)
    SFR[index_SF[0]] = SFR_SF
    return SFR, qflag, Starburst

def load_params(path, force_pysides_path = ''):

    file = open(path)

    params = {}
    for line in file:
        line = line.strip()
        if not line.startswith("#"):
            no_comment = line.split('#')[0]
            key_value = no_comment.split("=")
            if len(key_value) == 2:
                params[key_value[0].strip()] = key_value[1].strip()

    for key in params.keys():
        params[key] = eval(params[key])

    return params



#################################################################################
def GetJband(redshift, SMass, quenched):
    data_quenched = np.loadtxt("/mnt/data_cat3/yuka/repository/LIM_mock/quenched_ml.txt")
    redshift_bin = data_quenched[:,0]
    quenched_ml = data_quenched[:,1]
    quenched_scatter = data_quenched[:,2]

    quenched_med = interp1d(redshift_bin, quenched_ml, kind='linear')
    quenched_scat = interp1d(redshift_bin, quenched_scatter, kind='linear')
    quenched_ml = quenched_med(redshift) + np.random.normal(0, quenched_scat(redshift))

    data_SF = np.loadtxt("/mnt/data_cat3/yuka/repository/LIM_mock/SF_ml.txt")
    redshift_bin = data_SF[:,0]
    SF_ml = data_SF[:,1]
    SF_scatter = data_SF[:,2]

    SF_med = interp1d(redshift_bin, SF_ml, kind='linear')
    SF_scat = interp1d(redshift_bin, SF_scatter, kind='linear')
    SF_ml = SF_med(redshift) + np.random.normal(0, SF_scat(redshift))

    Mass_lumi = np.zeros_like(SMass)
    Mass_lumi[quenched] = quenched_ml[quenched]
    Mass_lumi[~quenched] = SF_ml[~quenched]
    Jband_norm = SMass / Mass_lumi # [Lsun]

    J_sun_mag = 4.55
    J_sun_lumi = 4.4653 * 10**20 * 10**(-0.4*J_sun_mag) # [erg/s/Hz]
    JLumi = Jband_norm * J_sun_lumi # [erg/s/Hz] luminosity at 10pc
    JLumi[JLumi<0] = 0.0
    return JLumi

###############################################################################
import numpy as np

def GetSED():
    """
    Returns
    -------
    nu : 1D ndarray
        Frequency grid [Hz], ascending.
    Lnu_rel : dict of str->1D ndarray
        Relative SEDs in Lν units for each template on `nu` (same length).
        Keys: 'ELL', 'SD', 'SB'
    LnuJ_template : dict of str->float
        J-band *band-averaged* Lν for each template (relative units).
        Keys: 'ELL', 'SD', 'SB'
    """
    # --- Constants ---
    c_A_per_s = 2.99792458e18  # speed of light [Angstrom/s]

    # --- Load J-band throughput (λ in μm? -> convert to Å) ---
    J_path = '/mnt/data_cat3/yuka/repository/LIM_mock/Jband_throughput.txt'
    wJ, tJ = np.loadtxt(J_path, unpack=True)
    # ここはあなたのファイル仕様に合わせて調整：
    # もし Jband_throughput の波長が μm なら ×1e4、もし Å ならこの行は外してください
    wave_J = wJ * 1e4   # [Å]
    T_J = tJ.copy()

    # --- Helper to load a SWIRE template and convert to relative Lν(ν) ---
    def load_template_as_Lnu_rel(path):
        """
        Parameters
        ----------
        path : str
            SWIRE SED file with columns: wavelength [Å], Fλ [erg s^-1 cm^-2 Å^-1] (relative/normalized)

        Returns
        -------
        nu       : 1D ndarray [Hz], ascending
        Lnu_rel  : 1D ndarray (relative units)
        """
        wl, Flam = np.loadtxt(path, unpack=True)  # wl in Å, Fλ in erg/s/cm^2/Å (relative)
        # Convert to Lν (relative): Lν ∝ λ^2 Fλ / c
        # 相対スケールなので定数因子は不要。比/平均で打ち消される。
        Lnu_rel = Flam * (wl**2) / c_A_per_s  # ~ erg/s/cm^2/Hz (relative)
        # Make frequency grid (ascending)
        nu = c_A_per_s / wl                   # [Hz], this is descending if wl ascending
        order = np.argsort(nu)                # ensure ascending
        nu = nu[order]
        Lnu_rel = Lnu_rel[order]
        return nu, Lnu_rel

    # --- Load three templates (正しいパスを使う！) ---
    path_ELL = "/mnt/data_cat3/yuka/data/LIM_mock/SWIRE/Ell5_template_norm.sed"
    path_SD  = "/mnt/data_cat3/yuka/data/LIM_mock/SWIRE/Sd_template_norm.sed"
    path_SB  = "/mnt/data_cat3/yuka/data/LIM_mock/SWIRE/M82_template_norm.sed"

    nu_ell, Lnu_rel_ell = load_template_as_Lnu_rel(path_ELL) # [um], [relative]
    nu_sd,  Lnu_rel_sd  = load_template_as_Lnu_rel(path_SD) # [um], [relative]
    nu_sb,  Lnu_rel_sb  = load_template_as_Lnu_rel(path_SB) # [um], [relative]

    # --- Make common frequency grid covering all templates ---
    nu_min = min(nu_ell.min(), nu_sd.min(), nu_sb.min())
    nu_max = max(nu_ell.max(), nu_sd.max(), nu_sb.max())
    N_nu = 5000
    nu = np.logspace(np.log10(nu_min), np.log10(nu_max), N_nu)  # [Hz], ascending

    # それぞれを共通gridへ補間 
    Lnu_rel_ell = np.interp(nu, nu_ell, Lnu_rel_ell, left=0.0, right=0.0) 
    Lnu_rel_sd = np.interp(nu, nu_sd, Lnu_rel_sd, left=0.0, right=0.0) 
    Lnu_rel_sb = np.interp(nu, nu_sb, Lnu_rel_sb, left=0.0, right=0.0)
    
    # --- J throughput を frequency へ補間（ascending ν 必須） ---
    # wave_J は Å → frequency に変換
    nu_J = c_A_per_s / wave_J
    orderJ = np.argsort(nu_J)
    nu_J = nu_J[orderJ]
    T_J  = T_J[orderJ]

    # 共通νグリッド上に J throughput を補間
    Tnu_J = np.interp(nu, nu_J, T_J, left=0.0, right=0.0)

    # --- Band-averaged Lν in J for each template (relative units) ---
    # <Lν>_J = ∫ Lν(ν) T(ν) dν / ∫ T(ν) dν
    def band_average_Lnu(Lnu_rel, Tnu):
        num = np.trapz(Lnu_rel * Tnu, nu)
        den = np.trapz(Tnu, nu)
        return num / den if den > 0 else 0.0

    LnuJ_ell = band_average_Lnu(Lnu_rel_ell, Tnu_J)
    LnuJ_sd  = band_average_Lnu(Lnu_rel_sd,  Tnu_J)
    LnuJ_sb  = band_average_Lnu(Lnu_rel_sb,  Tnu_J)

    Lnu_rel = {
        'ELL': Lnu_rel_ell,
        'SD' : Lnu_rel_sd,
        'SB' : Lnu_rel_sb
    }
    LnuJ_template = {
        'ELL': LnuJ_ell,
        'SD' : LnuJ_sd,
        'SB' : LnuJ_sb
    }
    return nu, Lnu_rel, LnuJ_template

##############################################################################

def make_mock(params):
    flist, dflist = frequency(params) # Hz, Hz

    catalog = Catalog()
    input_fname = '/mnt/data_cat3/moriwaki/Pinocchio/my_lightcone/output_large/pinocchio.r00001.plc.out'
    catalog.loadcatalog(input_fname)

    lumi_dist = cosmo.luminosity_distance(np.array(catalog.redshift_real)).to(u.cm)
    lumi_dist = lumi_dist.value
    
    ra = np.array(catalog.x) * 3600.0 #[arcsec]
    dec = np.array(catalog.y) * 3600.0 #[arcsec]
    
    ix = np.floor((ra + params.radius * 3600) / params.resolution).astype(np.int64)
    iy = np.floor((dec + params.radius * 3600) / params.resolution).astype(np.int64)
    print(ix)
    
    ##initialize
    Nx = int(2.0*3600*params.radius / params.resolution)
    Nz = len(flist)
    
    npix = np.array([Nx, Nx, Nz])

    ############################Stellar continuum############################
    valid_xy = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Nx)
    ix_v = ix[valid_xy]
    iy_v = iy[valid_xy]

    continuum = np.zeros([Nx, Nx, Nz], dtype=np.float32)

    nu_obs = flist  # shape (Nν,)
    zobs = np.asarray(catalog.redshift_obs)      # shape (Ngal,)
    nu_rest = (1.0 + zobs)[:, None] * nu_obs[None, :]   # shape (Ngal, Nν)

    mask_SB  = catalog.Starburst.astype(bool)
    mask_Q   = catalog.quenched.astype(bool) & (~mask_SB)
    mask_Sd  = (~mask_Q) & (~mask_SB)
    nu_tpl, Lnu_rel, LnuJ_template = GetSED() # [Hz]
    print(f"ELL:{LnuJ_template['ELL']} SD:{LnuJ_template['SD']}, SB:{LnuJ_template['SB']}")

    # スケール係数（J-band 平均 Lν の比）
    log_scale = np.zeros_like(catalog.JLumi, dtype=float)
    if np.any(mask_Q):
        log_scale[mask_Q]  = np.log10(catalog.JLumi[mask_Q])  - np.log10(LnuJ_template['ELL'])
        plt.hist(log_scale[mask_Q], bins=30, alpha=0.5, label='ELL')
        plt.yscale('log')
    if np.any(mask_Sd):
        log_scale[mask_Sd]  = np.log10(catalog.JLumi[mask_Sd])  - np.log10(LnuJ_template['SD'])
        plt.hist(log_scale[mask_Sd], bins=30, alpha=0.5, label='SD')
        plt.yscale('log')
    if np.any(mask_SB):
        log_scale[mask_SB]  = np.log10(catalog.JLumi[mask_SB])  - np.log10(LnuJ_template['SB'])
        plt.hist(log_scale[mask_SB], bins=30, alpha=0.5, label='SB')
        plt.yscale('log')
    plt.legend()

    def bandavg_group_top_hat(nu_grid, Lnu_grid, nu0_obs, dnu_obs, z):
        """
        テンプレ (nu_grid, Lnu_grid) を、各銀河の z に対して
        観測チャンネル (nu0_obs, dnu_obs) を top-hat として帯域平均 Lν を返す。
        返り値 shape = (N_gal, N_chan)
        """
        # 1) 累積積分 I(ν) を作成（テンプレ1本につき一回）
        # cumulative_trapezoid(x,y): ∫ y dx を返す（len=N-1）。先頭0を付けて長さを揃える
        I = cumulative_trapezoid(Lnu_grid, nu_grid, initial=0.0)  # shape = (N_grid,)

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


    # --- allocate Lν for all galaxies (物理単位) ---
    log_Lnu_all = -np.ones((len(zobs), Nz), dtype=np.float32)

    # 観測周波数 grid: flist (obs), 各銀河の z
    nu_rest = (1.0 + zobs)[:, None] * flist[None, :]
    nu0_obs = flist
    dnu_obs = dflist
    z = zobs

    # Q group (ELL)
    idx = np.where(mask_Q)[0]
    if idx.size:
        Lnu_band = bandavg_group_top_hat(nu_tpl, Lnu_rel['ELL'], nu0_obs, dnu_obs, z[idx])
        log_Lnu_all[idx, :] = (log_scale[idx, None]  + np.log10(Lnu_band)).astype(np.float32)

    # Sd group
    idx = np.where(mask_Sd)[0]
    if idx.size:
        Lnu_band = bandavg_group_top_hat(nu_tpl, Lnu_rel['SD'], nu0_obs, dnu_obs, z[idx])
        log_Lnu_all[idx, :] = (log_scale[idx, None]  + np.log10(Lnu_band)).astype(np.float32)

    # SB group
    idx = np.where(mask_SB)[0]
    if idx.size:
        Lnu_band = bandavg_group_top_hat(nu_tpl, Lnu_rel['SB'], nu0_obs, dnu_obs, z[idx])
        log_Lnu_all[idx, :] = (log_scale[idx, None]  + np.log10(Lnu_band)).astype(np.float32)

    # --- convert to observed fν (cgs) and Jy ---
    # fν = Lν / (4π DL^2) / (1+z)
    log_geom = 2* np.log10(lumi_dist)[:, None] + np.log10(4.0 * np.pi)                    # (Ngal,1)
    log_fnu_cgs = log_Lnu_all -log_geom - np.log10(1.0 + zobs)[:, None]            # erg/s/cm^2/Hz
    fnu_cgs = 10**log_fnu_cgs
    fnu_jy  = (fnu_cgs / Jy).astype(np.float32)/ (params.resolution * arcsec)**2              # Jy/sr

    # 🔴 valid銀河だけ抽出して add.at する
    fnu_valid = fnu_jy[valid_xy, :]                               # (Nvalid, Nz)

    # --- 3D accumulation: (x,y,ν) ---
    # add.at は重複 index の和を正しくとってくれる
    np.add.at(continuum, (ix_v, iy_v, slice(None)), fnu_valid)
    ########################### Ha intensity
    SFR = catalog.SFR
    Ha_intensity = np.zeros([Nx, Nx, Nz], dtype=np.float32)
    
    freq_obs = ha/(1 + np.array(zobs))*1e9 #[Hz]
    iz = np.searchsorted(flist, freq_obs, side='right')
    
    indices = np.array([ix, iy, iz]).T
    valid_xy = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Nx) & (iz >= 0) & (iz < Nz)
    is_SF = SFR > 0
    valid_mask = valid_xy & is_SF
    
    indices_valid = indices[valid_mask]
    lumi_dist_valid = lumi_dist[valid_mask]

    log_lumi = np.log10(SFR[valid_mask]) + np.log10( 8.3e40 )
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi))

    lumi_valid = 10**(log_lumi + log_lumi_scatter)
    flux_valid = lumi_valid / (4. * np.pi * lumi_dist_valid**2) #[erg/s/cm^2]
    intensity_valid = flux_valid / dflist[indices_valid[:,2]] /Jy / (params.resolution * arcsec)**2 #[Jy/sr]
    
    np.add.at(Ha_intensity, (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]), intensity_valid)
    ########################### Pa intensity
    Pa_intensity = np.zeros([Nx, Nx, Nz], dtype=np.float32)
    
    freq_obs = pa/(1 + np.array(zobs))*1e9 #[Hz]
    iz = np.searchsorted(flist, freq_obs, side='right')
    
    indices = np.array([ix, iy, iz]).T
    valid_xy = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Nx) & (iz >= 0) & (iz < Nz)
    valid_mask = valid_xy & is_SF
    
    indices_valid = indices[valid_mask]
    lumi_dist_valid = lumi_dist[valid_mask]

    log_lumi = np.log10(SFR[valid_mask]) + np.log10( 8.3e40 ) - 0.6
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi))

    lumi_valid = 10**(log_lumi + log_lumi_scatter) #[erg/s]
    flux_valid = lumi_valid / (4. * np.pi * lumi_dist_valid**2) #[erg/s/cm^2]
    intensity_valid = flux_valid / dflist[indices_valid[:,2]] /Jy / (params.resolution * arcsec)**2 #[Jy/sr]
    
    np.add.at(Pa_intensity, (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]), intensity_valid)
    
    ########################### [OIII] intensity
    OIII_intensity = np.zeros([Nx, Nx, Nz], dtype=np.float32)
    
    freq_obs = Oiii/(1 + np.array(zobs)) * 1e9 #[Hz]
    iz = np.searchsorted(flist, freq_obs, side='right')
    
    indices = np.array([ix, iy, iz]).T
    valid_xy = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Nx) & (iz >= 0) & (iz < Nz)
    valid_mask = valid_xy & is_SF
    
    indices_valid = indices[valid_mask]
    lumi_dist_valid = lumi_dist[valid_mask]

    log_lumi = np.log10(SFR[valid_mask]) + np.log10( 1.32e41 )
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi))

    lumi_valid = 10**(log_lumi + log_lumi_scatter) #[erg/s]
    flux_valid = lumi_valid / (4. * np.pi * lumi_dist_valid**2) #[erg/s/cm^2]
    intensity_valid = flux_valid / dflist[indices_valid[:,2]] /Jy / (params.resolution * arcsec)**2 #[Jy/sr]
    
    np.add.at(OIII_intensity, (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]), intensity_valid)
    
    total_intensity = Ha_intensity + OIII_intensity + Pa_intensity + continuum
    
    return total_intensity, Ha_intensity, OIII_intensity, Pa_intensity, continuum