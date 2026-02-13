import numpy as np
import os
os.environ['SPS_HOME'] = '/mnt/data_cat3/yuka/repository/fsps'
import fsps

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
from astropy.table import Table

import healpy as hp

import matplotlib.pyplot as plt

import datetime

from multiprocessing  import Pool
from functools import partial

import lim.ReadPinocchio5 as rp

import h5py

import zarr
from numcodecs import Blosc

Oiii = 5.997e5 #[GHz]
ha = 4.856e5 #[GHz]
pa = 1.599e5 #[GHz]
c = 2.9979e10 #[cm/s]
Jy = 1.0e-23
arcsec = 4.848136811094e-6
pc = 3085677857083127000
log_Lsun = 33.58297

sp = None

#crit_df = pd.read_csv('../TNG300_sSFR_classification_criteria.csv')
# スナップショット一覧
#df_z = crit_df[['snapshot_index', 'redshift']].drop_duplicates()
#df_z = df_z.sort_values('redshift').reset_index(drop=True)
#snap_ids = df_z['snapshot_index'].to_numpy()   # shape (Nsnap,)
#snap_z   = df_z['redshift'].to_numpy()         # shape (Nsnap,)

# 質量ビンの edges（全スナップで同じ前提）
#df0 = crit_df[crit_df['snapshot_index'] == snap_ids[0]].copy()
#df0 = df0.sort_values('logM_min')
#edges = np.r_[df0['logM_min'].values, df0['logM_max'].values[-1]]  # shape (Nbin+1,)
#Nbin  = len(edges) - 1
#Nsnap = len(snap_ids)

# low/high の 2D テーブル [Nsnap, Nbin]
#low_table  = np.full((Nsnap, Nbin), np.nan, dtype=float)
#high_table = np.full((Nsnap, Nbin), np.nan, dtype=float)

#for j, sid in enumerate(snap_ids):
    #d = crit_df[crit_df['snapshot_index'] == sid].copy()
    #d = d.sort_values('logM_min')
    #low_table[j, :]  = d['sSFR_low_thr'].to_numpy()
    #high_table[j, :] = d['sSFR_high_thr'].to_numpy()
    
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

class Params(object):
    def __init__(self):
        self.frequency = None
        self.R = None
        self.resolution = None
        self.radius = None
    def Default(self):
        self.frequency = [60000, 400000] #[GHz, GHz]
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

        self.F115W = None
        self.F150W = None
        self.F277W = None
        self.F444W = None
        self.Jmag = None

    def loadcatalogs(self, path, names):
        with Pool(processes=len(names)) as pool:
            results = pool.map(partial(LoadLightCone_SIDES, path), names)
        self.x = np.concatenate([result[0] for result in results])
        self.y = np.concatenate([result[1] for result in results])
        self.redshift_real = np.concatenate([result[2] for result in results])
        self.logm = np.concatenate([result[3] for result in results])
        self.Mstar = np.concatenate([result[4] for result in results])
        self.SFR = np.concatenate([result[5] for result in results])
        self.quenched = np.concatenate([result[6] for result in results])
        self.Starburst = np.concatenate([result[7] for result in results])
        self.redshift_obs = self.redshift_real

        self.JLumi = GetJband(self.redshift_real, self.Mstar, self.quenched, self.Starburst) # [erg/s/Hz] absolute luminosity

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



    def loadcatalog(self, path, name = 'Pinocchio'):
        if name=='Pinocchio':
            self.x, self.y, self.redshift_obs, self.redshift_real, self.logm = LoadLightCone(path)
            self.Mstar = GetStellarMass(self.logm, self.redshift_real)
            self.SFR, self.quenched, self.Starburst = get_SFR(self.redshift_real, np.log10(self.Mstar))
        elif name=='TNG':
            self.x, self.y, self.redshift_real, self.Mstar, self.SFR= LoadLightCone_TNG(path)
            self.redshift_obs = self.redshift_real
            self.quenched, self.Starburst, _ = get_label_array(self.redshift_real, self.Mstar, self.SFR)
            #self.SFR, self.quenched, self.Starburst = get_SFR(self.redshift_real, np.log10(self.Mstar))
        elif name=='JWST':
            self.x, self.y, self.F115W, self.F150W, self.F277W, self.F444W, self.Jmag, self.redshift_real, self.Mstar, self.SFR= LoadJWST(path)
            self.redshift_obs = self.redshift_real
            self.quenched, self.Starburst, _ = get_label_array(self.redshift_real, self.Mstar, self.SFR)

        elif 'SIDES' in name:
            self.x, self.y, self.redshift_real, self.logm, self.Mstar, self.SFR, self.quenched, self.Starburst = LoadLightCone_SIDES(path, name=name)
            self.redshift_obs = self.redshift_real
        
        self.JLumi = GetJband(self.redshift_real, self.Mstar, self.SFR, self.quenched, self.Starburst) # [erg/s/Hz] absolute luminosity

        remove = (self.redshift_real > 5.85) | (self.Mstar <= 0) | (self.SFR < 0) | (self.JLumi <= 0)
        self.x = self.x[~remove]
        self.y = self.y[~remove]
        self.redshift_obs = self.redshift_obs[~remove]
        self.redshift_real = self.redshift_real[~remove]
        if name=='Pinocchio':
            self.logm = self.logm[~remove]
        elif 'SIDES' in name:
            self.logm = self.logm[~remove]
        elif name=='JWST':
            self.F115W = self.F115W[~remove]
            self.F150W = self.F150W[~remove]
            self.F277W = self.F277W[~remove]
            self.F444W = self.F444W[~remove]
            self.Jmag = self.Jmag[~remove]

        self.Mstar = self.Mstar[~remove]
        self.SFR = self.SFR[~remove]
        self.JLumi = self.JLumi[~remove]
        self.quenched = self.quenched[~remove]
        self.Starburst = self.Starburst[~remove]

#def frequency(params):
    #fmin_GHz, fmax_GHz = params.frequency
    #R = params.R  # ← 分光のRにリネーム推奨
    #fmin = float(fmin_GHz) * 1e9  # Hz
    #fmax = float(fmax_GHz) * 1e9  # Hz

    # チャンネル数（上記の式）
    #N = int(np.floor(np.log(fmax/fmin) / np.log(1.0 + 1.0/R))) + 1

    # 等比列で中心周波数を生成
    #ratio = 1.0 + 1.0/R
    #flist = fmin * ratio ** np.arange(N, dtype=np.float64)  # Hz

    # Δf = f/R （定数Rなので各チャンネルごとに比例）
    #dflist = flist / R

    # float32に落とすなら最後に
    #return flist.astype(np.float32), dflist.astype(np.float32)

def frequency(params):
    fmin_GHz, fmax_GHz = params.frequency
    fmin = float(fmin_GHz) * 1e9  # Hz
    fmax = float(fmax_GHz) * 1e9  # Hz
    
    flist = []
    dflist = []
    
    fnow = fmin
    while(fnow < fmax):
        wavelength = c/fnow * 1e4 #um
        if ((wavelength <= 2.42)):
            R = 41
        elif((2.42 < wavelength)&(wavelength <= 3.83)):
            R = 35
        elif((3.83 < wavelength)&(wavelength < 4.42)):
            R=110
        elif(4.42 <= wavelength):
            R=130
        flist.append(fnow)
        df = fnow / R
        dflist.append(df)
        
        fnow += df
    return np.array(flist), np.array(dflist)
            

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

def LoadLightCone_TNG(path):
    data = np.loadtxt(path)

    theta = data[:,0]/3600.0 #[deg]
    phi = data[:,1]/3600.0
    redshift = data[:,2]

    log_mstar = data[:,7]
    log_SFR = data[:,6]
    mstar = 10**log_mstar/0.6774 #[Msun]
    SFR = 10**log_SFR #[Msun/yr]
    return theta[redshift<5.85], phi[redshift<5.85], redshift[redshift<5.85], mstar[redshift<5.85], SFR[redshift<5.85]

def LoadLightCone_SIDES(path, name):
    file = path + f'pySIDES_from_uchuu_tile_{name[5:-1]}_{name[-1]}.fits'

    ra, dec, redshift, Mhalo, Mstar, SFR, quenched, Starburst = Get_Sides(file)

    logm = np.log10(Mhalo)
    return ra[redshift<5.85], dec[redshift<5.85], redshift[redshift<5.85], logm[redshift<5.85], Mstar[redshift<5.85], SFR[redshift<5.85], quenched[redshift<5.85], Starburst[redshift<5.85]

def LoadJWST(path):
    file = path + 'COSMOSWeb_mastercatalog_v1_photom_primary.fits'
    hdu = fits.open(file)
    data = hdu[1].data
    ra = data['ra'] #[deg]
    dec = data['dec'] #[deg]
    f115w = data['mag_model_f115w']
    f150w = data['mag_model_f150w']
    f277w = data['mag_model_f277w']
    f444w = data['mag_model_f444w']
    Jmag = data['mag_model_uvista-j']
    hdu.close()

    file = path + 'COSMOSWeb_mastercatalog_v1_lephare.fits'
    hdu = fits.open(file)
    data = hdu[1].data
    redshift = data['zfinal']
    logMstar = data['mass_minchi2']
    Mstar = 10**logMstar #[Msun]
    logSFR = data['sfr_minchi2'] #[Msun/yr]
    SFR = 10**logSFR #[Msun/yr]
    type = data['type']
    hdu.close()

    selection = (redshift > 0) & (redshift<5.85) & (type == 0)

    return ra[selection], dec[selection], f115w[selection], f150w[selection], f277w[selection], f444w[selection], Jmag[selection], redshift[selection], Mstar[selection], SFR[selection]

def Get_Sides(file):
    hdu = fits.open(file)
    data = hdu[1].data

    ra = data['ra'] #[deg]
    dec = data['dec'] #[deg]
    redshift = data['redshift']
    Mhalo = np.log10(data['Mhalo']) #[Msun]
    Mstar = data['Mstar'] #[Msun]
    quenched = data['qflag'] #[bool]
    Starburst = data['issb'] #[bool]
    SFR = data['SFR'] #[Msun/yr]

    return ra[redshift<5.85], dec[redshift<5.85], redshift[redshift<5.85], Mhalo[redshift<5.85], Mstar[redshift<5.85], SFR[redshift<5.85], quenched[redshift<5.85], Starburst[redshift<5.85]


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

def get_label_array(redshift_arr, SMass_arr, SFR_arr):
    """
    redshift_arr : 各銀河の z 配列
    SMass_arr    : 各銀河の stellar mass (M_sun) 配列
    SFR_arr      : 各銀河の SFR (M_sun/yr) 配列

    戻り値:
      qflag  : quenched の bool 配列
      sbflag : starburst の bool 配列
      msflag : main_sequence の bool 配列
    """
    SMass_arr = np.asarray(SMass_arr)
    SFR_arr   = np.asarray(SFR_arr)
    redshift_arr = np.asarray(redshift_arr)

    logM = np.log10(SMass_arr)

    # sSFR としきい値
    with np.errstate(divide='ignore', invalid='ignore'):
        sSFR = SFR_arr / SMass_arr   # [1/yr]

    low, high = get_thresholds_array(redshift_arr, logM)
    #high = np.where(high>10**(-9.5), 10**(-9.5), high)
    low = 1e-11*np.ones(len(low))
    high = 10**(-9.5)*np.ones(len(high))

    # NaN 対応（しきい値が NaN のところは全部 False にしてしまう）
    valid_thr = np.isfinite(low) & np.isfinite(high) & np.isfinite(sSFR)

    qflag  = np.zeros_like(sSFR, dtype=bool)
    sbflag = np.zeros_like(sSFR, dtype=bool)
    msflag = np.zeros_like(sSFR, dtype=bool)

    qflag[valid_thr & (sSFR <  low)]  = True
    sbflag[valid_thr & (sSFR >  high)] = True
    msflag[valid_thr & (sSFR >= low) & (sSFR <= high)] = True

    return qflag, sbflag, msflag


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
def piecewise_linear(x, c, x0, m):
    """
    logsSFR <= x0: y = c  (定数)
    logsSFR >  x0: y = c + m * (x - x0)  (x0で連続な直線)
    """
    left = x <= x0
    right = ~left
    y = np.zeros_like(x)
    y[left] = x0[left]
    y[right] = c[right] + m[right] * (x[right] - x0[right])
    return y

def GetJband(redshift, SMass, SFR, quenched, starburst):
    #data = np.loadtxt('/mnt/data_cat3/yuka/repository/LIM_mock/fit_ML_J_params.txt')
    #params_all = np.zeros((len(redshift), 3))
    #params_all[redshift<0.5] = data[0]
    #params_all[(0.5<=redshift) & (redshift<1.0)] = data[1]
    #params_all[(1.0<=redshift) & (redshift<1.5)] = data[2]
    #params_all[(1.5<=redshift) & (redshift<2.0)] = data[3]
    #params_all[(2.0<=redshift) & (redshift<2.5)] = data[4]
    #params_all[(2.5<=redshift)] = data[5]


    #log_sSFR = np.log10(SFR/SMass)  # [1/yr]
    #ML = piecewise_linear(log_sSFR, params_all[:,1], params_all[:,0], params_all[:,2])
    #ML = 10**ML  # [M_sun / L_sun]
    


    data_quenched = np.loadtxt("/mnt/data_cat3/yuka/repository/LIM_mock/MLJ_Q.txt")
    z_q   = data_quenched[:,0]
    ml_q  = data_quenched[:,1]
    #sig_q = data_quenched[:,2]

    # ★ 昇順に並べ替える ★
    order_q = np.argsort(z_q)
    z_q   = z_q[order_q]
    ml_q  = ml_q[order_q]
    #sig_q = sig_q[order_q]

    quenched_med  = interp1d(z_q, ml_q,   kind='linear', bounds_error=False, fill_value='extrapolate')
    #quenched_scat = interp1d(z_q, sig_q, kind='linear', bounds_error=False, fill_value='extrapolate')
    #quenched_ml = quenched_med(redshift) + np.random.normal(0, quenched_scat(redshift))
    quenched_ml = quenched_med(redshift)

    data_SF = np.loadtxt("/mnt/data_cat3/yuka/repository/LIM_mock/MLJ_SF.txt")
    z_sf   = data_SF[:,0]
    ml_sf  = data_SF[:,1]
    #sig_sf = data_SF[:,2]

    # ★ こちらも昇順に並べ替える ★
    order_sf = np.argsort(z_sf)
    z_sf   = z_sf[order_sf]
    ml_sf  = ml_sf[order_sf]
    #sig_sf = sig_sf[order_sf]

    SF_med  = interp1d(z_sf, ml_sf,   kind='linear', bounds_error=False, fill_value='extrapolate')
   # SF_scat = interp1d(z_sf, sig_sf, kind='linear', bounds_error=False, fill_value='extrapolate')
    #SF_ml = SF_med(redshift) + np.random.normal(0, SF_scat(redshift))
    SF_ml = SF_med(redshift)

    data_SB = np.loadtxt("/mnt/data_cat3/yuka/repository/LIM_mock/MLJ_SB.txt")
    z_sb   = data_SB[:,0]
    ml_sb= data_SB[:,1]
    #sig_sf = data_SF[:,2]

    # ★ こちらも昇順に並べ替える ★
    order_sb = np.argsort(z_sb)
    z_sb   = z_sb[order_sb]
    ml_sb  = ml_sb[order_sb]
    #sig_sf = sig_sf[order_sf]

    SB_med  = interp1d(z_sb, ml_sb,   kind='linear', bounds_error=False, fill_value='extrapolate')
   # SF_scat = interp1d(z_sf, sig_sf, kind='linear', bounds_error=False, fill_value='extrapolate')
    #SF_ml = SF_med(redshift) + np.random.normal(0, SF_scat(redshift))
    SB_ml = SB_med(redshift)

    Mass_lumi = np.zeros_like(SMass)
    Mass_lumi[quenched]  = quenched_ml[quenched]
    Mass_lumi[starburst]  = SB_ml[starburst]
    MS = ~(quenched|starburst)
    Mass_lumi[MS] = SF_ml[MS]

    Jband_norm = SMass / Mass_lumi  # [Lsun] at 10pc

    J_sun_mag = 4.55
    J_flux_abs = Jband_norm * 3631.0 * Jy * 10**(-0.4*J_sun_mag) # [erg/s/cm^2/Hz] @10pc
    distance = 10.0 * pc # [cm]
    JLumi = J_flux_abs * 4.0 * np.pi * distance**2 # [erg/s/Hz]
    JLumi[JLumi<0] = 0.0
    return JLumi


def debug_ml_for_one_gal(idx, catalog, LnuJ_template):
    """
    idx 番目の銀河について：
      - Mstar, quenched/starburst フラグ
      - JLumi から J band の Lsun 換算 (Jband_norm)
      - 実効的な M/L_J
      - テンプレートの J-band Lν とのスケール
      - その JLumi から直接計算した "理論上の J AB 等級"
    を表示する。
    """
    z           = float(catalog.redshift_real[idx])
    Mstar       = float(catalog.Mstar[idx])     # [Msun]
    is_quenched = bool(catalog.quenched[idx])
    is_SB       = bool(catalog.Starburst[idx])
    JLumi       = float(catalog.JLumi[idx])     # [erg/s/Hz]

    # どのテンプレートを使っているか（あなたのマスクと合わせる）
    if is_quenched and (not is_SB):
        tpl = 'ELL'
    elif is_SB:
        tpl = 'SB'
    else:
        tpl = 'SD'

    # --- JLumi → Jband_norm [Lsun] → 実効 M/L_J ---
    J_sun_mag = 4.55        # あなたが GetJband で使っている値
    LnuJ_sun = 3631.0 * Jy * 10**(-0.4 * J_sun_mag) * 4.0 * np.pi * (10.0 * pc)**2
    #   JLumi = Jband_norm * LnuJ_sun なので
    Jband_norm = JLumi / LnuJ_sun           # [Lsun in J band at 10pc]
    ML_J       = Mstar / Jband_norm         # [Msun/Lsun]

    # --- SWIRE テンプレート側の J-band Lν と比較 ---
    LnuJ_tpl   = float(LnuJ_template[tpl])  # GetSED が返した <Lν>_J (relative→物理へのスケール用)
    scale_tpl  = JLumi / LnuJ_tpl           # テンプレ SED にかけているスケール係数

    # --- この JLumi から直接 J-band AB 等級を計算 ---
    # fν = Lν / (4π DL^2 (1+z))
    DL   = cosmo.luminosity_distance(z).to('cm').value
    fnuJ = JLumi / (4.0 * np.pi * DL**2 * (1.0 + z))   # [erg/s/cm^2/Hz]
    fnuJ = max(fnuJ, 1e-40)  # log が飛ばないように下限
    mJ_AB = -2.5 * np.log10(fnuJ) - 48.60

    print("=== debug_ml_for_one_gal ===")
    print(f"idx         = {idx}")
    print(f"z           = {z:.3f}")
    print(f"Mstar       = {Mstar:.3e} Msun")
    print(f"quenched    = {is_quenched}, starburst = {is_SB}")
    print(f"template    = {tpl}")
    print(f"JLumi       = {JLumi:.3e} erg/s/Hz")
    print(f"Jband_norm  = {Jband_norm:.3e} Lsun (J band at 10pc)")
    print(f"M/L_J_used  = {ML_J:.3f} Msun/Lsun")
    print(f"LnuJ_tpl({tpl}) = {LnuJ_tpl:.3e} (relative J-band <Lν>)")
    print(f"scale_tpl   = JLumi / LnuJ_tpl = {scale_tpl:.3e}")
    print(f"m_J(AB) from JLumi = {mJ_AB:.3f}")

###############################################################################
import numpy as np

def GetSED():
    #いくつか出力
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
    path_ELL = "/mnt/data_cat3/yuka/data/LIM_mock/SWIRE/Ell2_template_norm.sed"
    path_SD  = "/mnt/data_cat3/yuka/data/LIM_mock/SWIRE/Sdm_template_norm.sed"
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
        # AB系と整合する重み: T(ν)/ν dν
        num = np.trapz(Lnu_rel * Tnu / nu, nu)
        den = np.trapz(Tnu / nu, nu)
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

def make_fsps_sed(tage, tau=0, sfh=0, const=0.0, dust1 =1.0, dust2=0.0, logzsol=-0.1, fagn=0.0, agn_tau=10,agb_dust=1.0,add_neb_continuum=True):
    """
    Returns:
        nu: 1D array, ascending [Hz]
        Lnu: 1D array, relative Lν(ν)
    """
    c_A = c * 1e8
    sp = fsps.StellarPopulation(
        agb_dust=agb_dust,
        zcontinuous=1,
        imf_type=1,     # Chabrier
        sfh=sfh,
        const=const,
        tau=tau,
        dust_type=2,
        dust2=dust2,
        logzsol=logzsol,
        dust1=dust1,         # 若い星周りのバースト雲。あまり効かせたくなければ 0.5–1.0 程度
        add_neb_continuum=add_neb_continuum,
        gas_logu=-2.5,
        gas_logz=0.0,       # ガス金属量 ~ stellar と同じ
        add_igm_absorption=True,
        fagn = fagn,
        agn_tau = agn_tau
    )
    wave, spec = sp.get_spectrum(tage=tage, peraa=True) 
    # wave [Å], spec [L☉/Å] (per unit initial mass)
    plt.plot(wave, spec)
    plt.xlim(1e3, 1e4)
    #plt.ylim(1e-4, 10)
    plt.xscale('log')
    plt.yscale('log')
    # convert to Lnu: Lν = Lλ * λ^2 / c
    Lnu = spec * (wave**2) / c_A

    # frequency grid (descending → ascending)
    nu = c_A / wave
    order = np.argsort(nu)
    return nu[order], Lnu[order]


def GetSED_FSPS():
    """
    Return:
      nu : common frequency grid
      Lnu_rel : dict of template SEDs on the common grid
      LnuJ_template : dict of J-band averaged Lν for each template
    """
    # ==========================================================
    # 1) FSPS SED generator
    # ==========================================================
    def make(name, tage, sfh, tau=1.0, const=0.0, dust1=0.0, dust2=0.0, logzsol=-0.1, fagn=0.0, agn_tau=10.0, agb_dust=1.0, add_neb_continuum=True):
        nu_tmp, Lnu_tmp = make_fsps_sed(tage=tage, tau=tau, sfh=sfh, const=0.0, dust1 =dust1, dust2=dust2, agb_dust=agb_dust, logzsol=logzsol, add_neb_continuum=add_neb_continuum, fagn=fagn, agn_tau=agn_tau)
        return nu_tmp, Lnu_tmp

    # -------------------------
    # Quenched (ELL): 固定でOK
    # -------------------------
    #nu_Q, Lnu_Q = make("ELL", tage=2.0, sfh=0, tau=0.0, dust=0.0, zmet=5)
    nu_Q0, Lnu_Q0 = make("Q0", tage=5.0,  sfh=1, tau=0.3, logzsol=0.0, dust2=0.05, fagn=0.005, agn_tau=10)
    nu_Q1, Lnu_Q1 = make("Q1", tage=4.0,  sfh=1, tau=0.3, logzsol=0.0, dust2=0.05, fagn=0.01, agn_tau=10)
    nu_Q2, Lnu_Q2 = make("Q0", tage=2.5,  sfh=1, tau=0.3, logzsol=0.0, dust2=0.05, fagn=0.02, agn_tau=10)
    nu_Q3, Lnu_Q3 = make("Q0", tage=1.8,  sfh=1, tau=0.2, logzsol=0.0, dust2=0.05, fagn=0.03, agn_tau=10)
    nu_Q4, Lnu_Q4 = make("Q0", tage=1.0,  sfh=0, logzsol=0.0, dust2=0.03, fagn=0.05, agn_tau=10)
    nu_Q5, Lnu_Q5 = make("Q0", tage=0.6,  sfh=0, logzsol=-0.5, dust2=0.02, fagn=0.05, agn_tau=10)
    # -------------------------
    # Star-forming (Main sequence)
    # old：z≲2
    # young：z≳3 で混ぜる
    # -------------------------
    #nu_SD_old,   Lnu_SD_old   = make("SD_old", tage=4.0,  sfh=4, tau=3.0, dust1=1.0, dust2=0.3, zmet=3, add_neb_continuum=True)
    nu_MS0, Lnu_MS0 = make("MS0", tage=3.0,  sfh=4, tau=2.0, logzsol=-0.1, dust2=0.2, fagn=0.02, agn_tau=15, agb_dust=1.0)
    nu_MS1, Lnu_MS1 = make("MS1", tage=2.0,  sfh=4, tau=2.0, logzsol=-0.3, dust2=0.25, fagn=0.05, agn_tau=15, agb_dust=1.4)
    nu_MS2, Lnu_MS2 = make("MS2", tage=1.2,  sfh=4, tau=1.8, logzsol=-0.5, dust2=0.3, fagn=0.07, agn_tau=20, agb_dust=1.6)
    nu_MS3, Lnu_MS3 = make("MS3", tage=1.0,  sfh=4, tau=1.5, logzsol=-0.6, dust2=0.35, fagn=0.10, agn_tau=15, agb_dust=1.8)
    nu_MS4, Lnu_MS4 = make("MS4", tage=0.7,  sfh=4, tau=1.0, logzsol=-0.8, dust2=0.35, fagn=0.15, agn_tau=25, agb_dust=2.0)
    nu_MS5, Lnu_MS5 = make("MS5", tage=0.4,  sfh=4, tau=0.8, logzsol=-1.0, dust2=0.3, fagn=0.20, agn_tau=30, agb_dust=2.2)
    # -------------------------
    # Starburst
    # -------------------------
    #nu_SB_old,   Lnu_SB_old   = make("SB_old",   tage=0.3,  sfh=1, tau=1.0, dust1=1.0)
    nu_SB0, Lnu_SB0 = make("SB0", tage=0.15, sfh=1, tau=10.0, logzsol=-0.3, dust2=0.4, dust1=1.0, fagn=0.1, agn_tau=20, agb_dust=1.5)
    nu_SB1, Lnu_SB1 = make("SB1", tage=0.12, sfh=1, tau=10.0, logzsol=-0.5, dust2=0.5, dust1=1.2, fagn=0.12, agn_tau=25, agb_dust=1.7)
    nu_SB2, Lnu_SB2 = make("SB2", tage=0.10, sfh=1, tau=12.0, logzsol=-0.7, dust2=0.6, dust1=1.3, fagn=0.15, agn_tau=25, agb_dust=1.9)
    nu_SB3, Lnu_SB3 = make("SB3", tage=0.08, sfh=1, tau=15.0, logzsol=-0.8, dust2=0.7, dust1=1.5, fagn=0.19, agn_tau=27, agb_dust=2.0)
    nu_SB4, Lnu_SB4 = make("SB4", tage=0.06, sfh=1, tau=20.0, logzsol=-1.0, dust2=0.8, dust1=1.7, fagn=0.22, agn_tau=35, agb_dust=2.2)
    nu_SB5, Lnu_SB5 = make("SB5", tage=0.05, sfh=2, tau=5.0, logzsol=-1.3, dust2=1.0, dust1=2.0, fagn=0.25, agn_tau=50, agb_dust=2.3)

    # ==========================================================
    # 共通 frequency grid（log spacing）
    # ==========================================================
    nu_min = min(nu_Q0.min(), nu_MS0.min(), nu_SB0.min())
    nu_max = max(nu_Q0.max(), nu_MS0.max(), nu_SB0.max())
    nu = np.logspace(np.log10(nu_min), np.log10(nu_max), 5000)

    # ==========================================================
    # すべてのテンプレートを共通 ν grid に補間
    # ==========================================================
    Lnu_rel = {
        'Q0':      np.interp(nu, nu_Q0,        Lnu_Q0,        left=0, right=0),
        'Q1':      np.interp(nu, nu_Q1,        Lnu_Q1,        left=0, right=0),
        'Q2':      np.interp(nu, nu_Q2,        Lnu_Q2,        left=0, right=0),
        'Q3':      np.interp(nu, nu_Q3,        Lnu_Q3,        left=0, right=0),
        'Q4':      np.interp(nu, nu_Q4,        Lnu_Q4,        left=0, right=0),
        'Q5':      np.interp(nu, nu_Q5,        Lnu_Q5,        left=0, right=0),
        #'SD_old':   np.interp(nu, nu_SD_old,   Lnu_SD_old,   left=0, right=0),
        'MS0': np.interp(nu, nu_MS0, Lnu_MS0, left=0, right=0),
        'MS1': np.interp(nu, nu_MS1, Lnu_MS1, left=0, right=0),
        'MS2': np.interp(nu, nu_MS2, Lnu_MS2, left=0, right=0),
        'MS3': np.interp(nu, nu_MS3, Lnu_MS3, left=0, right=0),
        'MS4': np.interp(nu, nu_MS4, Lnu_MS4, left=0, right=0),
        'MS5': np.interp(nu, nu_MS5, Lnu_MS5, left=0, right=0),
        #'SB_old':   np.interp(nu, nu_SB_old,   Lnu_SB_old,   left=0, right=0),
        'SB0': np.interp(nu, nu_SB0, Lnu_SB0, left=0, right=0),
        'SB1': np.interp(nu, nu_SB1, Lnu_SB1, left=0, right=0),
        'SB2': np.interp(nu, nu_SB2, Lnu_SB2, left=0, right=0),
        'SB3': np.interp(nu, nu_SB3, Lnu_SB3, left=0, right=0),
        'SB4': np.interp(nu, nu_SB4, Lnu_SB4, left=0, right=0),
        'SB5': np.interp(nu, nu_SB5, Lnu_SB5, left=0, right=0),
    }

    # ==========================================================
    # J-band throughput 読み込み → νへ変換
    # ==========================================================
    c_A_per_s = 2.99792458e18  # Å/s
    J_path = '/mnt/data_cat3/yuka/repository/LIM_mock/Jband_throughput.txt'
    wJ, tJ = np.loadtxt(J_path, unpack=True)

    # wJ が μm の場合（あなたのコードにならう）
    wave_J = wJ * 1e4  # Å
    T_J = tJ.copy()

    nu_J = c_A_per_s / wave_J
    order = np.argsort(nu_J)
    nu_J = nu_J[order]
    T_J  = T_J[order]

    Tnu_J = np.interp(nu, nu_J, T_J, left=0.0, right=0.0)

    # ==========================================================
    # J-band averaged Lν を計算
    # ==========================================================
    def band_average_Lnu(Lnu_rel):
        num = np.trapz(Lnu_rel * Tnu_J / nu, nu)
        den = np.trapz(Tnu_J / nu, nu)
        return num / den if den > 0 else 0.0

    LnuJ_template = {
        key: band_average_Lnu(val) for key, val in Lnu_rel.items()
    }

    return nu, Lnu_rel, LnuJ_template


##############################################################################
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

def bandavg_Lnu_any_filter(nu_tpl, Lnu_tpl, nu_obs, Tnu_obs, z_gal):
    """
    テンプレ Lν(ν) を観測フィルタ (ν_obs, Tnu_obs) で、各銀河の赤方偏移 z_gal に対して
    帯域平均 <Lν> を返す（形状: [N_gal] または [N_gal, N_filt]）。
    """
    # 累積積分 I(ν)=∫Lν dν
    I = cumulative_trapezoid(Lnu_tpl, nu_tpl, initial=0.0)
    # 透過の積分 J(ν)=∫T dν（観測フィルタ側も積分）
    J = cumulative_trapezoid(Tnu_obs, nu_obs, initial=0.0)

    # 観測フィルタの有効周波数範囲（端）
    nu1_obs = nu_obs.min()
    nu2_obs = nu_obs.max()

    # 各銀河のレスト周波数範囲（(1+z)*ν）
    nu1_rest = (1.0 + z_gal) * nu1_obs
    nu2_rest = (1.0 + z_gal) * nu2_obs

    # I(ν) をレスト側端点に補間
    I1 = np.interp(nu1_rest, nu_tpl, I, left=I[0], right=I[-1])
    I2 = np.interp(nu2_rest, nu_tpl, I, left=I[0], right=I[-1])

    # 観測側の T の積分（分母）— これは z に依らず一定
    J_tot = J[-1] - J[0]
    # 単純に <Lν> ≈ ∫ Lν_rest(ν) T(ν/(1+z)) dν / ∫ T dν を
    # 近似的に I2-I1 と J_tot で表現（テンプレは十分細かい前提）。
    # より正確には (ν→ν'=(1+z)ν の写像で T をシフト) だが、実装簡素化版。
    Lnu_avg = (I2 - I1) / (nu2_rest - nu1_rest)  # 平均化
    # 透過の重みを厳密に入れたい場合は、各銀河ごとに T(ν_obs) を ν_rest=(1+z)ν_obs に
    # 引き伸ばして Lν×T を数値積分する実装に差し替えてください。

    return Lnu_avg

import numpy as np

def bandavg_Lnu_any_filter_exact(
    nu_tpl, Lnu_tpl,     # rest-frame Lnu(ν), ν ascending
    nu_obs, Tnu_obs,     # observed filter curve in ν_obs (ascending)
    z_array,
    batch_size=2000
):

    # AB 系なので重みは T(νobs)/νobs dνobs
    denom = np.trapz(Tnu_obs / nu_obs, nu_obs)

    out = np.empty_like(z_array, dtype=float)
    N = len(z_array)

    for s in range(0, N, batch_size):
        e = min(N, s + batch_size)
        z = z_array[s:e]                   # (B,)
        fac = 1.0/(1.0+z)                  # ν_obs = ν_rest/(1+z)

        # 観測系ν_obsにおけるフィルタを rest-frame ν に引き戻す
        nu_rest = nu_tpl[None,:]          # (1,Ntpl)
        nu_obs_here = nu_rest * fac[:,None]  # (B,Ntpl)

        # フィルタ値 T(ν_obs) を ν_rest grid 上に補間
        T_here = np.interp(nu_obs_here, nu_obs, Tnu_obs, left=0, right=0)

        # numerator = ∫ Lν(νrest) * T(νobs) / νrest dνrest
        num = np.trapz(Lnu_tpl[None,:] * T_here / nu_rest, nu_rest, axis=1)

        # AB 系では最終的に (1+z) で割る（dνobs = dνrest/(1+z)）
        out[s:e] = num / ((1.0+z) * denom)

    return out



def plot_Nz_for_magnitude_bins(
    redshift, 
    mag, 
    mag_bins=[(18,19),(19,20),(20,21),(21,22),(22,23),(23,24)],
    label_prefix="JWST_F115W"
):

    plt.figure(figsize=(8,6))
    bins = np.linspace(0.0, 3.0, 31)

    for (m1, m2) in mag_bins:
        sel = (mag >= m1) & (mag < m2)
        if sel.sum() == 0:
            continue

        z_sel = redshift[sel]

        plt.hist(z_sel, bins=bins, histtype='step',
                 label=f"{label_prefix}: {m1}–{m2} (N={sel.sum()})")

    plt.xlabel("Redshift")
    plt.ylabel("Normalized counts")
    plt.legend()
    plt.xlim(0.0, 3.0)
    plt.title(f"N(z) per magnitude bin — {label_prefix}")
    plt.grid()
    plt.show()

##############################################################################
def mix_factor(z, z1=1.5, z2=4.0):
    f = (z - z1) / (z2 - z1)
    return np.clip(f, 0, 1)

def Compare_mag(filter_names):
    catalog = Catalog()
    input_fname = '/mnt/data_cat3/yuka/data/Cosmos_Web/'
    catalog.loadcatalog(input_fname, name="JWST")
    lumi_dist = cosmo.luminosity_distance(np.array(catalog.redshift_real)).to(u.cm)
    lumi_dist = lumi_dist.value

    ############################Stellar continuum############################
    zobs = np.asarray(catalog.redshift_obs)      # shape (Ngal,)

    mask_SB  = catalog.Starburst.astype(bool)
    mask_Q   = catalog.quenched.astype(bool) & (~mask_SB)
    mask_Sd  = (~mask_Q) & (~mask_SB)
    #nu_tpl, Lnu_rel, LnuJ_template = GetSED() # [Hz]
    nu_tpl, Lnu_rel, LnuJ_template = GetSED_FSPS()
    idx_test = np.argmax(catalog.JLumi)
    #debug_ml_for_one_gal(idx_test, catalog, LnuJ_template)
    #print(f"ELL:{LnuJ_template['ELL']} SD:{LnuJ_template['SD']}, SB:{LnuJ_template['SB']}")

    # スケール係数（J-band 平均 Lν の比）
    log_scale = np.zeros_like(catalog.JLumi, dtype=float)
    f = mix_factor(catalog.redshift_real)

    redshift_bin0 = catalog.redshift_real < 0.5
    redshift_bin1 = (catalog.redshift_real >= 0.5) & (catalog.redshift_real < 1.0)
    redshift_bin2 = (catalog.redshift_real >= 1.0) & (catalog.redshift_real < 1.5)
    redshift_bin3 = (catalog.redshift_real >= 1.5) & (catalog.redshift_real < 2.0)
    redshift_bin4 = (catalog.redshift_real >= 2.0) & (catalog.redshift_real < 3.0)
    redshift_bin5 = (catalog.redshift_real >= 3.0)
    if np.any(mask_Q):
        #log_scale[mask_Q]  = np.log10(catalog.JLumi[mask_Q])  - np.log10(LnuJ_template['ELL'])
        log_scale[mask_Q&redshift_bin0]  = np.log10(catalog.JLumi[mask_Q&redshift_bin0])  - np.log10(LnuJ_template['Q0'])
        log_scale[mask_Q&redshift_bin1]  = np.log10(catalog.JLumi[mask_Q&redshift_bin1])  - np.log10(LnuJ_template['Q1'])
        log_scale[mask_Q&redshift_bin2]  = np.log10(catalog.JLumi[mask_Q&redshift_bin2])  - np.log10(LnuJ_template['Q2'])
        log_scale[mask_Q&redshift_bin3]  = np.log10(catalog.JLumi[mask_Q&redshift_bin3])  - np.log10(LnuJ_template['Q3'])
        log_scale[mask_Q&redshift_bin4]  = np.log10(catalog.JLumi[mask_Q&redshift_bin4])  - np.log10(LnuJ_template['Q4'])
        log_scale[mask_Q&redshift_bin5]  = np.log10(catalog.JLumi[mask_Q&redshift_bin5])  - np.log10(LnuJ_template['Q5'])
        #plt.hist(log_scale[mask_Q], bins=30, alpha=0.5, label='ELL')
        #plt.yscale('log')
    if np.any(mask_Sd):
        log_scale[mask_Sd&redshift_bin0]  = np.log10(catalog.JLumi[mask_Sd&redshift_bin0])  - np.log10(LnuJ_template['MS0'])
        log_scale[mask_Sd&redshift_bin1]  = np.log10(catalog.JLumi[mask_Sd&redshift_bin1])  - np.log10(LnuJ_template['MS1'])
        log_scale[mask_Sd&redshift_bin2]  = np.log10(catalog.JLumi[mask_Sd&redshift_bin2])  - np.log10(LnuJ_template['MS2'])
        log_scale[mask_Sd&redshift_bin3]  = np.log10(catalog.JLumi[mask_Sd&redshift_bin3])  - np.log10(LnuJ_template['MS3'])
        log_scale[mask_Sd&redshift_bin4]  = np.log10(catalog.JLumi[mask_Sd&redshift_bin4])  - np.log10(LnuJ_template['MS4'])
        log_scale[mask_Sd&redshift_bin5]  = np.log10(catalog.JLumi[mask_Sd&redshift_bin5])  - np.log10(LnuJ_template['MS5'])
        #plt.hist(log_scale[mask_Sd], bins=30, alpha=0.5, label='SD')
        #plt.yscale('log')
    if np.any(mask_SB):
        log_scale[mask_SB&redshift_bin0]  = np.log10(catalog.JLumi[mask_SB&redshift_bin0])  - np.log10(LnuJ_template['SB0'])
        log_scale[mask_SB&redshift_bin1]  = np.log10(catalog.JLumi[mask_SB&redshift_bin1])  - np.log10(LnuJ_template['SB1'])
        log_scale[mask_SB&redshift_bin2]  = np.log10(catalog.JLumi[mask_SB&redshift_bin2])  - np.log10(LnuJ_template['SB2'])
        log_scale[mask_SB&redshift_bin3]  = np.log10(catalog.JLumi[mask_SB&redshift_bin3])  - np.log10(LnuJ_template['SB3'])
        log_scale[mask_SB&redshift_bin4]  = np.log10(catalog.JLumi[mask_SB&redshift_bin4])  - np.log10(LnuJ_template['SB4'])
        log_scale[mask_SB&redshift_bin5]  = np.log10(catalog.JLumi[mask_SB&redshift_bin5])  - np.log10(LnuJ_template['SB5'])
        #plt.hist(log_scale[mask_SB], bins=30, alpha=0.5, label='SB')
        #plt.yscale('log')
    #plt.legend()

    ########################### Ha intensity
    SFR = catalog.SFR

    log_lumi_Ha = np.log10(SFR) + np.log10( 8.3e40 )
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_Ha))

    #lumi_Ha = 10**(log_lumi + log_lumi_scatter)
    log_flux_Ha = (log_lumi_Ha + log_lumi_scatter) - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_Ha = 10**log_flux_Ha

    ########################### Pa intensity
    log_lumi_Pa = np.log10(SFR) + np.log10( 8.3e40 ) - 0.6
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_Pa))

    log_flux_Pa = (log_lumi_Pa + log_lumi_scatter)  - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_Pa = 10**log_flux_Pa

    ########################### [OIII] intensity
    log_lumi_OIII = np.log10(SFR) + np.log10( 1.32e41 )
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_OIII))

    log_flux_OIII = (log_lumi_OIII + log_lumi_scatter)  - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_OIII = 10**log_flux_OIII

    zmin, zmax = zobs.min(), zobs.max()
    print(f'Redshift range: {zmin} - {zmax}')
    # z グリッドの点数は好みで (100–300 あれば十分なことが多い)
    z_grid = np.linspace(zmin, zmax, 200)

    # --- フィルタファイルのパス ---
    filter_files = {}
    for fname in filter_names:
        if 'Spitzer' in fname:
            filter_files[fname] = f"/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125{fname[-3:].lower()}trans_full.txt"
        else:
            filter_files[fname] = f"/mnt/data_cat3/yuka/data/nircam_throughputs/mean_throughputs/{fname}_May2024_mean_system_throughput.txt"

    filters = {}
    for fname, path in filter_files.items():
        nu_flt, Tnu_flt = load_filter_throughput_txt(fname, path)
        denom_T = np.trapz(Tnu_flt, nu_flt)
        filters[fname] = {
            "nu": nu_flt,
            "T":  Tnu_flt,
            "denom_T": max(denom_T, 1e-300),
        }

    # precomputed_bandavg[filter][template] = {"Lnu": array(z_grid), "nuLnu": array(z_grid)}
    precomputed_bandavg = {}
    for fname, fdat in filters.items():
        nu_flt = fdat["nu"]
        Tnu_flt = fdat["T"]

        precomputed_bandavg[fname] = {}
        for tpl_name, Lnu_tpl in Lnu_rel.items():   # Lnu_rel = {'ELL':..., 'SD':..., 'SB':...}
            # <Lnu>_band(z) を z_grid 上でだけ正確に計算
            Lnu_band_z = bandavg_Lnu_any_filter_exact(
                nu_tpl, Lnu_tpl, nu_flt, Tnu_flt, z_grid
            )
            nuLnu_band_z = bandavg_Lnu_any_filter_exact(
                nu_tpl, nu_tpl * Lnu_tpl, nu_flt, Tnu_flt, z_grid
            )
            precomputed_bandavg[fname][tpl_name] = {
                "Lnu":   Lnu_band_z,
                "nuLnu": nuLnu_band_z,
            }

    def Get_mag_all(filter_name):
        sSFR = np.log10(catalog.SFR/catalog.Mstar)
        fdat = filters[filter_name]
        nu_flt   = fdat["nu"]
        Tnu_flt  = fdat["T"]
        denom_T  = fdat["denom_T"]

        logLnu_band    = np.zeros(len(zobs), dtype=float)
        log_nuLnu_band = np.zeros(len(zobs), dtype=float)

        # ---------- 連続光：タイプ×z-binごとにテンプレートを選ぶ ----------

        # 便利なタプル（binごとに suffix を付けたキーを作る）
        zbin_masks = [redshift_bin0, redshift_bin1, redshift_bin2,
                      redshift_bin3, redshift_bin4, redshift_bin5]

        # Quenched: Q0..Q5
        for i_bin, zmask in enumerate(zbin_masks):
            sel = np.where(mask_Q & zmask)[0]
            if sel.size == 0:
                continue
            key = f"Q{i_bin}"
            Lnu_z   = np.interp(zobs[sel], z_grid,
                                precomputed_bandavg[filter_name][key]["Lnu"])
            nuLnu_z = np.interp(zobs[sel], z_grid,
                                precomputed_bandavg[filter_name][key]["nuLnu"])
            logLnu_band[sel]    = np.log10(Lnu_z)   + log_scale[sel]
            log_nuLnu_band[sel] = np.log10(nuLnu_z) + log_scale[sel]

        # Main sequence: MS0..MS5
        for i_bin, zmask in enumerate(zbin_masks):
            sel = np.where(mask_Sd & zmask)[0]
            if sel.size == 0:
                continue
            key = f"MS{i_bin}"
            Lnu_z   = np.interp(zobs[sel], z_grid,
                                precomputed_bandavg[filter_name][key]["Lnu"])
            nuLnu_z = np.interp(zobs[sel], z_grid,
                                precomputed_bandavg[filter_name][key]["nuLnu"])
            logLnu_band[sel]    = np.log10(Lnu_z)   + log_scale[sel]
            log_nuLnu_band[sel] = np.log10(nuLnu_z) + log_scale[sel]

        # Starburst: SB0..SB5
        for i_bin, zmask in enumerate(zbin_masks):
            sel = np.where(mask_SB & zmask)[0]
            if sel.size == 0:
                continue
            key = f"SB{i_bin}"
            Lnu_z   = np.interp(zobs[sel], z_grid,
                                precomputed_bandavg[filter_name][key]["Lnu"])
            nuLnu_z = np.interp(zobs[sel], z_grid,
                                precomputed_bandavg[filter_name][key]["nuLnu"])
            logLnu_band[sel]    = np.log10(Lnu_z)   + log_scale[sel]
            log_nuLnu_band[sel] = np.log10(nuLnu_z) + log_scale[sel]

        # --- 連続光の fν, νfν ----
        log_fnu_cont = (
            logLnu_band
            - np.log10(4.0*np.pi)
            - 2.0*np.log10(lumi_dist)
            - np.log10(1.0 + zobs)
        )
        fnu_cont = 10**log_fnu_cont

        log_nufnu_cont = (
            log_nuLnu_band
            - np.log10(4.0*np.pi)
            - 2.0*np.log10(lumi_dist)
        )
        nufnu_cont = 10**log_nufnu_cont

        # ===== 輝線部分（ここは元のままで OK） =====
        lam_Ha_um   = c/ha/1e5
        lam_Pa_um   = c/pa/1e5
        lam_OIII_um = c/Oiii/1e5

        c_um_s = 2.99792458e14  # [μm/s]
        nu_Ha_rest   = c_um_s / lam_Ha_um
        nu_Pa_rest   = c_um_s / lam_Pa_um
        nu_OIII_rest = c_um_s / lam_OIII_um

        nu_Ha_obs   = nu_Ha_rest   / (1.0 + zobs)
        nu_Pa_obs   = nu_Pa_rest   / (1.0 + zobs)
        nu_OIII_obs = nu_OIII_rest / (1.0 + zobs)

        T_Ha   = np.interp(nu_Ha_obs,   nu_flt, Tnu_flt, left=0.0, right=0.0)
        T_Pa   = np.interp(nu_Pa_obs,   nu_flt, Tnu_flt, left=0.0, right=0.0)
        T_OIII = np.interp(nu_OIII_obs, nu_flt, Tnu_flt, left=0.0, right=0.0)

        F_Ha   = 10**log_flux_Ha
        F_Pa   = 10**log_flux_Pa
        F_OIII = 10**log_flux_OIII

        fnu_line = (F_Ha*T_Ha + F_Pa*T_Pa + F_OIII*T_OIII) / denom_T
        fnu_total = fnu_cont + fnu_line

        nufnu_line = (nu_Ha_obs * F_Ha * T_Ha +
                    nu_Pa_obs * F_Pa * T_Pa +
                    nu_OIII_obs * F_OIII * T_OIII) / denom_T
        nufnu_total = nufnu_cont + nufnu_line

        # --- AB mag & LF ---
        m_ab = -2.5*np.log10(np.clip(fnu_total, 1e-300, None)) - 48.60

        plot_Nz_for_magnitude_bins(catalog.redshift_real, m_ab, label_prefix=f"{filter_name} distribution")
        if filter_name == 'F115W':
            quenched = catalog.quenched
            SB = catalog.Starburst
            MS = ~(quenched|SB)
            plot_Nz_for_magnitude_bins(catalog.redshift_real, catalog.F115W, label_prefix=f"{filter_name} distribution true")
            plt.figure()
            plt.title(filter_name)
            plt.scatter(
                m_ab,                      # x 軸
                catalog.F115W,             # y 軸
                s=1,
                c=catalog.redshift_real,   # 色：赤方偏移
                cmap='viridis',            # 好きなカラーマップに変更OK
                alpha=0.7
                )
            plt.colorbar(label='redshift (z)')
            x = np.linspace(15, 40, 50)
            y = x
            plt.plot(x, y, color='red', linestyle='dashed')
            plt.xlim(15, 40)
            plt.ylim(15, 40)
            plt.xlabel("Computed mag")
            plt.ylabel("Catalog mag")

            plt.figure()
            plt.title(filter_name)
            plt.scatter(m_ab,catalog.F115W, s=1, color='green', alpha=0.3)
            #plt.scatter(m_ab[SB],catalog.F115W[SB], s=1, color='blue', alpha=0.3)
            #plt.scatter(m_ab[quenched],catalog.F115W[quenched], s=1, color='red', alpha=0.3)
            plt.colorbar(label='redshift (z)')
            x = np.linspace(15, 40, 50)
            y = x
            plt.plot(x, y, color='red', linestyle='dashed')
            plt.xlim(15, 40)
            plt.ylim(15, 40)
            plt.xlabel("Computed mag")
            plt.ylabel("Catalog mag")

            plt.figure()
            bins_m = np.arange(15, 31, 0.5)
            hist_mock, _ = np.histogram(m_ab, bins=bins_m)
            hist_JWST, _ = np.histogram(catalog.F115W, bins=bins_m)
            m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
            plt.scatter(m_cent, hist_mock, color='blue')
            plt.scatter(m_cent, hist_JWST, color='red')
            plt.yscale('log')

        elif filter_name == 'F150W':
            quenched = catalog.quenched
            SB = catalog.Starburst
            MS = ~(quenched|SB)
            plot_Nz_for_magnitude_bins(catalog.redshift_real, catalog.F150W, label_prefix=f"{filter_name} distribution true")
            plt.figure()
            plt.title(filter_name)
            plt.scatter(
                m_ab,                      # x 軸
                catalog.F150W,             # y 軸
                s=1,
                c=catalog.redshift_real,   # 色：赤方偏移
                cmap='viridis',            # 好きなカラーマップに変更OK
                alpha=0.7
                )
            plt.colorbar(label='redshift (z)')
            x = np.linspace(15, 40, 50)
            y = x
            plt.plot(x, y, color='red', linestyle='dashed')
            plt.xlim(15, 40)
            plt.ylim(15, 40)
            plt.xlabel("Computed mag")
            plt.ylabel("Catalog mag")

            plt.figure()
            plt.title(filter_name)
            plt.scatter(m_ab,catalog.F150W, s=1, color='green', alpha=0.3)
            #plt.scatter(m_ab[SB],catalog.F150W[SB], s=1, color='blue', alpha=0.3)
            #plt.scatter(m_ab[quenched],catalog.F150W[quenched], s=1, color='red', alpha=0.3)
            plt.colorbar(label='redshift (z)')
            x = np.linspace(15, 40, 50)
            y = x
            plt.plot(x, y, color='red', linestyle='dashed')
            plt.xlim(15, 40)
            plt.ylim(15, 40)
            plt.xlabel("Computed mag")
            plt.ylabel("Catalog mag")

            plt.figure()
            bins_m = np.arange(15, 31, 0.5)
            hist_mock, _ = np.histogram(m_ab, bins=bins_m)
            hist_JWST, _ = np.histogram(catalog.F150W, bins=bins_m)
            m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
            plt.scatter(m_cent, hist_mock, color='blue')
            plt.scatter(m_cent, hist_JWST, color='red')
            plt.yscale('log')

        elif filter_name == 'F277W':
            quenched = catalog.quenched
            SB = catalog.Starburst
            MS = ~(quenched|SB)
            plot_Nz_for_magnitude_bins(catalog.redshift_real, catalog.F277W, label_prefix=f"{filter_name} distribution true")
            plt.figure()
            plt.title(filter_name)
            plt.scatter(
                m_ab,                      # x 軸
                catalog.F277W,             # y 軸
                s=1,
                c=catalog.redshift_real,   # 色：赤方偏移
                cmap='viridis',            # 好きなカラーマップに変更OK
                alpha=0.7
                )
            plt.colorbar(label='redshift (z)')
            x = np.linspace(15, 40, 50)
            y = x
            plt.plot(x, y, color='red', linestyle='dashed')
            plt.xlim(15, 40)
            plt.ylim(15, 40)
            plt.xlabel("Computed mag")
            plt.ylabel("Catalog mag")

            plt.figure()
            plt.title(filter_name)
            plt.scatter(m_ab,catalog.F277W, s=1, color='green', alpha=0.3)
            #plt.scatter(m_ab[SB],catalog.F277W[SB], s=1, color='blue', alpha=0.3)
            #plt.scatter(m_ab[quenched],catalog.F277W[quenched], s=1, color='red', alpha=0.3)
            plt.colorbar(label='redshift (z)')
            x = np.linspace(15, 40, 50)
            y = x
            plt.plot(x, y, color='red', linestyle='dashed')
            plt.xlim(15, 40)
            plt.ylim(15, 40)
            plt.xlabel("Computed mag")
            plt.ylabel("Catalog mag")

            plt.figure()
            bins_m = np.arange(15, 31, 0.5)
            hist_mock, _ = np.histogram(m_ab, bins=bins_m)
            hist_JWST, _ = np.histogram(catalog.F277W, bins=bins_m)
            m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
            plt.scatter(m_cent, hist_mock, color='blue')
            plt.scatter(m_cent, hist_JWST, color='red')
            plt.yscale('log')

        elif filter_name == 'F444W':
            quenched = catalog.quenched
            SB = catalog.Starburst
            MS = ~(quenched|SB)
            plot_Nz_for_magnitude_bins(catalog.redshift_real, catalog.F444W, label_prefix=f"{filter_name} distribution true")
            plt.figure()
            plt.title(filter_name)
            plt.scatter(
                m_ab,                      # x 軸
                catalog.F444W,             # y 軸
                s=1,
                c=catalog.redshift_real,   # 色：赤方偏移
                cmap='viridis',            # 好きなカラーマップに変更OK
                alpha=0.7
                )
            plt.colorbar(label='redshift (z)')
            x = np.linspace(15, 40, 50)
            y = x
            plt.plot(x, y, color='red', linestyle='dashed')
            plt.xlim(15, 40)
            plt.ylim(15, 40)
            plt.xlabel("Computed mag")
            plt.ylabel("Catalog mag")

            plt.figure()
            plt.title(filter_name)
            plt.scatter(m_ab, catalog.F444W, s=1, color='green', alpha=0.3)
            #plt.scatter(m_ab[SB],catalog.F444W[SB], s=1, color='blue', alpha=0.3)
            #plt.scatter(m_ab[quenched],catalog.F444W[quenched], s=1, color='red', alpha=0.3)
            plt.colorbar(label='redshift (z)')
            x = np.linspace(15, 40, 50)
            y = x
            plt.plot(x, y, color='red', linestyle='dashed')
            plt.xlim(15, 40)
            plt.ylim(15, 40)
            plt.xlabel("Computed mag")
            plt.ylabel("Catalog mag")

            plt.figure()
            bins_m = np.arange(15, 31, 0.5)
            hist_mock, _ = np.histogram(m_ab, bins=bins_m)
            hist_JWST, _ = np.histogram(catalog.F444W, bins=bins_m)
            m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
            plt.scatter(m_cent, hist_mock, color='blue')
            plt.scatter(m_cent, hist_JWST, color='red')
            plt.yscale('log')

        return nufnu_total
    
    results = [Get_mag_all(f) for f in filter_names]
    return results
    
    

def make_mock_spitzer(names, filter_names):
    Ch1_resolution = 1.66 #arcsec
    Ch2_resolution = 1.72 #arcsec

    Ch1_path = "/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125ch1trans_full.txt"
    Ch2_path = "/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125ch2trans_full.txt"

    Ch1_data = np.loadtxt(Ch1_path)
    Ch2_data = np.loadtxt(Ch2_path)

    Ch1_wavelength = Ch1_data[:,0] #um
    Ch1_wavelength = Ch1_wavelength[::-1]
    Ch1_nu = c / (Ch1_wavelength * 1e-4)  # Hz

    Ch1_throughput = Ch1_data[:,1]
    Ch1_throughput = Ch1_throughput[::-1]

    Ch2_wavelength = Ch2_data[:,0] #um
    Ch2_wavelength = Ch2_wavelength[::-1]
    Ch2_nu = c / (Ch2_wavelength * 1e-4)  # Hz

    Ch2_throughput = Ch2_data[:,1]
    Ch2_throughput = Ch2_throughput[::-1]

    catalog = Catalog()
    input_fname = '/mnt/data_cat3/yuka/data/lightcone_5400sec.txt'
    catalog.loadcatalog(input_fname, name="TNG")
    #input_fname = '/mnt/data_cat3/yuka/data/SIDES/'
    #catalog.loadcatalogs(input_fname, names=names)

    lumi_dist = cosmo.luminosity_distance(np.array(catalog.redshift_real)).to(u.cm)
    lumi_dist = lumi_dist.value

    ra = np.array(catalog.x) * 3600.0 #[arcsec]
    ra -= np.min(ra)
    dec = np.array(catalog.y) * 3600.0 #[arcsec]
    dec -= np.min(dec)

    Ch1_ix = np.floor(ra / Ch1_resolution).astype(np.int64)
    Ch1_iy = np.floor(dec / Ch1_resolution).astype(np.int64)

    Ch2_ix = np.floor(ra / Ch2_resolution).astype(np.int64)
    Ch2_iy = np.floor(dec / Ch2_resolution).astype(np.int64)
    #print(ix)
    
    ##initialize
    #Nx = int(2.0*3600*params.radius / params.resolution)
    ra_size = np.max(ra) - np.min(ra)
    dec_size = np.max(dec) - np.min(dec)
    Ch1_Nx = int(np.ceil( ra_size / Ch1_resolution ))
    Ch2_Nx = int(np.ceil( ra_size / Ch2_resolution ))
    #Ny = int(np.ceil( (np.max(dec) - np.min(dec)) / Ch1_resolution ))
    
    Ch1_npix = np.array([Ch1_Nx, Ch1_Nx])
    Ch2_npix = np.array([Ch2_Nx, Ch2_Nx])

    ############################Stellar continuum############################
    zobs = np.asarray(catalog.redshift_obs)      # shape (Ngal,)

    mask_SB  = catalog.Starburst.astype(bool)
    mask_Q   = catalog.quenched.astype(bool) & (~mask_SB)
    mask_Sd  = (~mask_Q) & (~mask_SB)
    nu_tpl, Lnu_rel, LnuJ_template = GetSED_FSPS() # [Hz]
    idx_test = np.argmax(catalog.JLumi)
    debug_ml_for_one_gal(idx_test, catalog, LnuJ_template)
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

    ########################### Ha intensity
    SFR = catalog.SFR

    log_lumi_Ha = np.log10(SFR) + np.log10( 8.3e40 )
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_Ha))

    #lumi_Ha = 10**(log_lumi + log_lumi_scatter)
    log_flux_Ha = (log_lumi_Ha + log_lumi_scatter) - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_Ha = 10**log_flux_Ha

    ########################### Pa intensity
    log_lumi_Pa = np.log10(SFR) + np.log10( 8.3e40 ) - 0.6
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_Pa))

    log_flux_Pa = (log_lumi_Pa + log_lumi_scatter)  - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_Pa = 10**log_flux_Pa

    ########################### [OIII] intensity
    log_lumi_OIII = np.log10(SFR) + np.log10( 1.32e41 )
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_OIII))

    log_flux_OIII = (log_lumi_OIII + log_lumi_scatter)  - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_OIII = 10**log_flux_OIII

    zmin, zmax = zobs.min(), zobs.max()
    print(f'Redshift range: {zmin} - {zmax}')
    # z グリッドの点数は好みで (100–300 あれば十分なことが多い)
    z_grid = np.linspace(zmin, zmax, 200)

    # --- フィルタファイルのパス ---
    filter_files = {}
    for fname in filter_names:
        if 'Spitzer' in fname:
            filter_files[fname] = f"/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125{fname[-3:].lower()}trans_full.txt"
        else:
            filter_files[fname] = f"/mnt/data_cat3/yuka/data/nircam_throughputs/mean_throughputs/{fname}_May2024_mean_system_throughput.txt"

    filters = {}
    for fname, path in filter_files.items():
        nu_flt, Tnu_flt = load_filter_throughput_txt(fname, path)
        denom_T = np.trapz(Tnu_flt, nu_flt)
        filters[fname] = {
            "nu": nu_flt,
            "T":  Tnu_flt,
            "denom_T": max(denom_T, 1e-300),
        }

    # precomputed_bandavg[filter][template] = {"Lnu": array(z_grid), "nuLnu": array(z_grid)}
    precomputed_bandavg = {}
    for fname, fdat in filters.items():
        nu_flt = fdat["nu"]
        Tnu_flt = fdat["T"]

        precomputed_bandavg[fname] = {}
        for tpl_name, Lnu_tpl in Lnu_rel.items():   # Lnu_rel = {'ELL':..., 'SD':..., 'SB':...}
            # <Lnu>_band(z) を z_grid 上でだけ正確に計算
            Lnu_band_z = bandavg_Lnu_any_filter_exact(
                nu_tpl, Lnu_tpl, nu_flt, Tnu_flt, z_grid
            )
            nuLnu_band_z = bandavg_Lnu_any_filter_exact(
                nu_tpl, nu_tpl * Lnu_tpl, nu_flt, Tnu_flt, z_grid
            )
            precomputed_bandavg[fname][tpl_name] = {
                "Lnu":   Lnu_band_z,
                "nuLnu": nuLnu_band_z,
            }
    def Get_LF_all(filter_name):
        fdat = filters[filter_name]
        nu_flt   = fdat["nu"]
        Tnu_flt  = fdat["T"]
        denom_T  = fdat["denom_T"]

        # --- 連続光：z の補間を使う ---
        idx_Q  = np.where(mask_Q)[0]
        idx_Sd = np.where(mask_Sd)[0]
        idx_SB = np.where(mask_SB)[0]

        logLnu_band    = np.zeros(len(zobs), dtype=float)
        log_nuLnu_band = np.zeros(len(zobs), dtype=float)

        # それぞれのテンプレートについて，z_grid 上で precompute した値を
        # zobs に線形補間して使う
        if idx_Q.size:
            Lnu_z    = np.interp(zobs[idx_Q],    z_grid, precomputed_bandavg[filter_name]['ELL']['Lnu'])
            nuLnu_z  = np.interp(zobs[idx_Q],    z_grid, precomputed_bandavg[filter_name]['ELL']['nuLnu'])
            logLnu_band[idx_Q]    = np.log10(Lnu_z)   + log_scale[idx_Q]
            log_nuLnu_band[idx_Q] = np.log10(nuLnu_z) + log_scale[idx_Q]

        if idx_Sd.size:
            Lnu_z    = np.interp(zobs[idx_Sd],   z_grid, precomputed_bandavg[filter_name]['SD']['Lnu'])
            nuLnu_z  = np.interp(zobs[idx_Sd],   z_grid, precomputed_bandavg[filter_name]['SD']['nuLnu'])
            logLnu_band[idx_Sd]    = np.log10(Lnu_z)   + log_scale[idx_Sd]
            log_nuLnu_band[idx_Sd] = np.log10(nuLnu_z) + log_scale[idx_Sd]

        if idx_SB.size:
            Lnu_z    = np.interp(zobs[idx_SB],   z_grid, precomputed_bandavg[filter_name]['SB']['Lnu'])
            nuLnu_z  = np.interp(zobs[idx_SB],   z_grid, precomputed_bandavg[filter_name]['SB']['nuLnu'])
            logLnu_band[idx_SB]    = np.log10(Lnu_z)   + log_scale[idx_SB]
            log_nuLnu_band[idx_SB] = np.log10(nuLnu_z) + log_scale[idx_SB]

        # --- 連続光の fν, νfν ----
        log_fnu_cont = (
            logLnu_band
            - np.log10(4.0*np.pi)
            - 2.0*np.log10(lumi_dist)
            - np.log10(1.0 + zobs)
        )
        fnu_cont = 10**log_fnu_cont

        log_nufnu_cont = (
            log_nuLnu_band
            - np.log10(4.0*np.pi)
            - 2.0*np.log10(lumi_dist)
        )
        nufnu_cont = 10**log_nufnu_cont

        # ===== 輝線部分（ここは元のままで OK） =====
        lam_Ha_um   = c/ha/1e5
        lam_Pa_um   = c/pa/1e5
        lam_OIII_um = c/Oiii/1e5

        c_um_s = 2.99792458e14  # [μm/s]
        nu_Ha_rest   = c_um_s / lam_Ha_um
        nu_Pa_rest   = c_um_s / lam_Pa_um
        nu_OIII_rest = c_um_s / lam_OIII_um

        nu_Ha_obs   = nu_Ha_rest   / (1.0 + zobs)
        nu_Pa_obs   = nu_Pa_rest   / (1.0 + zobs)
        nu_OIII_obs = nu_OIII_rest / (1.0 + zobs)

        T_Ha   = np.interp(nu_Ha_obs,   nu_flt, Tnu_flt, left=0.0, right=0.0)
        T_Pa   = np.interp(nu_Pa_obs,   nu_flt, Tnu_flt, left=0.0, right=0.0)
        T_OIII = np.interp(nu_OIII_obs, nu_flt, Tnu_flt, left=0.0, right=0.0)

        F_Ha   = 10**log_flux_Ha
        F_Pa   = 10**log_flux_Pa
        F_OIII = 10**log_flux_OIII

        fnu_line = (F_Ha*T_Ha + F_Pa*T_Pa + F_OIII*T_OIII) / denom_T
        fnu_total = fnu_cont + fnu_line

        nufnu_line = (nu_Ha_obs * F_Ha * T_Ha +
                    nu_Pa_obs * F_Pa * T_Pa +
                    nu_OIII_obs * F_OIII * T_OIII) / denom_T
        nufnu_total = nufnu_cont + nufnu_line

        # --- AB mag & LF ---
        m_ab = -2.5*np.log10(np.clip(fnu_total, 1e-300, None)) - 48.60

        area_deg2 = ra_size/3600.0 * dec_size/3600.0
        bins_m = np.arange(15, 31, 0.5)
        hist_m, _ = np.histogram(m_ab, bins=bins_m)
        dNdm   = hist_m / area_deg2 / np.diff(bins_m)
        m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
        err_dNdm = np.sqrt(hist_m) / area_deg2 / np.diff(bins_m)

        plot_Nz_for_magnitude_bins(catalog.redshift_real, m_ab, label_prefix=f"{filter_name} distribution")


        data = np.loadtxt(f'../{filter_name}.txt', delimiter=',')
        plt.figure()
        plt.title(filter_name)
        if 'W' in filter_name:
            data_y = data[:,1] / (0.5)  # W1, W2 are per 0.5 mag
        else :
            data_y = data[:,1]
        plt.errorbar(m_cent, dNdm, yerr=err_dNdm, fmt='o', ms=3, label='mock (cont+lines)')
        plt.scatter(data[:,0], data_y, label='observation', color='red', s=10)
        plt.yscale('log'); plt.xlabel(f'm_AB ({filter_name})'); plt.ylabel('dN/dm [deg$^{-2}$ mag$^{-1}$]')
        plt.legend(); plt.tight_layout()

        return nufnu_total


    def _Get_LF_all(filter):
        if "Spitzer" in filter:
            if "Ch1" in filter:
                nu_flt, Tnu_flt = load_filter_throughput_txt(filter,
                    f"/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125ch1trans_full.txt")
            elif "Ch2" in filter:
                nu_flt, Tnu_flt = load_filter_throughput_txt(filter, 
                    f"/mnt/data_cat3/yuka/data/Spitzer_throughputs/201125ch2trans_full.txt")

        # --- フィルタ読み込み（ν昇順 & T(ν)） ---
        else:
            nu_flt, Tnu_flt = load_filter_throughput_txt(
                f"/mnt/data_cat3/yuka/data/nircam_throughputs/mean_throughputs/{filter}_May2024_mean_system_throughput.txt",
                wave_unit='um')
        denom_T = np.trapz(Tnu_flt, nu_flt)  # ∫T dν （分母）

        # ===== 連続光（既存と同じ：<Lν> -> fν_cont）=====
        idx_Q  = np.where(mask_Q)[0]
        idx_Sd = np.where(mask_Sd)[0]
        idx_SB = np.where(mask_SB)[0]

        logLnu_band = np.zeros(len(zobs), dtype=float)  # log10 <Lν>
        log_nuLnu_band = np.zeros(len(zobs), dtype=float)  # log10 <Lν>
        if idx_Q.size:
            logLnu_band[idx_Q] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, Lnu_rel['ELL'], nu_flt, Tnu_flt, zobs[idx_Q]
            )) + log_scale[idx_Q]
            log_nuLnu_band[idx_Q] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, nu_tpl*Lnu_rel['ELL'], nu_flt, Tnu_flt, zobs[idx_Q]
            )) + log_scale[idx_Q]
        if idx_Sd.size:
            logLnu_band[idx_Sd] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, Lnu_rel['SD'], nu_flt, Tnu_flt, zobs[idx_Sd]
            )) + log_scale[idx_Sd]
            log_nuLnu_band[idx_Sd] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, nu_tpl*Lnu_rel['SD'], nu_flt, Tnu_flt, zobs[idx_Sd]
            )) + log_scale[idx_Sd]
        if idx_SB.size:
            logLnu_band[idx_SB] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, Lnu_rel['SB'], nu_flt, Tnu_flt, zobs[idx_SB]
            )) + log_scale[idx_SB]
            log_nuLnu_band[idx_SB] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, nu_tpl*Lnu_rel['SB'], nu_flt, Tnu_flt, zobs[idx_SB]
            )) + log_scale[idx_SB]

        # 連続光の fν（cgs, erg s^-1 cm^-2 Hz^-1）
        log_fnu_cont = (
            logLnu_band
            - np.log10(4.0*np.pi)
            - 2.0*np.log10(lumi_dist)
            - np.log10(1.0 + zobs)
        )
        fnu_cont = 10**log_fnu_cont  # (Ngal,)

        log_nufnu_cont = (
            log_nuLnu_band
            - np.log10(4.0*np.pi)
            - 2.0*np.log10(lumi_dist)
        )
        nufnu_cont = 10**log_nufnu_cont  # (Ngal,) erg/s/cm^2

        # ===== 輝線の寄与をバンド平均 fν に足す =====
        # rest 波長 [μm]（必要に応じて調整）
        lam_Ha_um   = c/ha/1e5
        lam_Pa_um   = c/pa/1e5
        lam_OIII_um = c/Oiii/1e5

        c_um_s = 2.99792458e14  # [μm/s]
        nu_Ha_rest   = c_um_s / lam_Ha_um
        nu_Pa_rest   = c_um_s / lam_Pa_um
        nu_OIII_rest = c_um_s / lam_OIII_um

        # 観測周波数（各銀河）
        nu_Ha_obs   = nu_Ha_rest   / (1.0 + zobs)
        nu_Pa_obs   = nu_Pa_rest   / (1.0 + zobs)
        nu_OIII_obs = nu_OIII_rest / (1.0 + zobs)

        # フィルタ透過 T(ν_obs) を補間
        T_Ha   = np.interp(nu_Ha_obs,   nu_flt, Tnu_flt, left=0.0, right=0.0)
        T_Pa   = np.interp(nu_Pa_obs,   nu_flt, Tnu_flt, left=0.0, right=0.0)
        T_OIII = np.interp(nu_OIII_obs, nu_flt, Tnu_flt, left=0.0, right=0.0)

        # フラックス（erg/s/cm^2）
        F_Ha   = 10**log_flux_Ha
        F_Pa   = 10**log_flux_Pa
        F_OIII = 10**log_flux_OIII

        # バンド平均 fν_line = F_line * T(ν_obs) / ∫T dν
        fnu_line = (F_Ha*T_Ha + F_Pa*T_Pa + F_OIII*T_OIII) / max(denom_T, 1e-300)

        # 総フラックス密度（連続光＋輝線）
        fnu_total = fnu_cont + fnu_line

        nufnu_line = (nu_Ha_obs * F_Ha * T_Ha +
                  nu_Pa_obs * F_Pa * T_Pa +
                  nu_OIII_obs * F_OIII * T_OIII) / max(denom_T, 1e-300)

        nufnu_total = nufnu_cont + nufnu_line

        # AB 等級
        m_ab = -2.5*np.log10(np.clip(fnu_total, 1e-300, None)) - 48.60

        # === number counts（dN/dm）を作って比較 ===
        #area_deg2 = (np.pi * (params.radius**2))
        #area_deg2 = 4.0 * (params.radius**2)
        area_deg2 = ra_size/3600.0 * dec_size/3600.0
        bins_m = np.arange(15, 31, 0.5)
        hist_m, _ = np.histogram(m_ab, bins=bins_m)
        dNdm = hist_m / area_deg2 / np.diff(bins_m)

        m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
        err_dNdm = np.sqrt(hist_m) / area_deg2 / np.diff(bins_m)

        data = np.loadtxt(f'../{filter}.txt', delimiter=',')
        plt.figure()
        plt.title(filter)
        plt.errorbar(m_cent, dNdm, yerr=err_dNdm, fmt='o', ms=3, label='mock (cont+lines)')
        plt.scatter(data[:,0], data[:,1], label='observation', color='red', s=10)
        plt.yscale('log'); plt.xlabel(f'm_AB ({filter})'); plt.ylabel('dN/dm [deg$^{-2}$ mag$^{-1}$]')
        plt.legend(); plt.tight_layout()
        return nufnu_total
    
    results = [Get_LF_all(f) for f in filter_names]

    Ch1_flux = results[filter_names.index('SpitzerCh1')]
    Ch2_flux = results[filter_names.index('SpitzerCh2')]
 
    intensity_Ch1    = 1e6* Ch1_flux / (Ch1_resolution * arcsec)**2  # [nW/m^2/sr]
    intensity_Ch2    = 1e6* Ch2_flux / (Ch2_resolution * arcsec)**2  # [nW/m^2/sr]

    Ch1_image = np.zeros(Ch1_npix)
    Ch2_image = np.zeros(Ch2_npix)

    indices_Ch1 = np.array([Ch1_ix, Ch1_iy]).T
    valid_xy = (Ch1_ix >= 0) & (Ch1_ix < Ch1_Nx) & (Ch1_iy >= 0) & (Ch1_iy < Ch1_Nx)
    valid_mask = valid_xy
    
    indices_valid = indices_Ch1[valid_mask]
    
    np.add.at(Ch1_image,
              (indices_valid[:,0], indices_valid[:,1]),
              intensity_Ch1[valid_mask])
    
    indices_Ch2 = np.array([Ch2_ix, Ch2_iy]).T
    valid_xy = (Ch2_ix >= 0) & (Ch2_ix < Ch2_Nx) & (Ch2_iy >= 0) & (Ch2_iy < Ch2_Nx)
    valid_mask = valid_xy

    indices_valid = indices_Ch2[valid_mask]
    
    np.add.at(Ch2_image,
              (indices_valid[:,0], indices_valid[:,1]),
              intensity_Ch2[valid_mask])
    
    return Ch1_image, Ch2_image



def make_mock_old(params, name='Pinocchio', mag_lim=None):
    #SFR history, stellar mass function for redshift bins, JWST wide band luminosity function
    #spectrum of some galaxies
    flist, dflist = frequency(params) # Hz, Hz

    catalog = Catalog()
    if name=='Pinocchio':
        input_fname='/mnt/data_cat3/moriwaki/Pinocchio/my_lightcone/output_large/pinocchio.r00001.plc.out'
    elif name=='TNG':
        input_fname = '/mnt/data_cat3/yuka/data/lightcone_5400sec.txt'
    elif 'SIDES' in name:
        input_fname = '/mnt/data_cat3/yuka/data/SIDES/'

    catalog.loadcatalog(input_fname, name=name)
    lumi_dist = cosmo.luminosity_distance(np.array(catalog.redshift_real)).to(u.cm)
    lumi_dist = lumi_dist.value
    
    #ra = (np.array(catalog.x) - params.radius) * 3600.0 #[arcsec]
    #dec = (np.array(catalog.y) - params.radius) * 3600.0 #[arcsec]
    
    #ix = np.floor((ra + params.radius * 3600) / params.resolution).astype(np.int64)
    #iy = np.floor((dec + params.radius * 3600) / params.resolution).astype(np.int64)
    ra = np.array(catalog.x) * 3600.0 #[arcsec]
    ra -= np.min(ra)
    dec = np.array(catalog.y) * 3600.0 #[arcsec]
    dec -= np.min(dec)

    ix = np.floor(ra / params.resolution).astype(np.int64)
    iy = np.floor(dec / params.resolution).astype(np.int64)
    #print(ix)
    
    ##initialize
    #Nx = int(2.0*3600*params.radius / params.resolution)
    Nx = int(np.ceil( (np.max(ra) - np.min(ra)) / params.resolution ))
    Ny = int(np.ceil( (np.max(dec) - np.min(dec)) / params.resolution ))
    Nz = len(flist)
    
    npix = np.array([Nx, Ny, Nz])

    ############################Stellar continuum############################
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

    ########################### Ha intensity
    SFR = catalog.SFR
    Ha_intensity = np.zeros([Nx, Ny, Nz], dtype=np.float32)

    log_lumi_Ha = np.log10(SFR) + np.log10( 8.3e40 )
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_Ha))

    #lumi_Ha = 10**(log_lumi + log_lumi_scatter)
    log_flux_Ha = (log_lumi_Ha + log_lumi_scatter) - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_Ha = 10**log_flux_Ha

    ########################### Pa intensity
    Pa_intensity = np.zeros([Nx, Ny, Nz], dtype=np.float32)

    log_lumi_Pa = np.log10(SFR) + np.log10( 8.3e40 ) - 0.6
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_Pa))

    log_flux_Pa = (log_lumi_Pa + log_lumi_scatter)  - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_Pa = 10**log_flux_Pa

    ########################### [OIII] intensity
    OIII_intensity = np.zeros([Nx, Ny, Nz], dtype=np.float32)

    log_lumi_OIII = np.log10(SFR) + np.log10( 1.32e41 )
    log_lumi_scatter = 0.3 * np.random.randn(len(log_lumi_OIII))

    log_flux_OIII = (log_lumi_OIII + log_lumi_scatter)  - np.log10(4. * np.pi) -2* np.log10(lumi_dist) #[erg/s/cm^2]
    flux_OIII = 10**log_flux_OIII

    def Get_LF(filter):
        # --- フィルタ読み込み（ν昇順 & T(ν)） ---
        nu_flt, Tnu_flt = load_filter_throughput_txt(
            f"/mnt/data_cat3/yuka/data/nircam_throughputs/mean_throughputs/{filter}_May2024_mean_system_throughput.txt",
            wave_unit='um'
        )
        denom_T = np.trapz(Tnu_flt, nu_flt)  # ∫T dν （分母）

        # ===== 連続光（既存と同じ：<Lν> -> fν_cont）=====
        idx_Q  = np.where(mask_Q)[0]
        idx_Sd = np.where(mask_Sd)[0]
        idx_SB = np.where(mask_SB)[0]

        logLnu_band = np.zeros(len(zobs), dtype=float)  # log10 <Lν>
        if idx_Q.size:
            logLnu_band[idx_Q] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, Lnu_rel['ELL'], nu_flt, Tnu_flt, zobs[idx_Q]
            )) + log_scale[idx_Q]
        if idx_Sd.size:
            logLnu_band[idx_Sd] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, Lnu_rel['SD'], nu_flt, Tnu_flt, zobs[idx_Sd]
            )) + log_scale[idx_Sd]
        if idx_SB.size:
            logLnu_band[idx_SB] = np.log10(bandavg_Lnu_any_filter_exact(
                nu_tpl, Lnu_rel['SB'], nu_flt, Tnu_flt, zobs[idx_SB]
            )) + log_scale[idx_SB]

        # 連続光の fν（cgs, erg s^-1 cm^-2 Hz^-1）
        log_fnu_cont = (
            logLnu_band
            - np.log10(4.0*np.pi)
            - 2.0*np.log10(lumi_dist)
            - np.log10(1.0 + zobs)
        )
        fnu_cont = 10**log_fnu_cont  # (Ngal,)

        # ===== 輝線の寄与をバンド平均 fν に足す =====
        # rest 波長 [μm]（必要に応じて調整）
        lam_Ha_um   = c/ha/1e5
        lam_Pa_um   = c/pa/1e5
        lam_OIII_um = c/Oiii/1e5

        c_um_s = 2.99792458e14  # [μm/s]
        nu_Ha_rest   = c_um_s / lam_Ha_um
        nu_Pa_rest   = c_um_s / lam_Pa_um
        nu_OIII_rest = c_um_s / lam_OIII_um

        # 観測周波数（各銀河）
        nu_Ha_obs   = nu_Ha_rest   / (1.0 + zobs)
        nu_Pa_obs   = nu_Pa_rest   / (1.0 + zobs)
        nu_OIII_obs = nu_OIII_rest / (1.0 + zobs)

        # フィルタ透過 T(ν_obs) を補間
        T_Ha   = np.interp(nu_Ha_obs,   nu_flt, Tnu_flt, left=0.0, right=0.0)
        T_Pa   = np.interp(nu_Pa_obs,   nu_flt, Tnu_flt, left=0.0, right=0.0)
        T_OIII = np.interp(nu_OIII_obs, nu_flt, Tnu_flt, left=0.0, right=0.0)

        # フラックス（erg/s/cm^2）
        F_Ha   = 10**log_flux_Ha
        F_Pa   = 10**log_flux_Pa
        F_OIII = 10**log_flux_OIII

        # バンド平均 fν_line = F_line * T(ν_obs) / ∫T dν
        fnu_line = (F_Ha*T_Ha + F_Pa*T_Pa + F_OIII*T_OIII) / max(denom_T, 1e-300)

        # 総フラックス密度（連続光＋輝線）
        fnu_total = fnu_cont + fnu_line

        # AB 等級
        m_ab = -2.5*np.log10(np.clip(fnu_total, 1e-300, None)) - 48.60

        # === number counts（dN/dm）を作って比較 ===
        #area_deg2 = (np.pi * (params.radius**2))
        #area_deg2 = 4.0 * (params.radius**2)
        area_deg2 = (Nx * params.resolution / 3600.0) * (Ny * params.resolution / 3600.0)
        bins_m = np.arange(15, 31, 0.5)
        hist_m, _ = np.histogram(m_ab, bins=bins_m)
        dNdm = hist_m / area_deg2 / np.diff(bins_m)

        m_cent = 0.5*(bins_m[:-1] + bins_m[1:])
        err_dNdm = np.sqrt(hist_m) / area_deg2 / np.diff(bins_m)

        data = np.loadtxt('../J150W.txt', delimiter=',')
        plt.figure()
        plt.errorbar(m_cent, dNdm, yerr=err_dNdm, fmt='o', ms=3, label='mock (cont+lines)')
        plt.scatter(data[:,0], data[:,1], label='observation', color='red', s=10)
        plt.yscale('log'); plt.xlabel(f'm_AB ({filter})'); plt.ylabel('dN/dm [deg$^{-2}$ mag$^{-1}$]')
        plt.legend(); plt.tight_layout()
        return m_ab


    m_ab = Get_LF('F150W')
    if mag_lim==None:
        limit = (m_ab < float('inf'))
    else:
        limit = (m_ab < mag_lim)

    #######################################################################################
    #####add continuum###########################
    valid_xy = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Ny) & limit

    # 🔴 valid銀河だけ抽出して add.at する
    fnu_valid = fnu_jy[valid_xy, :]                               # (Nvalid, Nz)

    # --- 3D accumulation: (x,y,ν) ---
    # add.at は重複 index の和を正しくとってくれる
    ix_v = ix[valid_xy]
    iy_v = iy[valid_xy]

    continuum = np.zeros([Nx, Ny, Nz], dtype=np.float32)
    np.add.at(continuum, (ix_v, iy_v, slice(None)), fnu_valid)
    
    ###############################################
    fedges = np.empty(len(flist) + 1, dtype=float)
    fedges[1:-1] = 0.5*(flist[:-1] + flist[1:])
    fedges[0]    = flist[0]  - 0.5*(dflist[0])
    fedges[-1]   = flist[-1] + 0.5*(dflist[-1])

    def bin_index_from_freq_edges(fedges, nu_obs):
        iz = np.searchsorted(fedges, nu_obs, side='right') - 1
        return np.clip(iz, 0, len(fedges)-2)
    ######add Ha##################################
    freq_obs = ha/(1 + np.array(zobs))*1e9 #[Hz]
    iz_Ha = bin_index_from_freq_edges(fedges, freq_obs)

    indices = np.array([ix, iy, iz_Ha]).T
    valid_xy = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Ny) & (iz_Ha >= 0) & (iz_Ha < Nz)
    is_SF = SFR > 0
    valid_mask = valid_xy & is_SF & limit
    
    indices_valid = indices[valid_mask]
    flux_Ha_valid   = flux_Ha[valid_mask]
    dnu_valid       = dflist[indices_valid[:, 2]]
    intensity_Ha    = flux_Ha_valid / dnu_valid / Jy / (params.resolution * arcsec)**2  # [Jy/sr]
    intensity_Ha_pergal = flux_Ha / dflist[iz_Ha] / Jy / (params.resolution * arcsec)**2 #[Jy/sr]
    
    np.add.at(Ha_intensity,
              (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]),
              intensity_Ha)
    ########add Pa#############################
    freq_obs = pa/(1 + np.array(zobs))*1e9 #[Hz]
    iz_Pa = bin_index_from_freq_edges(fedges, freq_obs)
    
    indices = np.array([ix, iy, iz_Pa]).T
    valid_xy = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Ny) & (iz_Pa >= 0) & (iz_Pa < Nz)
    valid_mask = valid_xy & is_SF & limit
    
    indices_valid = indices[valid_mask]
    flux_Pa_valid   = flux_Pa[valid_mask]
    dnu_valid       = dflist[indices_valid[:, 2]]
    intensity_Pa    = flux_Pa_valid / dnu_valid / Jy / (params.resolution * arcsec)**2  # [Jy/sr]
    intensity_Pa_pergal = flux_Pa / dflist[iz_Pa] / Jy / (params.resolution * arcsec)**2 #[Jy/sr]

    np.add.at(Pa_intensity,
              (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]),
              intensity_Pa)
    ########add OIII############################
    freq_obs = Oiii/(1 + np.array(zobs)) * 1e9 #[Hz]
    iz_OIII = bin_index_from_freq_edges(fedges, freq_obs)
    
    indices = np.array([ix, iy, iz_OIII]).T
    valid_xy = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Ny) & (iz_OIII >= 0) & (iz_OIII < Nz)
    valid_mask = valid_xy & is_SF & limit
    
    indices_valid   = indices[valid_mask]
    flux_OIII_valid = flux_OIII[valid_mask]
    dnu_valid       = dflist[indices_valid[:, 2]]
    intensity_OIII  = flux_OIII_valid / dnu_valid / Jy / (params.resolution * arcsec)**2
    intensity_OIII_pergal = flux_OIII / dflist[iz_OIII] / Jy / (params.resolution * arcsec)**2 #[Jy/sr]

    np.add.at(OIII_intensity,
              (indices_valid[:,0], indices_valid[:,1], indices_valid[:,2]),
              intensity_OIII)

    total_intensity = Ha_intensity + OIII_intensity + Pa_intensity + continuum
    
    #outdir = "../output/catalog"
    outdir = input_fname + 'spectra/' + name
    
    save_meta_parquet(
        outdir, catalog, ix, iy, limit,
        iz_Ha, iz_Pa, iz_OIII,
        intensity_Ha_pergal, intensity_Pa_pergal, intensity_OIII_pergal,
        flist, dflist
        )

    # 2) 輝線スペクトル 2D を作成（点置き or LSF分配）
    fnu_lines_jy = build_line_spectra_jy(
        Ngal=fnu_jy.shape[0], flen=fnu_jy.shape[1],
        iz_Ha=iz_Ha, I_Ha=intensity_Ha_pergal,
        iz_Pa=iz_Pa, I_Pa=intensity_Pa_pergal,
        iz_OIII=iz_OIII, I_OIII=intensity_OIII_pergal,
        gaussian_sigma_bins=None        # 例：None=点置き、0.7等でガウシアン分配
        )

    # 3) HDF5 保存（連続光＋輝線＋合成）
    save_spectra_h5(outdir, flist, dflist, fnu_cont_jy=fnu_jy, fnu_lines_jy=fnu_lines_jy)
    
    return total_intensity, Ha_intensity, OIII_intensity, Pa_intensity, continuum

def save_meta_parquet(outdir, catalog, ix, iy, limit,
                      iz_Ha, iz_Pa, iz_OIII,
                      I_Ha, I_Pa, I_OIII, flist, dflist):
    os.makedirs(outdir, exist_ok=True)

    Mhalo = 10**np.asarray(catalog.logm)

    df = pd.DataFrame({
        "gal_id"       : np.arange(len(catalog.Mstar), dtype=np.int64),
        "x"            : np.asarray(catalog.x),
        "y"            : np.asarray(catalog.y),
        "z_obs"        : np.asarray(catalog.redshift_obs),
        "z_real"       : np.asarray(catalog.redshift_real),
        "Mhalo"        : np.asarray(Mhalo),
        "Mstar"        : np.asarray(catalog.Mstar),
        "SFR"          : np.asarray(catalog.SFR),
        "JLumi"        : np.asarray(catalog.JLumi),
        "is_quenched"  : np.asarray(catalog.quenched).astype(bool),
        "is_starburst" : np.asarray(catalog.Starburst).astype(bool),
        "ix"           : np.asarray(ix, dtype=np.int32),
        "iy"           : np.asarray(iy, dtype=np.int32),
        "in_limit"     : np.asarray(limit, dtype=bool),
        # line info per galaxy（落ちたbinと面輝度）
        "Ha_iz"        : np.asarray(iz_Ha,   dtype=np.int32),
        "Ha_I_Jy_sr"   : np.asarray(I_Ha,    dtype=np.float32),
        "Pa_iz"        : np.asarray(iz_Pa,   dtype=np.int32),
        "Pa_I_Jy_sr"   : np.asarray(I_Pa,    dtype=np.float32),
        "OIII_iz"      : np.asarray(iz_OIII, dtype=np.int32),
        "OIII_I_Jy_sr" : np.asarray(I_OIII,  dtype=np.float32),
        # 参照
        "Nz"           : int(len(flist)),
        "nu_min_Hz"    : float(flist.min()),
        "nu_max_Hz"    : float(flist.max()),
    })
    path = outdir + "_galaxies_meta.parquet"
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
    path = outdir + "_spectra.h5"
    #path = os.path.join(outdir, "spectra_TNG.h5")
    with h5py.File(path, "w") as h5:
        h5.create_dataset("flist_Hz",  data=flist,  compression="gzip")
        h5.create_dataset("dflist_Hz", data=dflist, compression="gzip")
        h5.create_dataset(
            "fnu_cont_jy", data=fnu_cont_jy.astype(np.float32),
            compression="gzip", compression_opts=4,
            chunks=(min(1024, fnu_cont_jy.shape[0]), min(2048, fnu_cont_jy.shape[1]))
        )
        if fnu_lines_jy is not None:
            h5.create_dataset(
                "fnu_lines_jy", data=fnu_lines_jy.astype(np.float32),
                compression="gzip", compression_opts=4,
                chunks=(min(1024, fnu_lines_jy.shape[0]), min(2048, fnu_lines_jy.shape[1]))
            )
            fnu_total = fnu_cont_jy + fnu_lines_jy
            h5.create_dataset(
                "fnu_total_jy", data=fnu_total.astype(np.float32),
                compression="gzip", compression_opts=4,
                chunks=(min(1024, fnu_total.shape[0]), min(2048, fnu_total.shape[1]))
            )

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

    return wave, log_spec

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
        
def _main():
    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.SHUFFLE)

    params = Params()
    params.Default()
    params.frequency = [60000, 400000]
    params.radius = 9.0
    params.resolution = 6.5

    flist, dflist = frequency(params)

    init_worker_fsps()
    
    meta = pd.read_parquet("/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES_galaxies_all_meta.parquet")
    idx = meta["gal_id"].to_numpy()
    ix = meta["ix"].to_numpy()
    iy = meta["iy"].to_numpy()
    redshift = meta["z_obs"].to_numpy()
    
    Ha_iz = meta["Ha_iz"].to_numpy()
    Pa_iz = meta["Pa_iz"].to_numpy()
    OIII_iz = meta["OIII_iz"].to_numpy()
    
    Ha_F = meta["Ha_F"].to_numpy()
    Pa_F = meta["Pa_F"].to_numpy()
    OIII_F = meta["OIII_F"].to_numpy()
    
    unique_idx = np.unique(idx)
    Ntemplate_used = len(unique_idx)
    print(f"number of template is {Ntemplate_used}")
    
    f_sample = "/mnt/data_cat3/yuka/repository/LIM_mock/theta_samples_100000_seed1_mag.npz"
    data = np.load(f_sample, allow_pickle=True)
    samples = data["samples"]
    logM_samples = data["logM_samples"]
    logM_formed = data['logM_formed']
    tage = data['tage']
    #sfr = data['SFR']
    log10SFR_bins = data['log10SFR_bins']
    t_edge_bins = data['t_edge_bins']
    redshift_bin = samples[:,-1]

    print("Load {}".format(f_sample))

    worker_iter = zip(samples[unique_idx], logM_formed[unique_idx], tage[unique_idx], log10SFR_bins[unique_idx], t_edge_bins[unique_idx])
    with Pool(processes=20) as pool:
        results = pool.map(get_true_spectrum, worker_iter)
        #results = pool.map(get_true_spectrum, samples)
    waves, log_abs_spectrum = zip(*results)
    wavelength = np.array(waves[0])
    log_abs_spectrum = np.array(log_abs_spectrum) # (N_samples, N_lambda) rest frame luminosity[erg/s/Hz]
    
    Nx = int(params.radius * 3600.0 / params.resolution)
    Npix = Nx*Nx
    Nspec = len(flist)
    
    z = zarr.open(
        "/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES/continuum_all.zarr",
        mode="w",
        shape=(Nx, Nx, Nspec),
        chunks=(256, 256, 1),
        dtype="float32",
        compressor=compressor,
        zarr_format=2
    )
    
    #400のredshift binに分ける
    Nbin = 400
    z_edges = np.linspace(0, 4, Nbin)
    redshift_idz =  np.digitize(redshift, z_edges) - 1
    redshift_idz = np.clip(redshift_idz, 0, Nbin - 1)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    
    Nbin_wide = 20
    z_edge_wide = np.linspace(0, 4, Nbin_wide)
    for i in range(Nbin_wide - 1):
        continuum_flat = np.zeros((Npix, Nspec), dtype=np.float32)
        
        zmin = z_edge_wide[i]
        zmax = z_edge_wide[i+1]
        
        selection = (redshift>=zmin)&(redshift<zmax)&(ix>=0)&(ix<Nx)&(iy>=0)&(iy<Nx)
        
        target_idz = redshift_idz[selection] #赤方偏移のbin
        target_idx = idx[selection] #template bin
        target_ix = ix[selection]
        target_iy = iy[selection]
        
        pairs = np.stack([target_idz, target_idx], axis=1)
        pairs_unique, inv = np.unique(pairs, axis=0, return_inverse=True)
        templates = log_abs_spectrum[pairs_unique[:,1]]
        redshifts = z_centers[pairs_unique[:,0]]
        
        obs_spectrum = bandavg_group_top_hat(wavelength, templates, flist, dflist, redshifts)
        
        pix = target_ix + Nx * target_iy
        
        for ch in range(Nspec):
            continuum_flat[:, ch] += np.bincount(
                pix,
                weights=obs_spectrum[inv, ch],
                minlength=Npix
                )
    
        continuum_box = continuum_flat.reshape(Nx, Nx, Nspec)
        
        z[:, :, :] += continuum_box

        # メモリ解放
        del continuum_flat, continuum_box
    #########################あとは輝線################################
    Ha_z = zarr.open(
        "/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES/Ha.zarr",
        mode="w",
        shape=(Nx, Nx, Nspec),
        chunks=(256, 256, 1),
        dtype="float32",
        compressor=compressor,
        zarr_format=2
    )

    Pa_z = zarr.open(
        "/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES/Pa.zarr",
        mode="w",
        shape=(Nx, Nx, Nspec),
        chunks=(256, 256, 1),
        dtype="float32",
        compressor=compressor,
        zarr_format=2
    )

    OIII_z = zarr.open(
        "/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES/OIII.zarr",
        mode="w",
        shape=(Nx, Nx, Nspec),
        chunks=(256, 256, 1),
        dtype="float32",
        compressor=compressor,
        zarr_format=2
    )
    ##########Ha###########
    sel = (Ha_iz >= 0) & (Ha_iz < Nspec)&(ix>=0)&(ix<Nx)&(iy>=0)&(iy<Nx)
    Ha_z[ix[sel],iy[sel],Ha_iz[sel]] += Ha_F[sel].astype(np.float32)
    
    selection = (Pa_iz>=0)&(Pa_iz<Nspec)&(ix>=0)&(ix<Nx)&(iy>=0)&(iy<Nx)
    Pa_z[ix[selection], iy[selection], Pa_iz[selection]] += Pa_F[selection]
    
    selection = (OIII_iz>=0)&(OIII_iz<Nspec)&(ix>=0)&(ix<Nx)&(iy>=0)&(iy<Nx)
    OIII_z[ix[selection], iy[selection], OIII_iz[selection]] += OIII_F[selection]
    
def main():
    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.SHUFFLE)

    params = Params()
    params.Default()
    params.frequency = [60000, 400000]
    params.radius = 9.0
    params.resolution = 6.5

    flist, dflist = frequency(params)

    init_worker_fsps()
    
    meta = pd.read_parquet("/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES_galaxies_all_meta.parquet")
    idx = meta["gal_id"].to_numpy()
    ix = meta["ix"].to_numpy()
    iy = meta["iy"].to_numpy()
    redshift = meta["z_obs"].to_numpy()
    
    Ha_iz = meta["Ha_iz"].to_numpy()
    Pa_iz = meta["Pa_iz"].to_numpy()
    OIII_iz = meta["OIII_iz"].to_numpy()
    
    Ha_F = meta["Ha_F"].to_numpy()
    Pa_F = meta["Pa_F"].to_numpy()
    OIII_F = meta["OIII_F"].to_numpy()
    
    Nx = int(params.radius * 3600.0 / params.resolution)
    Npix = Nx*Nx
    Nspec = len(flist)
    
    Ha_z = zarr.open(
        "/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES/Ha.zarr",
        mode="w",
        shape=(Nx, Nx, Nspec),
        chunks=(256, 256, 1),
        dtype="float32",
        compressor=compressor,
        zarr_format=2
    )

    Pa_z = zarr.open(
        "/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES/Pa.zarr",
        mode="w",
        shape=(Nx, Nx, Nspec),
        chunks=(256, 256, 1),
        dtype="float32",
        compressor=compressor,
        zarr_format=2
    )

    OIII_z = zarr.open(
        "/mnt/data_cat3/yuka/repository/LIM_mock/output/catalog/spectra/SIDES/OIII.zarr",
        mode="w",
        shape=(Nx, Nx, Nspec),
        chunks=(256, 256, 1),
        dtype="float32",
        compressor=compressor,
        zarr_format=2
    )
    
    sel = (Ha_iz >= 0) & (Ha_iz < Nspec)&(ix>=0)&(ix<Nx)&(iy>=0)&(iy<Nx)
    Ha_z[ix[sel],iy[sel],Ha_iz[sel]] += Ha_F[sel].astype(np.float32)
    
    selection = (Pa_iz>=0)&(Pa_iz<Nspec)&(ix>=0)&(ix<Nx)&(iy>=0)&(iy<Nx)
    Pa_z[ix[selection], iy[selection], Pa_iz[selection]] += Pa_F[selection]
    
    selection = (OIII_iz>=0)&(OIII_iz<Nspec)&(ix>=0)&(ix<Nx)&(iy>=0)&(iy<Nx)
    OIII_z[ix[selection], iy[selection], OIII_iz[selection]] += OIII_F[selection]
    
    

if __name__ == "__main__":
    main()
    