"""
Load text data like this

---
# redshift: 0.2 - 12.0
# tan_LoS: 0.06 0.10
# area: 5400.0 arcsec x 5400.0 arcsec
# theta_x[arcsec] theta_y[arcsec] redshift log_lumi_dis vz[km/s] rmass_half log_sfr log_mstar log_Z
5.062337e+03 2.392666e+02 0.1893 27.4628 2.098443e+02 1.197039e+00 0.7396 10.3443 -1.7383
---

"""

cspeed = 2.9979e10  # speed of light in cm/s
Jy = 1.0e-23        # jansky (erg/s/cm2/Hz)
arcsec = 4.848136811094e-6 # [rad] ... arcmin / 60 //


import numpy as np
import os
import argparse
from astropy.cosmology import FlatLambdaCDM, z_at_value
from astropy.cosmology import Planck15
import astropy.units as u
import scipy
from multiprocessing import Pool
from functools import partial
import pandas as pd

import h5py

from scipy.interpolate import interp1d
from scipy.special import erf

from line_model import calc_line_luminosity, line_dict
class Args(object):
    def __init__(self):
        self.sigma = 0.20

class Param(object):
    def __init__(self):
        self.redshift = None
        self.BoxLen = None
        self.cosmo = None

class Boxes(object):
    def __init__(self):
        self.position = None
        self.velocity = None
        self.log_mass = None
        self.log_Mstar = None
        self.log_SFR = None
        self.Z = None
        self.param = None

    def loadUCHUU(self, path, redshift):
        param = Param()
        param.redshift = redshift
        param.BoxLen = 400.0
        param.cosmo = Planck15
        _redshift = str(redshift).replace('.', 'p')
        filename = f"MiniUchuu_halolist_z{_redshift}.h5"
        self.redshift = redshift
        with h5py.File(os.path.join(path, filename), "r") as f:
            x = f["x"][:]
            y = f["y"][:]
            z = f["z"][:]

            vx = f["vx"][:]
            vy = f["vy"][:]
            vz = f["vz"][:]

            Mvir = f["Mvir"][:]

        self.position = np.array([x, y, z]).T
        self.velocity = np.array([vx, vy, vz]).T
        self.mass = Mvir
        self.param = param
    
    def loadTNG(self, path, redshift):
        param = Param()
        param.redshift = redshift
        param.BoxLen = 205.0
        param.cosmo = Planck15
        
        _redshift = str(redshift).replace('.', 'p')
        filename = f"TNG-300_z{_redshift}.txt"
        data = np.loadtxt(os.path.join(path, filename), skiprows=4)
        
        x = data[:,0]
        y = data[:,1]
        z = data[:,2]
        
        vx = data[:,3]
        vy = data[:,4]
        vz = data[:,5]
        
        log_sfr = data[:,7]
        log_mstar = data[:,8]
        Z = data[:,9]
        
        log_mass = data[:,10]
        
        self.position = np.array([x, y, z]).T
        self.velocity = np.array([vx, vy, vz]).T
        
        self.log_mass = log_mass
        self.log_Mstar = log_mstar
        self.log_SFR = log_sfr
        self.Z = Z
        self.param = param
        
    def loadTNG_Hydro(self, path, redshift):
        param = Param()
        param.redshift = redshift
        param.BoxLen = 205.0
        param.cosmo = Planck15
        
class Slices(object):
    def __init__(self):
        self.ra = None
        self.dec = None
        self.z = None
        self.obs_z = None
        self.log_lumi_dist = None
        self.intensity = None
    
def get_slices(boxes, line_list, count, R=41.0, area=5.0):
    redshift = boxes.param.redshift
    BoxLen = boxes.param.BoxLen
    cosmo = boxes.param.cosmo
    h = cosmo.H(0).value/100.0
        
    z_samples = np.linspace(0, 2.0, 2000)  
    d_samples = cosmo.comoving_distance(z_samples).value*h  # Mpc/h 単位
    d_z = interp1d(d_samples, z_samples, kind='cubic', fill_value="extrapolate")
        
    min_redshift = redshift - (1.0 + redshift)/2.0/R
    max_redshift = redshift + (1.0 + redshift)/2.0/R
    min_dist = cosmo.comoving_distance(min_redshift).value*h  # Mpc/h 単位
    max_dist = cosmo.comoving_distance(max_redshift).value*h  # Mpc/h 
    slice_width = max_dist - min_dist
    print(f'slice width for redshift {redshift} will be {slice_width} Mpc/h')
    if (count*slice_width>BoxLen):
        print('There will be overlaps in the galaxies')
        
    index = np.arange(count)
        
    for line in line_list:
        rest_frequency = line_dict[line][0]
        rest_wavelength = line_dict[line][1]
        with Pool(processes=count) as pool:
            func = partial(split, boxes=boxes, R=R,  count=count, d_z=d_z, line=line, area=area)
            results = pool.map(func, index)
        return results
    
def split(index, boxes, R, count, d_z, line, area):
    param = boxes.param
    redshift = param.redshift
    cosmo = param.cosmo
    h = cosmo.H(0).value/100.0
    BoxLen = param.BoxLen
    offset = index * BoxLen/count
  
    log_mass = boxes.log_mass
    log_Mstar = boxes.log_Mstar
    log_SFR = boxes.log_SFR
    Z = boxes.Z
    
    min_redshift = redshift - (1.0 + redshift)/2.0/R
    max_redshift = redshift + (1.0 + redshift)/2.0/R
    
    ra, dec, obs_redshift, true_redshift, d_L = get_coordinate(boxes, d_z, offset, pec_vel=True)
    
    obs = (min_redshift<obs_redshift)&(max_redshift>obs_redshift)&(ra>-area/2.0)&(ra<area/2.0)&(dec>-area/2.0)&(dec<area/2.0)
    
    slice = Slices()
    slice.ra = ra[obs]
    slice.dec = dec[obs]
    slice.z = true_redshift[obs]
    slice.obs_z = obs_redshift[obs]
    slice.log_lumi_dist = np.log10(d_L[obs])
    
    if (log_Mstar is None):
        log_Mstar = get_Smass(true_redshift[obs], log_mass[obs])
        
    if (log_SFR is None):
        log_sfr = get_SFR(true_redshift[obs], log_Mstar)
        
    args = Args()
    luminosity = calc_line_luminosity(args, log_sfr, log_Mstar, line)
    
    slice.intensity = {line:luminosity}
        
    
        
    return slice

def get_coordinate(boxes, d_z, offset, pec_vel=True):
    param = boxes.param
    redshift = param.redshift
    cosmo = param.cosmo
    h = cosmo.H(0).value/100.0
    BoxLen = param.BoxLen
    
    distance_box = cosmo.comoving_distance(redshift).value*h
    
    x = boxes.position[:,0]
    x -= (x>BoxLen/2.0)*BoxLen
    
    y = boxes.position[:,1]
    y -= (y>BoxLen/2.0)*BoxLen
    
    z = boxes.position[:,2]
    z = z - offset
    z = z - BoxLen*(z>BoxLen/2.0) + distance_box
    
    vx = boxes.velocity[:,0]
    vy = boxes.velocity[:,1]
    vz = boxes.velocity[:,2]
    
    #convert z coordinate to redshift
    print(f"distance to center is {distance_box} Mpc/h")
    distance = (x**2 + y**2 + z**2)**0.5 #Mpc/h
    redshift = d_z(distance)
    
    d_L = cosmo.luminosity_distance(redshift).value
    
    ra = np.degrees(np.arctan2(x, distance))
    dec = np.degrees(np.arctan2(y, distance))

    c = cspeed /1e2/1e3
    
    if (pec_vel):
        v_parallel = (x*vx + y*vy + (z+distance_box)*vz)/distance
        obs_redshift =  (1 + redshift) * (1 + v_parallel / c) - 1
        
    return ra, dec, obs_redshift, redshift, d_L

def get_Smass(redshift, log_mass):
    base_dir = os.path.dirname(__file__)  # この.pyファイルのディレクトリ
    full_path = os.path.join(base_dir, 'params/STARS.par')
    star_params = load_params(full_path)
    
    log_MA = star_params['B'] + redshift*star_params['mu']
    A = star_params['C'] * (1+redshift)**star_params['nu']
    gamma = star_params['D'] * (1+redshift)**star_params['eta']
    beta = star_params['F'] * redshift + star_params['E']
    log_hA = log_mass - log_MA
    hA = 10.0**log_hA
    
    log_ratio = np.log10(2.0*A) -np.log10(hA**(-beta) + hA**gamma)
    log_smass = log_mass + log_ratio
    return log_smass

def get_SFR(redshift, log_smass):
    base_dir = os.path.dirname(__file__)  # この.pyファイルのディレクトリ
    full_path = os.path.join(base_dir, 'params/SIDES.par')
    params = load_params(full_path)
    print('Generate the star-formation properties...')

    print('Draw quenched galaxies...')

    Ngal = len(redshift)

    #Draw randomly which galaxies are quenched using the recipe from Bethermin+17
    Mtz = params['Mt0'] + params['alpha1'] * redshift + params['alpha2'] * redshift**2
    sigmaz =  params['sigma0'] +  params['beta1'] * redshift +  params['beta2'] * redshift**2
    qfrac0z = params['qfrac0'] * (1.+redshift)**params['gamma']
    
    Prob_SF = (1.-qfrac0z) * 0.5 * (1. - erf( ( log_smass - Mtz) /sigmaz ) )

    Xuni = np.random.rand(Ngal)

    qflag = Xuni > Prob_SF
    #Generate SFR for non-quenched objects

    print('Generate SFRs...')

    index_SF = np.where(qflag == False)

    m_SF = np.array(log_smass[index_SF[0]]+ np.log10(params['Chab2Salp'] / 1.e9 ))
    z_SF = np.array(redshift[index_SF[0]])
    r = np.log10(1.+z_SF)
    expr = np.maximum(m_SF - params['m1'] - params['a2'] * r, 0.)

    logSFRms_SF = m_SF - params['m0'] + params['a0'] * r - params['a1'] * expr**2 - np.log10(params['Chab2Salp'])

    logSFRms_SF += params['corr_zmean_lowzcorr'] * (params['zmax_lowzcorr'] - np.minimum(z_SF, params['zmax_lowzcorr'])) / (params['zmax_lowzcorr'] - params['zmean_lowzcorr'])

    Psb = params['Psb_hz'] + params['slope_Psb'] * (params['z_Psb_knee'] - np.minimum(z_SF, params['z_Psb_knee']))

    Xuni = np.random.rand(np.size(Psb))

    issb = (Xuni < Psb )

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
    return SFR
    


        
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
        
    

            
            
        

parser = argparse.ArgumentParser(description="Create mock data for SUBLIME")

parser.add_argument("--seed", type=int, default=123, help="random seed")

### I/O parameters
parser.add_argument("--path", type=str, default="./output/lightcone.txt")
parser.add_argument("--output_dir", type=str, default="output")

### area and angular resolution
parser.add_argument("--side_length", type=float, default=100.0, help="side length in arcsec")
parser.add_argument("--angular_resolution", type=float, default=30, help="angular resolution in arcsec")

### spectral range and resolution
parser.add_argument("--fmin", type=float, default=10.0, help="minimum frequency in GHz")
parser.add_argument("--fmax", type=float, default=100.0, help="maximum frequency in GHz")
parser.add_argument("--R", type=float, default=100, help="spectral resolution R")

### line model paramters 
parser.add_argument("--z_dependence_alpha", type=float, default=1.0, help="The SFR at very high redshift is multiplied by this alpha.")
parser.add_argument("--z_middle", type=float, default=8.0, help="If z_dependence_alpha != 1, the SFR deviates from the original value as redshift becomes higher; z_middle is the redshift when the factor is 0.5(alpha + 1)")
parser.add_argument("--sigma", type=float, default=0.2, help="The scatter in log [dex] added to luminosity - SFR relation")

### galaxy catalog parameters
parser.add_argument("--threshold", type=float, default=-1, help="intensity map is generated when threshold < 0. If threshold > 0, galaxy catalogs with flux [erg/s/cm2] > threshold is generated.")

args = parser.parse_args()

np.random.seed(args.seed)   

print("# frequency: {:.4f} - {:.4f} [GHz]".format(args.fmin, args.fmax))
print("# spectral resolution R: {:.4f}".format(args.R))
print("# area : {:.4f} arcsec x {:.4f} arcsec".format(args.side_length, args.side_length))
print("# angular resolution : {:.4f} arcsec".format(args.angular_resolution))

if args.threshold < 0:
    print("# Create intensity map")
else:
    print("# Create galaxy number density map with flux > {:e} erg/s/cm2".format(args.threshold))

### Define the spectral bins
flist = []
dflist = []

fnow = args.fmin * 1e9
while fnow <= args.fmax * 1e9:
    flist.append(fnow)
    dflist.append(fnow / args.R)
    fnow += fnow / args.R
flist = np.array(flist, dtype=np.float32)
dflist = np.array(dflist, dtype=np.float32)
print(f"frequency range is {flist[0]/1e9} - {flist[-1]/1e9} [GHz]")

### Load data
print("# Load data from {}".format(args.path))

df = pd.read_csv(args.path, comment="#", header=None, delim_whitespace=True, dtype=np.float32)
lc_data = df.values
x = lc_data[:, 0]
y = lc_data[:, 1]
z = lc_data[:, 2]
log_lumi_dis = lc_data[:, 3]
vz = lc_data[:, 4] * 1e5 # [cm/s]
rmass_half = lc_data[:, 5]
log_sfr = lc_data[:, 6]
log_mstar = lc_data[:, 7]
log_Z = lc_data[:, 8]

ix = np.floor(x / args.angular_resolution).astype(np.int32)
iy = np.floor(y / args.angular_resolution).astype(np.int32)

print("# Number of galaxies: {}".format(len(x)))

z_obs = z + vz / cspeed * (1 + z)
sfr_factor = (args.z_dependence_alpha - 1) * (0.5 * ( np.tanh( 0.5 * (z - args.z_middle) ) + 1 ) ) + 1
log_sfr += np.log10(sfr_factor)
log_mstar[log_mstar < -10] = -10
log_Z[log_Z < -10] = -10

### Initialize intensity map
Nx = int(args.side_length / args.angular_resolution)
Nz = len(flist) - 1
total_intensity = np.zeros((Nx, Nx, Nz), dtype=np.float32)
npix = np.array([Nx, Nx, Nz])

### Save intensity maps
output_file=f"{args.output_dir}/map_{args.side_length}sec_fmin{args.fmin}_fmax{args.fmax}_angular{args.angular_resolution}_spectral{args.R}_sigma{args.sigma}_alpha{args.z_dependence_alpha}_th{args.threshold}.h5"

with h5py.File(output_file, "w") as f:
    args_dict = vars(args)
    args_dict = {k: (v if v is not None else "None") for k, v in args_dict.items()}
    for key, value in args_dict.items():
        f.attrs[key] = value

    f.create_dataset("frequency", data=flist, compression="gzip")

    for line_name in list(line_dict.keys()):
        freq_obs = line_dict[line_name][0] / ( 1. + z_obs )
        iz = np.searchsorted(flist, freq_obs, side="right") - 1

        indices = np.array([ix, iy, iz]).T # (num_galaxies, 3)

        valid_mask = np.all((indices >= 0) & (indices < npix), axis=1)
        
        intensity_line = np.zeros((Nx, Nx, Nz), dtype=np.float32)

        if np.sum(valid_mask) > 0:
            print(f"observed frequency of {line_name} is {freq_obs[0]/1e9} - {freq_obs[-1]/1e9} [Hz]")
            print("# Found {} valid galaxies for {}".format(np.sum(valid_mask), line_name))
            print(f"observed redshift range of {line_name} is {np.min(z_obs[valid_mask])} - {np.max(z_obs[valid_mask])}")
            indices_valid = indices[valid_mask]
            log_lumi_dis_valid = log_lumi_dis[valid_mask]
            log_sfr_valid = log_sfr[valid_mask]
            log_mstar_valid = log_mstar[valid_mask]
            log_Z_valid = log_Z[valid_mask]
            z_valid = z[valid_mask]

            log_lumi_valid = calc_line_luminosity(args, z_valid, log_sfr_valid, log_mstar_valid, log_Z_valid, line_name)
            flux_valid = 10 ** ( log_lumi_valid - 2 * log_lumi_dis_valid ) / ( 4. * np.pi ) # [erg/s/cm2]

            if args.threshold < 0:
                intensity_valid = flux_valid / dflist[indices_valid[:, 2]] / Jy / ( args.angular_resolution * arcsec )**2 # [Jy/sr]
            else:
                # 1 if flux_valid > args.threshold else 0
                intensity_valid = (flux_valid > args.threshold).astype(np.float32) 
                
            np.add.at(intensity_line, (indices_valid[:, 0], indices_valid[:, 1], indices_valid[:, 2]), intensity_valid)

        total_intensity += intensity_line

        f.create_dataset("intensity_{}".format(line_name), data=intensity_line, compression="gzip")

    f.create_dataset("total_intensity", data=total_intensity, compression="gzip")

print("Intensity map saved to {}".format(output_file))