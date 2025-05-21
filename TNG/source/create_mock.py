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

import pandas as pd

import h5py

from line_model import calc_line_luminosity, line_dict

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