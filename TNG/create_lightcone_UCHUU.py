import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
import scipy
from multiprocessing import Pool
from functools import partial
from scipy.interpolate import interp1d
import h5py
import os
import sys
from scipy.special import erf


class Sides(object):
    def __init__(self):
        self.position = None
        self.velocity = None
        self.mass = None
        self.Mstar = None
        self.SFR = None
        self.metalicity = None
        self.redshift = None

    def loadUCHUU(self, path, redshift):
        _redshift = str(redshift).replace('.', 'p')
        filename = f"MiniUchuu_halolist_z{_redshift}.h5"
        self.redshift = redshift
        with h5py.File(os.path.join(path, "MiniUchuu_halolist_z0p09.h5"), "r") as f:
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

        Sparam = {
        'A':0.0353,
        'MA':12.05,
        'beta':0.88,
        'gamma':0.599
        }

        self.Mstar = get_starM(self.mass, 0.6774, scatter=0.2, **Sparam)
        #self.SFR = get_SFR(self)


def get_starM(halo_mass, h, scatter = 0.2, **param):
    halo = halo_mass/h
    MA = np.power(10.0, param['MA'])
    ratio = 2*param['A']*np.power((np.power((halo/MA), -param['beta'])+np.power((halo/MA), param['gamma'])), -1.0)
    stellar = ratio * halo
    stellar_scattered = stellar * np.power(10.0, (np.random.normal(0, scatter, size=len(stellar))))
    return stellar_scattered*h

def get_SFR(sides):
    base_dir = os.path.dirname(__file__)  # この.pyファイルのディレクトリ
    full_path = os.path.join(base_dir, 'params/SIDES.par')
    params = load_params(full_path)
    print('Generate the star-formation properties...')

    print('Draw quenched galaxies...')

    Ngal = len(sides.mass)

    #Draw randomly which galaxies are quenched using the recipe from Bethermin+17
    Mtz = params['Mt0'] + params['alpha1'] * sides.redshift + params['alpha2'] * sides.redshift**2
    sigmaz =  params['sigma0'] +  params['beta1'] * sides.redshift +  params['beta2'] * sides.redshift**2
    qfrac0z = params['qfrac0'] * (1.+sides.redshift)**params['gamma']
    
    Prob_SF = (1.-qfrac0z) * 0.5 * (1. - erf( ( np.log10(sides.Mstar) - Mtz) /sigmaz ) )

    Xuni = np.random.rand(Ngal)

    qflag = Xuni > Prob_SF
    #Generate SFR for non-quenched objects

    print('Generate SFRs...')

    index_SF = np.where(qflag == False)

    m_SF = np.array(np.log10(sides.Mstar[index_SF[0]] * params['Chab2Salp'] / 1.e9 ))
    z_SF = np.array(sides.redshift[index_SF[0]])
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
