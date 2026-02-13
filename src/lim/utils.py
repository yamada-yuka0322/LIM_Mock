import numpy as np
from numpy.fft import fftn, fftfreq

import sys

from scipy.ndimage import gaussian_filter

import astropy.cosmology
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u

import os
from functools import lru_cache

from Corrfunc.mocks import DDtheta_mocks

import h5py
import pandas as pd

arcmin = 2.908882086656e-4 # [rad] ... PI / 180 / 60 //
arcsec = 4.848136811094e-6 # [rad] ... arcmin / 60 //
GHz = 1e9
Jy = 1.0e-23

# 可能ならスレッド数を明示（サーバで効く）
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")


# --- rfft2 の互換ラッパ（SciPy→NumPy workers→NumPy互換） ---
try:
    import scipy.fft as sfft  # あれば workers 対応
    def _rfft2(x, axes=(-2, -1), workers=None):
        return sfft.rfft2(x, axes=axes, workers=workers)
except Exception:
    def _rfft2(x, axes=(-2, -1), workers=None):
        try:
            # pocketfft なら workers が使える
            return np.fft.rfft2(x, axes=axes, workers=workers)
        except TypeError:
            # mkl-fft 等：workers 非対応
            return np.fft.rfft2(x, axes=axes)


def _as_f32_contig(a):
    return np.ascontiguousarray(a, dtype=np.float32)

@lru_cache(maxsize=64)
def _cached_bins(Ny, Nx, Lx, Ly, dlogk):
    """
    rfft2 の形状 (Ny, Nx//2+1) に対応する |k| のビン情報をキャッシュ
    戻り値:
      bin_index: (Ny, Nx//2+1) 各モードのビン番号（[-1: 無効]）
      k_center: (Nbins,)       ビン中心
      counts:   (Nbins,)       各ビンのモード数（0含む）
    """
    dx = Lx / Nx
    dy = Ly / Ny
    kx = 2*np.pi * np.fft.rfftfreq(Nx, d=dx)
    ky = 2*np.pi * np.fft.fftfreq(Ny, d=dy)
    KX, KY = np.meshgrid(kx, ky, indexing='xy')
    kabs = np.sqrt(KX**2 + KY**2)

    # 物理的に使う k 範囲
    k_min = 2.0*np.pi / max(Lx, Ly)
    k_nyq_x = np.pi / dx
    k_nyq_y = np.pi / dy
    k_max = 0.8 * min(k_nyq_x, k_nyq_y)

    # 対数ビン
    lo, hi = np.log10(k_min), np.log10(k_max)
    # bin edges（左閉右開にしたいので、edges を一つ多く）
    edges = 10**np.arange(lo, hi + dlogk, dlogk)
    # 中心
    k_center = np.sqrt(edges[1:] * edges[:-1])

    # digitize でビン番号（0..Nbins-1）。範囲外を -1 に
    with np.errstate(invalid='ignore'):
        idx = np.digitize(kabs.ravel(), edges, right=False) - 1
    Nbins = edges.size - 1
    bad = (idx < 0) | (idx >= Nbins) | ~np.isfinite(kabs.ravel())
    idx[bad] = -1
    idx = idx.reshape(kabs.shape)

    # counts だけ先に（0埋め）
    counts = np.bincount(idx[idx >= 0], minlength=Nbins).astype(np.int64)

    return idx.astype(np.int32), k_center.astype(np.float32), counts

def my_fft2d_real(X, Lxy=None, b=1.0, workers=None):
    """
    2D FFT（実入力）を物理的に素直な正規化で返す。

    Parameters
    ----------
    X : 2D array
        実空間マップ（例: intensity [Jy/sr]）
    Lxy : float or (Lx, Ly)
        領域サイズ（X, Y 方向）。ここでは [deg] を想定。
    b : float
        ここではスケーリングに使わない（保持だけ）。
    workers : int or None
        numpy FFT の並列数。

    Returns
    -------
    Ft : 2D complex array (Ny, Nx//2+1)
        連続極限に対応する FT:
        Ft(k) ≈ ∫ I(θ) e^{-ik·θ} d^2θ
    kx, ky : 1D arrays
        各軸の波数 [1/deg]（2π/Lx スケール）
    V : float
        全領域の面積（Lx * Ly） [deg^2]
    """
    X = np.ascontiguousarray(X, dtype=np.float32) #Jy/sr
    Ny, Nx = X.shape

    if Lxy is None:
        Lx = float(Nx)
        Ly = float(Ny)
    else:
        if np.isscalar(Lxy):
            Lx = Ly = float(Lxy)
        else:
            Lx, Ly = float(Lxy[0])*(np.pi/180.0), float(Lxy[1])*(np.pi/180.0) #radian

    # ピクセルスケールとセル面積
    dx, dy = Lx / Nx, Ly / Ny   # size of 1 pixel [rad]
    Vcell = dx * dy             # size of 1 pixel [sr]
    V = Lx * Ly                 # total area [sr]

    # k の定義（単位は 1/deg）
    kx = 2 * np.pi * np.fft.rfftfreq(Nx, d=dx/(np.pi/180.0))  # [1/deg]
    ky = 2 * np.pi * np.fft.fftfreq(Ny,  d=dy/(np.pi/180.0))  # [1/deg]

    # FFT: ∑ I_ij e^{-ik·θ_ij} Δθ_x Δθ_y を近似
    Ft = _rfft2(X, axes=(-2, -1), workers=workers) * Vcell * (1/np.abs(b))**0.5 # [Jy]

    # b は今のところスケーリングに使わない（b≠1で何かしたいならここに入れる）
    return Ft, kx, ky, V


def calc_power_2d_fast(params, intensity, intensity1=None, dlogk=0.15, return_k_edges=False, workers=None):
    """
    高速版：rfft2 + 一括ビニング。
    workers: numpy FFT の並列数（例: os.cpu_count()）
    """
    # 領域サイズ（あなたの定義に合わせる）
    Lx = Ly = params.radius  # [deg] など

    # FFT（左側）
    Ft1, kx, ky, V = my_fft2d_real(intensity, Lxy=(Lx, Ly), b=1.0, workers=workers)

    if intensity1 is None:
        # Auto power：|F1|^2 / V
        # rfft の対称性に伴う2倍則は、実空間の定義/正規化に依る。
        # 相対比較が主ならそのままでOK。厳密な規格化は分析系に合わせて調整。
        Pfield = (Ft1 * np.conj(Ft1)).real / V
    else:
        # Cross power: Re(F1 * F2*) / V
        Ft2, _, _, _ = my_fft2d_real(intensity1, Lxy=(Lx, Ly), b=1.0, workers=workers)
        Pfield = (Ft1 * np.conj(Ft2)).real / V #[Jy^2/sr] 

    Ny, Nx_half = Pfield.shape
    bin_index, k_center, counts = _cached_bins(Ny, (Nx_half-1)*2, Lx, Ly, dlogk)

    # 無効を除去して 1D に
    valid = bin_index >= 0
    idx = bin_index[valid].ravel()
    vals = Pfield[valid].ravel().astype(np.float64, copy=False) #[Jy^2/sr] 

    Nbins = k_center.size
    # 合計と二乗和を一気に
    sum_w  = np.bincount(idx, weights=vals, minlength=Nbins) #[Jy^2/sr] 
    sum_w2 = np.bincount(idx, weights=vals*vals, minlength=Nbins) #[(Jy^2/sr)**2] 
    nmode  = np.bincount(idx, minlength=Nbins).astype(np.int64)

    # 平均
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = sum_w / nmode
        # 不偏分散 ddof=1
        var  = (sum_w2 - (sum_w*sum_w)/nmode) / np.maximum(nmode - 1, 1)
        stderr = np.sqrt(var) / np.sqrt(nmode)

    # ビンに1点も無ければ NaN
    mean[nmode == 0] = np.nan
    stderr[nmode <= 1] = np.nan

    if return_k_edges:
        # エッジはキャッシュ内で閉じているので再構成もできるが、必要なら _cached_bins を拡張
        pass

    return k_center, mean, stderr

def c(X, Lxy=None, b=1.0, workers=None):
    X = np.ascontiguousarray(X, dtype=np.float32)
    Ny, Nx = X.shape
    if Lxy is None:
        Lx = float(Nx); Ly = float(Ny)
    else:
        if np.isscalar(Lxy):
            Lx = Ly = float(Lxy)
        else:
            Lx, Ly = float(Lxy[0]), float(Lxy[1])

    dx, dy = Lx / Nx, Ly / Ny
    Vcell = dx * dy
    V = Lx * Ly

    kx = 2*np.pi * np.fft.rfftfreq(Nx, d=dx)
    ky = 2*np.pi * np.fft.fftfreq(Ny,  d=dy)

    # ← ここだけラッパで呼ぶ
    Ft = Vcell * _rfft2(X, axes=(-2, -1), workers=workers) * (1/np.abs(b))**0.5
    return Ft, kx, ky, V



def prepare_data(data, params, half=True, noise=0, psf=0.0):
    pixel = params.resolution
    FWHM = 2*np.sqrt(2*np.log(2))
    if (psf>0.0):
        Sigma = psf/FWHM/ pixel
        data = gaussian_filter(data, sigma=(Sigma, Sigma, 0))
    
    if noise>0:
        if half:
            f_nu = 0.2*10**(-0.4*(noise+48.6))*np.sqrt(2)
        else:
            f_nu = 0.2*10**(-0.4*(noise+48.6))
        sigma_noise = f_nu / Jy/ (pixel * arcsec)**2 # [Jy/sr]
        gauss_noise = np.random.normal(0, sigma_noise, data.shape)
        data += gauss_noise
    return data

def arcsec_to_cMpc(l_arcsec, z):
    l_rad = l_arcsec * u.arcsec / u.radian
    l_cMpc = ( cosmo.comoving_transverse_distance(z) * l_rad ).to(u.Mpc)
    return l_cMpc.value 

def freq_to_comdis(nu_obs, nu_rest):
    z = nu_rest / nu_obs - 1
    if z < 0:
        print("Error: z < 0")
        sys.exit(1)
    return cosmo.comoving_distance(z).to(u.Mpc).value

def my_fft(X, L=None, b=1.): # b = 2 * np.pi for inverse FT
    """
    Actually whatever dimension is ok
    input:
        X: (Nx, Ny, Nz) data 
        L: (3,) or float, box size
        b: normalization factor
    output:
        ft: (Nx, Ny, Nz) Fourier transform
        freq: (Nx, Ny, Nz) frequency
    """

    dim = len(X.shape)
    N = np.array(X.shape)

    if L is None:
        L = N
    elif np.isscalar(L):
        L = L * np.ones(dim)

    dx = np.array([float(l)/float(n) for l, n in zip(L, N)])
    Vx = np.prod(dx) # volume of a cell

    freq = [fftfreq(n, d=d)*2*np.pi for n, d in zip(N, dx)]
    ft = Vx * fftn(X) * np.sqrt(1/np.abs(b)) ** dim
    
    return ft, freq

def calc_power_2d(params, intensity, intensity1=None, dlogk=0.15, return_k_edges=False):
    """
    input:
        intensity:  (Nx, Ny)      2D map
        intensity1: (Nx, Ny) or None
                     None なら auto power, 与えられたら cross power を計算
        dlogk:      対数kのビン幅
    output:
        k_center: (Nk,)  ビン中心の |k|（単位は my_fft の定義に依存）
        Pk:       (Nk,)  角平均パワー（auto: |F1|^2/V, cross: Re(F1*F2*)/V）
        Pk_err:   (Nk,)  各ビンの標準誤差（std/√Nmode）
    """
    side_length = params.radius  # [deg] （あなたの my_fft がこの単位に対応している前提）

    # Fourier transform
    L = np.array([side_length, side_length])
    V = np.prod(L)

    ft1, freq = my_fft(intensity, L=L)        # freq = (kx, ky)
    if intensity1 is None:
        # Auto power
        power_field = (ft1 * np.conj(ft1)).real / V
    else:
        # Cross power
        if intensity1.shape != intensity.shape:
            raise ValueError(f"intensity1 shape {intensity1.shape} != intensity shape {intensity.shape}")
        ft2, _ = my_fft(intensity1, L=L)
        power_field = (ft1 * np.conj(ft2)).real / V  # 実部のみ採用

    kx, ky = np.meshgrid(freq[0], freq[1], indexing='ij')
    k_abs = np.sqrt(kx**2 + ky**2)

    # 対数ビン（k=0 は log 取れないので自然に除外される）
    logk = np.log10(k_abs, where=(k_abs > 0), out=np.full_like(k_abs, -np.inf))

    k_min = 1.0 / side_length                          # [deg^-1]
    dx = params.resolution / 3600.0
    k_nyq = 1.0 / (2.0 * dx)                 # [deg^-1]
    k_max = 0.8 * k_nyq                      # やや余裕を持たせる
    lo, hi = np.log10(k_min), np.log10(k_max)
    logk_bins = np.arange(lo, hi + dlogk + 1e-12, dlogk)

    Nk = len(logk_bins) - 1
    power1d = np.full(Nk, np.nan, dtype=float)
    power1d_err = np.full(Nk, np.nan, dtype=float)

    for i in range(Nk):
        mask = (logk >= logk_bins[i]) & (logk < logk_bins[i+1])
        if not np.any(mask):
            continue
        vals = power_field[mask]
        power1d[i]     = np.mean(vals)
        # 標準誤差（ビン内 scatter/√Nmode）
        if vals.size > 1:
            power1d_err[i] = np.std(vals, ddof=1) / np.sqrt(vals.size)
        else:
            power1d_err[i] = np.nan
            
    counts = []
    for i in range(len(edges)-1):
        m = (logk >= edges[i]) & (logk < edges[i+1])
        counts.append(np.count_nonzero(m))
    print("nmode per bin:", counts[:10], "…")

    k_center = 10**(0.5 * (logk_bins[1:] + logk_bins[:-1]))
    if return_k_edges:
        return k_center, power1d, power1d_err, 10**logk_bins
    return k_center, power1d, power1d_err


def wtheta_landy_szalay(ra, dec, ra_rand, dec_rand, bins_deg, nthreads=8):
    # Corrfunc expects degrees for RA/DEC and theta bins in degrees for DDtheta_mocks
    # autocorr=1 for auto-correlation; pass RA2/DEC2 only for cross-corr.
    DD = DDtheta_mocks(autocorr=1, nthreads=nthreads, binfile=bins_deg,
                       RA1=ra, DEC1=dec, RA2=None, DEC2=None)
    RR = DDtheta_mocks(autocorr=1, nthreads=nthreads, binfile=bins_deg,
                       RA1=ra_rand, DEC1=dec_rand, RA2=None, DEC2=None)
    DR = DDtheta_mocks(autocorr=0, nthreads=nthreads, binfile=bins_deg,
                       RA1=ra, DEC1=dec, RA2=ra_rand, DEC2=dec_rand)

    # Extract pair counts per bin
    dd = DD['npairs'].astype(float)
    rr = RR['npairs'].astype(float)
    dr = DR['npairs'].astype(float)

    # Normalize counts (important when ND != NR)
    ND = len(ra)
    NR = len(ra_rand)

    # Number of possible pairs
    #norm_DD = ND*(ND-1)/2.0
    #norm_RR = NR*(NR-1)/2.0
    #norm_DR = ND*NR

    #DDn = dd / norm_DD
    #RRn = rr / norm_RR
    #DRn = dr / norm_DR

    #w = (DDn - 2.0*DRn + RRn) / RRn
    w = 1.0 + (NR/ND)**2*(dd/rr) - 2.0*(NR/ND)*(dr/rr)

    # handy x-axis: bin midpoints
    theta_lo = DD['thetamin']
    theta_hi = DD['thetamax']
    theta_mid = (theta_lo * theta_hi)**0.5

    return theta_mid, w

def build_box_stars(ra, dec, spectrum, params, center_ra, center_dec, flist):
    Nx = int(3600 * params.radius / params.resolution)
    Nz = len(flist)
    star_intensity = np.zeros((Nx, Nx, Nz), dtype=np.float32)

    ra_min = center_ra - params.radius/2.0
    dec_min = center_dec - params.radius/2.0
    
    iy = ((ra - ra_min) / params.resolution * 3600.0).astype(np.int32)   # x in your code
    ix = ((dec - dec_min) / params.resolution * 3600.0).astype(np.int32) # y in your code
    
    inside = (ix >= 0) & (ix < Nx) & (iy >= 0) & (iy < Nx)

    ix = ix[inside]
    iy = iy[inside]
    spectrum = spectrum[inside]
    print(f'number of star is {np.sum(inside)}')
    
    pid = (ix.astype(np.int64) * Nx + iy.astype(np.int64))
    Npix = Nx * Nx
    
    for k in range(Nz):
        tmp = np.bincount(pid, weights=spectrum[:, k], minlength=Npix)
        star_intensity[:,:,k] = tmp.reshape(Nx, Nx)
        
    star_intensity /= (params.resolution * arcsec)**2  # to [erg/s/cm^2/Hz/sr]
        
    return star_intensity

def get_EL(path, name):
    meta = pd.read_parquet(path)
    F = meta[name+'_F'].to_numpy()
    iz = meta[name+'_iz'].to_numpy()
    return F, iz

def build_box_galaxies(ra, dec, spec_path, meta_path, params, flist, type='total'):
    Nx = int(2.0 * 3600 * params.radius / params.resolution)
    Nz = len(flist)
    galaxy_intensity = np.zeros((Nx, Nx, Nz), dtype=np.float32)

    ra_min = np.min(ra)
    dec_min = np.min(dec)
    
    iy = ((ra - ra_min) / params.resolution * 3600.0).astype(np.int32)   # x in your code
    ix = ((dec - dec_min) / params.resolution * 3600.0).astype(np.int32) # y in your code
    
    with h5py.File(spec_path, "r") as h5:
        if type == 'total':
            SED = h5["fnu_total"][:]
        elif type == 'line':
            SED = h5["fnu_lines"][:]
        elif type == 'cont':
            SED = h5["fnu_cont"][:]
        else:
            flux, iz = get_EL(meta_path, type)
            SED = np.zeros((len(ra), Nz), dtype=np.float32)
            SED[np.arange(len(ra)), iz] = flux
             
    
    pid = (ix.astype(np.int64) * Nx + iy.astype(np.int64))
    Npix = Nx * Nx
    
    for k in range(Nz):
        tmp = np.bincount(pid, weights=SED[:, k], minlength=Npix)
        galaxy_intensity[:,:,k]= tmp.astype(np.float32)
        
    galaxy_intensity /= (params.resolution * arcsec)**2  # to [erg/s/cm^2/Hz/sr]
        
    return galaxy_intensity
    
