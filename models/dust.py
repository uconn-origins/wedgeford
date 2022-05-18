from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata

from . units import *




def run_optool(model, fluid=1,na=''):
    lam_min = xray_min #10 keV in microns
    lam_max = 1.0e5 #longest wavelength
    nlam = 500
    amin = model.dust['amin'][fluid-1] #microns
    amax = model.dust['amax'][fluid-1] #microns
    apow = model.dust['apow'][fluid-1]
    odir = model.outdir
    model.dust['rho_si'] = 3.1518 #combination graphite + astrosilicate grains Draine+03
    cmd = '~/bin/optool -mie -c astrosil 0.87 -c gra 0.13 -p 0.25 -o {} -radmc dust-{} -l {} {} {} -a {} {} {} {}'.format(odir,fluid, lam_min, lam_max, nlam, amin, amax, apow,na)
    os.system(cmd)

def kappa_xray(model,fluid=1):
    def calc_fb(model,fluid=1):
        amax = model.dust['amax'][fluid-1]
        fb_file = '{}templates/fb.npy'.format(model.models_dir)
        Fb = np.load(fb_file)
        amax_array = np.linspace(0,3,100)
        index = np.argmin(np.abs(10**amax_array - amax))
        return Fb[index,:]
    kappa_file = '{}templates/kappa_xray.txt'.format(model.models_dir)
    loge,logk_g,logk_d = np.loadtxt(kappa_file)
    Ekev = 10**loge
    lam = 1e4*(h*c)/(Ekev*keV)
    d2g = model.disk['Mfrac'][fluid-1]
    kappa_gas = 10**(logk_g)*(1./(d2g)) #conversion to per grams of dust
    kappa_dust = 10**(logk_d)*(1./(d2g)) #converstion to per grams of dust
    fb = calc_fb(model,fluid=fluid)
    return lam, kappa_gas, fb*kappa_dust

def amend_kappa(model,fluid=1,filename=''):
    from scipy import interpolate
    if filename == '':
        kappa_file = '{}dustkappa_dust-{}.inp'.format(model.outdir,fluid)
    else:
        kappa_file = '{}dustkappa_{}.inp'.format(model.outdir,filename)
    lam,kabs,kscat,g= np.loadtxt(kappa_file,skiprows=26,unpack=True)
    lam_pe,kappa_g,kappa_d = kappa_xray(model,fluid=fluid)
    lam_thr = 0.001 #cutoff for optool calculation
    lam_max = np.amax(lam_pe)
    lam_1 = lam[(lam <= lam_max)]
    lam_2 = lam[(lam <= lam_thr)]
    f_1 = interpolate.interp1d(np.log10(lam_pe),np.log10(kappa_g),fill_value='extrapolate')
    f_2 = interpolate.interp1d(np.log10(lam_pe),np.log10(kappa_d),fill_value='extrapolate')
    kabs_gas = 10**f_1(np.log10(lam_1))
    kabs_dust = 10**f_2(np.log10(lam_2))
    kabs[(lam <= lam_thr)] = kabs_dust
    if fluid == 1:
        kabs[(lam <= lam_max)] += kabs_gas
    kabs[np.isnan(kabs)] = kabs[~np.isnan(kabs)][0]
    return lam,kabs,kscat,g
