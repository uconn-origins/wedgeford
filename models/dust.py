from pylab import *
#import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata

from . units import *


"""
Functions related to modifying dust optical properties

"""
def read_wavelength(fname=None):
    """Reads the wavelength grid

    Parameters
    ----------
    fname : str, optional
                    File name from which the spatial grid should be read. If omitted 'wavelength_micron.inp' will be used.
    """
    if fname is None:
        fname = 'wavelength_micron.inp'
        
    if os.path.exists(fname) != True:
        print('no wavelength file found at path:', fname, 'reading from template')
        templates_dir = model.models_dir+'templates/'
        fname = templates_dir +'wavelength_micron.inp'
    
    print('Reading ' + fname)
    
    data = np.fromfile(fname, count=-1, sep=" ", dtype=np.float64)
    wav = data[1:]
    freq = c / wav * 1e4
    return wav, freq

def run_optool(model, fluid=1,na='70'):
    """ runs the command line options for optool from model inputs as a python script
    using Draine+03 grain compositions 
    
    Parameters:
    -------------
    model: model_ class object 
    
    fluid: int > 0, to generate opacities by dust population number
    
    """
    lam_min = xray_min 
    lam_max = 1.0e4 #longest wavelength
    nlam = 500
    amin = model.dust['amin'][fluid-1] #microns
    amax = model.dust['amax'][fluid-1] #microns
    apow = model.dust['apow'][fluid-1]
    odir = model.outdir
    model.dust['rho_si'] = 3.1518 #combination graphite + astrosilicate grains Draine+03
    cmd = '~/bin/optool -mie -c astrosil 0.87 -c gra 0.13 -p 0.25 -o {} -radmc dust-{} -l {} {} {} -a {} {} {} {}'.format(odir,fluid, lam_min, lam_max, nlam, amin, amax, apow,na)
    os.system(cmd)

def model_kappa_xray(model,wav = None, fluid=1):
    """ calculates gas and dust x-ray absorption opacities using Bethell+Bergin 2011
    if fluid = 1, (i.e. dust follows gas, we add the gas opacities scaled by dust to gas ratio)
    
    Parameters:
    -------------
    model: model_ class object 
    
    wav: wavelength in microns to calculate kappa over, if None, reads in local file
    
    fluid: int > 0, to generate opacities by dust population number
    
    Returns: wav, kappa_tot, the total calculated absorption opacities 
    """
    def calc_fb(model,fluid=1):
        amax = model.dust['amax'][fluid-1]
        fb_file = '{}templates/fb.npy'.format(model.models_dir)
        Fb = np.load(fb_file)
        amax_array = np.linspace(0,3,100)
        index = np.argmin(np.abs(10**amax_array - amax))
        return Fb[index,:]
    
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
    
    kappa_file = '{}templates/kappa_xray.txt'.format(model.models_dir)
    loge,logk_g,logk_d = np.loadtxt(kappa_file)
    Ekev = 10**loge
    lam = 1e4*(h*c)/(Ekev*keV)
    d2g = model.env['d2g']
    if fluid < 2:
        kappa_gas = 10**(logk_g)*(1./(d2g)) #conversion to per grams of dust
        k_g = np.interp(wav,lam[::-1],kappa_gas[::-1],right=0.0,left = 0.0)
    else:
        k_g = np.zeros_like(wav)
    kappa_dust = 10**(logk_d)*(1./(d2g)) #conversion to per grams of dust
    fb = calc_fb(model,fluid=fluid)
    k_d = np.interp(wav,lam[::-1],(kappa_dust*fb)[::-1],right=0.0,left=0.0)
    k_d[wav> 0.001] = 0.0 #only for where optool sucks
    
    kappa_tot = k_g + k_d
    
    return wav, kappa_tot

def read_kappa(model,fluid=1,filename=''):
    """
    counts header length then reads in the columns from dustkappa..inp file, assumes a 4 column file
    
    Parameters:
    -----------
    model: model_ class object
    
    fluid: int, which dust population you want
    
    filename: str, optional, if not specified file is of the form "dustkappa_dust-X" where X is the fluid number
    
    Returns: wav, kabs, kscat, g for selected file
    """
    from scipy import interpolate
    if filename == '':
        kappa_file = '{}dustkappa_dust-{}.inp'.format(model.outdir,fluid)
    else:
        kappa_file = '{}dustkappa_{}.inp'.format(model.outdir,filename)
    
    if os.path.exists(kappa_file):
        with open(kappa_file, 'r') as f:
            numrows = 0
            for line in f:
                if line.startswith('#'):
                    numrows +=1
                else:
                    pass
        numrows += 2
        wav,kabs,kscat,g = np.loadtxt(kappa_file,unpack=True,skiprows=numrows)
        return wav, kabs, kscat, g
    else:
        print('Error: could not find dust opacity file:', kappa_file)
        return None

