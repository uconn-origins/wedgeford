from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from models.units import *

###########################################################################################
# models that output the flux density for different kinds of high energy radiation sources
# needs the out class object to be passed through to read in the star information
# 
############################################################################################

def UV_powerlaw(outputs,puv=1.1,LUV=None,flya=0.5,nuv=10):
    opts = outputs
    if 'star' in opts.rad.keys():
        lam = opts.rad['star']['lam']
        star_flux = opts.rad['star']['fnu']
    else:
        opts.rad_source(field='star')
        lam = opts.rad['star']['lam']
        star_flux = opts.rad['star']['fnu']
    star_flux = np.reshape(star_flux,np.shape(star_flux)[0])
    lya_lam = 1216*1e-4 #microns from angstroms
    uv_lam = np.append([0.1],np.linspace(lya_lam,0.205,nuv-1))
    I_nu = (uv_lam/uv_lam[-1])**(puv+2)
    nu_uv = c/(uv_lam*1e-4)# Hz
    nu = c/(lam*1e-4) #Hz
    if LUV != None: # scale by fraction of stellar luminosity from UV band
        Fstar = np.trapz(star_flux,x=nu) #total flux
        Ftot = np.trapz(pi*I_nu,x=nu_uv) #flux in UV band
        norm = ((LUV*Fstar)/Ftot)
    else: #scale to match the stellar luminosity at the edge of the UV band
        norm = np.amax(star_flux[lam<=np.amax(uv_lam)])
    uv_flux = I_nu*norm
    Ftot = np.trapz(uv_flux[::-1],x=nu_uv[::-1])
    dnu_lya = -np.gradient(nu_uv)[1]
    uv_flux[1] = Ftot/((1./flya) - 1.)/dnu_lya #set the lyman alpha flux to be flya of the total
    return uv_lam, uv_flux

##### to-do better Lyman alpha with the line emission from 
##### https://home.strw.leidenuniv.nl/~ewine/photo/radiation_fields.html

def Xray_powerlaw(outputs,TX=10e6,nx=20,fX=0.01):
    opts = outputs
    opts.rad_source(field='star')
    lam = opts.rad['star']['lam']
    star_flux = opts.rad['star']['fnu']
    star_flux = np.reshape(star_flux,np.shape(star_flux)[0])
    x_lam = np.logspace(-1,1,nx)*1e-3 #nm to microns
    nu_x = c/(x_lam*1e-4)# Hz
    nu = c/(lam*1e-4) #Hz  
    I_nu = (nu_x[-1]/nu_x)*np.exp(-(h*nu_x)/(kb*TX))
    Fstar = np.trapz(star_flux,x=nu) #total flux
    Ftot = np.trapz(I_nu,x=nu_x) #flux in Xray
    norm = (Fstar*fX)/(Ftot) # Flux in the X-ray is fX of total flux
    return x_lam, I_nu*norm


###### TO-D0: self-consistent UV + X-Ray emission using the Calvet + Gullbring Accretion Shock models
###### The UV/optical component of this already in radmc3d's "accreting spot spectrum" method, just need to tie it into the 
###### X-ray scaling of the bremsstrahlung ---> UV and X-ray gets done at the same time

def write_radiation_field(outputs,lam,fnu): #writes the wavelength and spectrum info for the new fields
    opts = outputs
    opts.m.write_wavelength(fname='wavelength_micron_he.inp',lam=lam)
    with open(opts.m.outdir+'stars.inp',"r") as f:
        f.readline()
        f.readline()
        line = f.readline()
    f.close()
    with open(opts.m.outdir+'stars_he.inp','w') as f:
        f.write('2 \n')
        f.write('1 \t')
        f.write(str(len(lam))+ '\n')
        f.write(line)
        f.write('\n')
        lam.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')
        f.write('\n')
        fnu.tofile(f, sep='\n', format="%13.6e")
    f.close()
    

