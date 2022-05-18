from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from . units import *

###########################################################################################
# models that output the flux density for different kinds of radiation sources
# takes the model class as input
# 
############################################################################################

def UV_powerlaw(model,puv=1.1,LUV=None,flya=0.5,nuv=10): ### power-law, not currently used
    grid = rpy.reggrid.radmc3dGrid()
    grid.readWavelengthGrid()
    rad_data = rpy.radsources.radmc3dRadSources(ppar=model.radmcpar,grid=grid)
    rad_data.getStarSpectrum(grid=opts.data.grid,ppar=model.radmcpar)
    star_flux = rad_data.fnustar
    lam = grid.wav
    star_flux = np.reshape(star_flux,np.shape(star_flux)[0])
    lya_lam = 1216*1e-4 #microns from angstroms
    uv_lam = np.append([uv_min],np.linspace(lya_lam,uv_max,nuv-1))
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

def Xray_powerlaw(model,TX=10e6,nx=20,fX=0.01): ### power-law X-ray emission, not currently used
    grid = rpy.reggrid.radmc3dGrid()
    grid.readWavelengthGrid()
    rad_data = rpy.radsources.radmc3dRadSources(ppar=model.radmcpar,grid=grid)
    rad_data.getStarSpectrum(grid=grid,ppar=model.radmcpar)
    star_flux = rad_data.fnustar
    lam = grid.wav
    star_flux = np.reshape(star_flux,np.shape(star_flux)[0])
    x_lam = np.logspace(np.log10(xray_min),np.log10(xray_max),nx) #nm to microns
    nu_x = c/(x_lam*1e-4)# Hz
    nu = c/(lam*1e-4) #Hz  
    I_nu = (nu_x[-1]/nu_x)*np.exp(-(h*nu_x)/(kb*TX))
    Fstar = np.trapz(star_flux,x=nu) #total flux
    Ftot = np.trapz(I_nu,x=nu_x) #flux in Xray
    norm = (Fstar*fX)/(Ftot) # Flux in the X-ray is fX of total flux
    return x_lam, I_nu*norm

##### model for X-ray emission from accretion shock onto the star
def acc_onto_star(model,lam_Xray):
    accrate = model.star['accrate']
    f = model.star['f']
    Area = f*4.*pi*(model.star['Rs']*Rsun)**2
    Area_norm = 4.*pi*(pc**2)
    vs = 3.1e7*np.sqrt(model.star['Ms']/(model.star['Rs']))*np.sqrt(4/5.)
    rhos = (accrate*Msun/yr)/(vs*Area)
    Ts = 8.6e5*(model.star['Ms']/model.star['Rs'])
    vps = vs*0.25
    rhops = rhos*4.0
    Pram = 0.5*rhos*(vs**2)
    C = (93*np.sqrt(3)-40*pi)*np.sqrt(mh**5/kb)*np.sqrt(2**5)/(2*(np.sqrt(2+0.25+1)))/1.64e-18
    ds_cool = C*(vs**4)/(rhos)
    rho_cool = Pram*(mu*mh)/(0.8*kb*model.star['Ts'])
    v_cool = rhos*vs/rho_cool
    T_cool = Ts*np.sqrt(v_cool/vps)
    #lam_Xray = np.logspace(np.log10(xray_min),np.log10(xray_max),nx) #microns
    nu_Xray = c/(lam_Xray*1e-4) 
    def jff(rho,T): #bremmstrahlung emission
        gff = np.abs(np.sqrt(3)*np.log(kb*T/(h*nu_Xray))/pi)
        n = rho/(mu*mh)
        return 5.44e-39*gff*(n**2)*np.sqrt(1./T)*exp(-h*nu_Xray/(kb*T))
    F_1 = (jff(rhops,Ts)*Area*ds_cool*0.5)/Area_norm
    F_2 = (jff(rho_cool,T_cool)*Area*ds_cool*0.5)/Area_norm
    
    xmodel = {}
    xmodel['shock'] = {'T': Ts,'rho':rhops,'ds':ds_cool/2.,'F':F_1,'lam':lam_Xray}
    xmodel['cool'] = {'T':T_cool,'rho':rho_cool,'ds':ds_cool/2,'F':F_2,'lam':lam_Xray}
    return xmodel
    
def Xray_accretion(model,fX=None,nx=20):
    lam_Xray = np.logspace(np.log10(xray_min),np.log10(xray_max),nx)
    xmodel = acc_onto_star(model,lam_Xray)
    if fX != None:
        lam, fnu = Xray_powerlaw(xmodel,TX=xmodel['shock']['T'],nx=nx,fX=fX) #will scale hottest model component based on desired Xray flux fraction
    else:
        lam = xmodel['shock']['lam']
        fnu = xmodel['shock']['F'] + xmodel['cool']['F'] #otherwise will add up cooling and shock layer
    return lam, fnu

###### model for inserting a Lyman alpha line:

 
def Lya_line(model, fLya=1e-2):
    lam_nm, flya = np.loadtxt(model.models_dir + 'templates/spec_lya.txt',skiprows=4, unpack=True)
    lam_mu = lam_nm*1e-3 #microns
    nu = c/(lam_mu*1e-4)
    fnu = flya*(h*nu)* (c*1e7) / (nu**2)
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    grid = rpy.reggrid.radmc3dGrid()
    grid.readWavelengthGrid()
    rad_data = rpy.radsources.radmc3dRadSources(ppar=model.radmcpar,grid=grid)
    rad_data.nstar = 1
    rad_data.getStarSpectrum(grid=grid,ppar=model.radmcpar)
    fnu_s = rad_data.fnustar.squeeze()
    rad_data.getSpotSpectrum(grid=grid,ppar=model.radmcpar)
    fnu_uv = rad_data.fnuspot.squeeze()
    lam_s = rad_data.grid.wav
    nu_s = c/(lam_s*1e-4)
    
    Fstar = np.trapz(fnu_s,x=nu_s) #total stellar flux
    Ftot = np.trapz(fnu,x=nu) #total lyman alpha flux
    norm = (Fstar*fLya)/(Ftot) #scaling so that the total in the lya line is fLya of the stellar flux
    fpeak_lya = np.amax(norm*fnu) #peak lyman alpha flux to insert into the RT calculation
    #fnu_lya = interpolate.interp1d(x= nu, y = fnu*norm, fill_value=np.amin(fnu*norm),bounds_error=False) #function to put in a line profile
    return fpeak_lya

