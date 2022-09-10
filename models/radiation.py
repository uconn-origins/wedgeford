from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from . units import *

"""
Functions for calculations related to radmc3d radiation sources



"""
def read_wavelength(fname=None):
    """Reads the wavelength grid

    Parameters
    ----------
    fname : str, optional
            Filename from which the wavelength grid should be read. If omitted 'wavelength_micron.inp' will be used.
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

def calc_blackbody_spectrum(model,wav = None):
    """
    calculates the blackbody spectrum of star at input temperature and radius
    
    Parameters:
    -----------
    model: model_ object class 
    
    wav: wavelength in microns, if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
            
    Returns: wav, fnu, ndarrays for the frequencies and blackbody spectrum at every input wavelength in ergs/cm^2/s/Hz 
    """
    
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
    Tstar = model.star['Ts']
    Rstar = model.star['Rs']*Rsun
    
    fnu =  2. * h * freq**3. / c**2 / (np.exp(h * freq / (kb * Tstar) )- 1.0) * np.pi * Rstar**2. / pc**2.
    
    return wav, fnu

def model_accrate_infall(model):
    """ calculates a value for the expected disk-wide accretion rate due to infall from Kuznetsova+2022
    Parameters:
    ----------
    model: instance of model_ class
    
    Returns: value for the accretion rate in Msun/yr
    """
    Rc = model.env['Rc']
    Min = model.env['Min']
    Rin = Rc*np.sin(np.radians(model.env['theta_min']))**2
    tauc = Rc**(3/2.)
    sigc = model.sig_profile(model.r)[model.r <= Rc][-1]
    sigin = (Min*Msun)/(4*pi*((Rc*AU)**2-(Rin*AU)**2)) * tauc
    logalpha = min(2*(np.log10(sigin) - np.log10(sigc)) - 0.5, -1.75)
    model.disk['alpha'] = 10**logalpha
    
    Ms = model.star['Ms']
    R0 = model.disk['R0'][0]
    
    vc = np.sqrt(Gconv*Ms*Msun/(R0*AU))
    h0 = model.H(R0)
    rho0 = model.rho_disk_profile(r=np.array([R0,R0]),z=np.array([0,0]))[0]
    Macc = sqrt(18*pi**3)*(10**logalpha)*vc*rho0*(h0*AU)**3/(R0*AU)/(Msun/yr)

    return Macc

def calc_accretion_spectrum_uv(model,wav = None,accrate=None):
    """ calculates the uv continuum spectrum due to accretion onto the stellar surface
    Parameters:
    ------------
    model: model_ object class 
    
    wav: wavelenghth in microns if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
    
    accrate: acccretion rate in Msun/yr, if None, model input is used
            
    Returns: wav, fnu, ndarrays for the frequencies and spectrum at every input wavelength in ergs/cm^2/s/Hz 
    """
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
        
    Mstar = model.star['Ms']*Msun
    Rstar = model.star['Rs']*Rsun
    
    surffrac = model.star['f']
    
    if accrate is None:
        accrate = model.star['accrate']

    tot_acclum = 0.5 * Gconv * Mstar * accrate * (Msun / yr) / Rstar
    spotsurf = 4. * np.pi * Rstar**2 * surffrac
    tspot = (tot_acclum / spotsurf / sigsb)**0.25
    fnu = (np.pi * Rstar**2 * surffrac / (pc**2) ) * (2. * h * freq**3 / c**2) / (np.exp(h * freq / (kb * tspot)) - 1.0)
    
    return wav, fnu

def calc_accretion_spectrum_xray(model,wav=None,accrate=None,fX=None):
    """ calculates the xray spectrum due to shocks from accretion onto the stellar surface
    Parameters:
    ------------
    model: model_ object class 
    
    wav: wavelenghth in microns if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
    
    accrate: accretion rate in Msun/yr, if None, uses model accretion rate
    
    fX: optional, if not None, then used to scale a power-law matching the xray temperature of the shocked component
            
    Returns: wav, fnu, ndarrays for the frequencies and spectrum at every input wavelength in ergs/cm^2/s/Hz 
    """
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
    
    if accrate is None:
        accrate = model.star['accrate']
        
    xmodel = model_acc_onto_star(model,wav=wav,accrate=accrate)
    
    if fX is None:
        fnu = xmodel['shock']['F'] + xmodel['cool']['F'] #otherwise will add up cooling and shock layer
    else:
        wav, fnu = calc_powerlaw_xray(model,wav=wav,TX=xmodel['shock']['T'],fX=fX) #will scale hottest model component based on desired Xray flux fraction
    return wav, fnu


def model_acc_onto_star(model,wav=None,accrate=None):
    """ calculates the xray spectrum due to shocks from accretion onto the stellar surface
    based on Calvet+Gullbring 1998/Brickhouse 2013
    
    two component accretion shock model: shock surface + cooling layer
    flux is bremmstrahlung power-law emission profile
    
    Parameters:
    ------------
    model: model_ object class 
    
    wav: wavelength in microns if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
            
    accrate: acccretion rate in Msun/yr, if None, model input is used
            
    Returns: xmodel, dict, parameters of two component shock model
            T, rho, ds (width), flux, wavelength
    """
    if accrate is None:
        accrate = model.star['accrate']
        
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
        
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
    
    def jff(rho,T,freq): #bremmstrahlung emission
        gff = np.abs(np.sqrt(3)*np.log(kb * T/(h * freq))/pi)
        n = rho/(mu*mh)
        return 5.44e-39 * gff * (n**2)*np.sqrt(1./T)*exp(- h * freq/(kb*T))
    
    F_1 = (jff(rhops,Ts,freq)*Area*ds_cool*0.5)/Area_norm
    F_2 = (jff(rho_cool,T_cool,freq)*Area*ds_cool*0.5)/Area_norm
    
    nixray = (wav > xray_max)
    F_1[nixray] = 0.0
    F_2[nixray] = 0.0
    
    xmodel = {}
    xmodel['shock'] = {'T': Ts,'rho':rhops,'ds':ds_cool/2.,'F':F_1,'lam':wav}
    xmodel['cool'] = {'T':T_cool,'rho':rho_cool,'ds':ds_cool/2,'F':F_2,'lam':wav}
    
    return xmodel
    

def model_Lya_line(model,wav=None,fLya=1e-2):
    """ calculates the peak flux of a Lyman alpha line based on a template spectrum,
    scaling to requested fLya of black body flux
    Parameters:
    ----------
    model: model_ object class 
    
    wav: wavelength in microns if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
            
    fLya: proportion of stellar flux that lyman alpha line should contain
    
    Returns: value of the peak Lyman alpha flux 
    """
    
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
    
    lam_nm, flya = np.loadtxt(model.models_dir + 'templates/spec_lya.txt',skiprows=4, unpack=True)
    lam_mu = lam_nm*1e-3 #microns
    nu_lya = c/ lam_mu * 1e4
    
    fnu_lya = np.interp(wav, lam_mu, flya*(h*nu_lya)* (c*1e7) / (nu_lya**2))
    
    wav, fnu_bb = calc_blackbody_spectrum(model,wav=wav)
    
    Fstar = np.trapz(fnu_bb,x=freq) #total stellar flux
    Ftot = np.trapz(fnu_lya,x=freq) #total lyman alpha flux
    norm = (Fstar*fLya)/(Ftot) #scaling so that the total in the lya line is fLya of the stellar flux
    
    fpeak_lya = np.amax(norm*fnu_lya) #peak lyman alpha flux to insert into the RT calculation

    return fpeak_lya

def calc_ISRF(model,wav=None, gnorm=None):
    """ using template ISRF measurement, normalized field by number of G0 
    
    Parameters:
    ----------
    model: model_ object class 
    
    wav: wavelength in microns if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
            
    gnorm: scaling for isrf energy by G0, if None, input model G0 is used
    
    Returns: wav, fnu of isrf in ergs/cms2/s/Hz/str
    """
    
    from scipy.interpolate import interp1d
    
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
        
    if gnorm is None:
        gnorm = model.rad['G0']
    
    templates_dir = model.models_dir+'templates/'
    fname= templates_dir+'ISRF.csv'
    
    ISRF = np.loadtxt(fname,skiprows=1,delimiter=',')
    lami = ISRF[:,0] #micron
    flam = ISRF[:,1]/(4*pi) #ergs/cm2/s/micron/str
    
    fnu_ = flam*(lami*(lami*1e-4))/c #conversion to ergs/cms2/s/Hz/str
    
    Inu = interp1d(lami, fnu_,fill_value='extrapolate')
    Inu_isrf = np.clip(Inu(wav),a_min=np.amin(fnu_)*1e-2,a_max=None)
    
    isrf_index = (wav > 0.0912) & (wav < 2.4) # wavelengths over which G0 is measured
    
    Ftot = np.trapz(Inu_isrf[isrf_index][::-1],x=freq[isrf_index][::-1])
    
    norm = gnorm*G0/Ftot
    Inu_isrf *= norm
    
    return wav, Inu_isrf

def model_viscous_dissipation(model,accrate = None):
    """ calculates the viscous dissipation due to disk accretion in the midplane
    Parameters:
    ----------
    model: model_ object class 
    
    accrate: accretion rate in Msun/yr, if None, uses model accretion rate
    
    Returns: ndarray of dissipation rate across the model grid in ergs/cm^3/s
    """
    if accrate is None:
        accrate = model.star['accrate']

    Mstar = model.star['Ms']*Msun
    Rstar = model.star['Rs']*Rsun
    
    R_CYL,Z_CYL = model.make_rz()
    rcm = R_CYL * AU
    D_disk = np.zeros_like(R_CYL)
    
    # assuming energy is generated in the midplane layers only
    midplane_index = (Z_CYL/model.H(R_CYL) < 1.)
    
    D_disk[midplane_index] = ((3.0 * Gconv * Mstar * accrate * Msun / yr ) / (4.0 * np.pi * rcm[midplane_index]**3)
                     * (1.0 - (Rstar / rcm[midplane_index])**0.5))
    return D_disk
    

###### function for binning the high energy radiation fields onto the wavelengths for chemistry 
def calc_input_spectrum(model,wav=None):
    """ computes the input spectrum in stars.inp
    on other wavelengths
    
    Parameters:
    -----------
    model: model_ class object for parameters
    
    wav: wavelength in microns if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
    
    Returns: total flux in units of ergs/s/cm^2/Hz
    
    """
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
        
    accrate = model.star['accrate']
    
    wav,fnu_bb = calc_blackbody_spectrum(model,wav = wav)
    
    wav,fnu_uv = calc_accretion_spectrum_uv(model,wav=wav,accrate=accrate)
    flya_peak = model_Lya_line(model,wav=wav,fLya=model.rad['fLya'])
    fnu_uv[np.argmin(np.abs(wav - lam_lya))] += flya_peak
    fnu_tot = fnu_uv + fnu_bb
        
    if model.rad['xray'] == True:
        wav,fnu_xray = calc_accretion_spectrum_xray(model,wav=wav,accrate=accrate)
        fnu_tot += fnu_xray
        
    return fnu_tot



def calc_powerlaw_uv(model,wav = None, puv=1.1, LUV=None,flya=0.5): ### power-law, not currently used
    """ calculates a spectrum for an uv power-law emission profile
    *** not in use ***
    
    Parameters:
    -----------
    model: model_ class object for parameters
    
    wav: wavelenghth in microns if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
            
    puv: power-law coefficient
    
    LUV: fraction of stellar flux in uv to normalize power-law, if None, then power-law will be scaled to match 
            the b-body flux at edge of uv band
            
    flya: fraction of uv emission to scale for Lya (rudimentary)
    
    
    Returns: wav, fnu, ndarrays for the frequencies and spectrum at every input wavelength in ergs/cm^2/s/Hz
    """
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
        
    freq, fnu_bb = calc_blackbody_spectrum(model,wav=wav)
    
    lya_lam = np.argmin(np.abs(wav - 1216*1e-4))
    uv_lam = (wav > uv_min) & (wav < uv_max)
    fnu_uv = (wav/wav[uv_lam][-1])**(puv+2)
    
    if LUV != None: # scale by fraction of stellar luminosity from UV band
        Fstar = np.trapz(fnu_bb,x=freq) #total flux
        Ftot = np.trapz(fnu_uv[uv_lam],x=freq[uv_lam]) #flux in UV band
        norm = ((LUV*Fstar)/Ftot)
    else: #scale to match the stellar luminosity at the edge of the UV band
        norm = np.amax(fnu_bb[uv_lam])
        
    uv_flux = fnu_uv*norm
    uv_flux[wav>uv_max] = 0.0
    
    Ftot = np.trapz(uv_flux[uv_lam][::-1],x=freq[uv_lam][::-1])
    
    dnu_lya = np.abs(np.gradient(freq))[lya_lam]
    
    uv_flux[lya_lam] = Ftot/((1./flya) - 1.)/dnu_lya #set the lyman alpha flux to be flya of the total
    
    return wav, uv_flux

def calc_powerlaw_xray(model,wav=None,TX=10e6,fX=0.01): ### power-law X-ray emission
    """ calculates a spectrum for an x-ray power-law emission profile
    Parameters:
    -----------
    model: model_ class object for parameters
    
    wav: wavelenghth in microns if None, wavelength is read in from the model's outdir
            otherwise, provided wavelengths are used for calculation
            
    TX: X-ray temperature for power-law 
    
    fX: fraction of stellar flux in x-ray to normalize power-law
    
    Returns: wav, fnu, ndarrays for the frequencies and spectrum at every input wavelength in ergs/cm^2/s/Hz
    """
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
        
    freq, fnu_bb = calc_blackbody_spectrum(model,wav=wav)

    x_lam = (wav > xray_min) & (wav < xray_max)
 
    fnu_xray = (freq[x_lam][-1]/freq)*np.exp(-(h * freq)/(kb*TX))
    fnu_xray[wav>xray_max] = 0.0
    
    Fstar = np.trapz(fnu_bb,x=freq) #total flux
    Ftot = np.trapz(fnu_xray[x_lam],x=freq[x_lam]) #flux in Xray
    
    norm = (Fstar*fX)/(Ftot) # Flux in the X-ray is fX of total flux
    
    return wav, fnu_xray*norm