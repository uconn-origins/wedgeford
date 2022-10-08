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
    
    wav: wavelength in microns if None, wavelength is read in from the model's outdir
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
        
    if model.star['xmodel'] is None:
        model_acc_onto_star(model,wav=wav,accrate=accrate)
    elif model.star['xmodel']['lam'] != wav:
        model_acc_onto_star(model,wav=wav,accrate=accrate)
    xmodel = model.star['xmodel']
    
    if fX is None:
        fnu = xmodel['F']
    else:
        wav, fnu = calc_powerlaw_xray(model,wav=wav,TX=xmodel['Tx'],fX=fX)
    fnu[wav > xray_max] = 0.0
    return wav, fnu


def model_acc_onto_star(model,wav=None,accrate=None):
    """ calculates the xray spectrum due to shocks from accretion onto the stellar surface
    based on Calvet+Gullbring 1998/Brickhouse 2013 -- models emission from cooling layer
    
    flux is bremmstrahlung power-law emission profile from APEC calculations
    
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
    Ri = 5 #magnetospheric truncation radius in units of stellar radii
    # set to 5 from Calvet&Gullbring 1998
    vs = 3.07e7*np.sqrt(model.star['Ms']/0.5)*np.sqrt(2.0/(model.star['Rs']))*np.sqrt(1.-(1./Ri))
    rhos = (accrate*Msun/yr)/(vs*Area)
    #Ts = 8.6e5*(model.star['Ms']/0.5)*(2.0/model.star['Rs'])
    Ts = (3/16.) * (mu*mh/kb) * (vs**2)
    vps = vs*0.25
    rhops = rhos*4.0
    #C = ((93*np.sqrt(3)-40*pi)/10)*np.sqrt(mh**5/kb)*(np.sqrt(2**5)/(2*(np.sqrt(2+0.25+1))))/1.64e-19 # Eq 9 Feldmeier+1997
    #based on radiative cooling for low density plasma in this regime
    #ds_cool = C*(vps**4)/(rhops)
    #from hartmann 2016
    ds_cool = vps*6.7e3*(Ts**(1.5))/(rhops/(mu*mh))
    def shock_cool(ds_cool, rho_ps, T_ps):
        x = np.logspace(-2,0,10)
        s = 1 - x
        a = 0.87225
        h = a * x**(2/7.) * ( 1. + ((1/a) - 1.) * x**(2/7.))
        rho_cool = (1./h)*rho_ps
        T_cool = ((1/3.) * h * ( 4. - h ))*T_ps
        ds_x = np.abs(np.gradient(s))*ds_cool # for volume calculation
        return ds_x, T_cool, rho_cool
    
    #def jff(rho,T,freq): # optically thin bremmstrahlung emission
        #gff = np.abs(np.sqrt(3)*np.log(kb * T/(h * freq))/pi)
        #n = rho/(mu*mh)
        #return 5.44e-39 * gff * (n**2)*np.sqrt(1./T)*exp(- h * freq/(kb*T))
    ds_x_grid, T_cool_grid, rho_cool_grid = shock_cool(ds_cool, rhops, Ts)
    F_grid = np.array([calc_ee_brems(freq,T_cool,rho_cool) * ds_x 
              for T_cool, rho_cool,ds_x in zip(T_cool_grid,rho_cool_grid,ds_x_grid)])
    F_total = np.sum(np.array(F_grid),axis=0)*Area/Area_norm
    #F_1 = (jff(rhops,Ts,freq)*Area*ds_cool*0.5)/Area_norm
    #F_2 = (jff(rho_cool,T_cool,freq)*Area*ds_cool*0.5)/Area_norm
    
    xmodel = {'Tx': np.amax(T_cool_grid), 'F': F_total, 'ds': ds_cool, 'lam':wav}
    model.star['xmodel'] = xmodel
    

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
    """ calculates a spectrum for an x-ray power-law emission profile, 
    adapted from https://github.com/AtomDB/pyatomdb/blob/master/pyatomdb/pyatomdb/apec.py
    
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


def calc_ee_brems(freq, Te, rho):
    """
    calculate the electron-electron bremsstrahlung 
    (taken straight from APEC models)
    Parameters
    ----------
    freq : array (float)
    frequencies in Hz
    T : float
    Electron temperature in K
    rho : float
    gas density (cm^-3)
    Returns
    -------
    array(float)
    ee_brems in photons cm^3 s^-1 keV-1 at each point E.
    This should be multiplied by the bin width to get flux per bin.
    References
    ----------
    Need to check this!
    """
    #
    #  T is the electron temperature (in keV)
    #  N is the electron density (in cm^-3)
    #  E is the list of energies (in keV)
    #
    from scipy import integrate
    
    T = (kb*Te)/keV #electron temperature in keV
    E = h*freq/keV # energy vals in keV
    N = rho/(mu*mh) # electron number density
    
    aI1 = np.array([(3.15847E+0, -2.52430E+0, 4.04877E-1, 6.13466E-1, 6.28867E-1, 3.29441E-1),
             (2.46819E-2, 1.03924E-1, 1.98935E-1, 2.18843E-1, 1.20482E-1, -4.82390E-2),
             (-2.11118E-2, -8.53821E-2, -1.52444E-1, -1.45660E-1, -4.63705E-2, 8.16592E-2),
             (1.24009E-2, 4.73623E-2, 7.51656E-2, 5.07201E-2, -2.25247E-2, -8.17151E-2),
             (-5.41633E-3, -1.91406E-2, -2.58034E-2, -2.23048E-3, 5.07325E-2, 5.94414E-2),
             (1.70070E-3, 5.39773E-3, 4.13361E-3, -1.14273E-2, -3.23280E-2, -2.19399E-2),
             (-3.05111E-4, -7.26681E-4, 4.67015E-3, 1.24789E-2, -1.16976E-2, -1.13488E-2),
             (-1.21721E-4, -7.47266E-4, -2.20675E-3, -2.74351E-3, -1.00402E-3, -2.38863E-3),
             (1.77611E-4, 8.73517E-4, -2.67582E-3, -4.57871E-3, 2.96622E-2, 1.89850E-2),
             (-2.05480E-5, -6.92284E-5, 2.95254E-5, -1.70374E-4, -5.43191E-4, 2.50978E-3),
             (-3.58754E-5, -1.80305E-4, 1.40751E-3, 2.06757E-3, -1.23098E-2, -8.81767E-3)])

    # column j=0-5; line i=0-10
    aI2 = np.array([(-1.71486E-1, -3.68685E-1, -7.59200E-2, 1.60187E-1, 8.37729E-2),
             (-1.20811E-1, -4.46133E-4, 8.88749E-2, 2.50320E-2, -1.28900E-2),
             (9.87296E-2, -3.24743E-2, -8.82637E-2, -7.52221E-3, 1.99419E-2),
             (-4.59297E-2, 5.05096E-2, 5.58818E-2, -9.11885E-3, -1.71348E-2),
             (-2.11247E-2, -5.05387E-2, 9.20453E-3, 1.67321E-2, -3.47663E-3),
             (1.76310E-2, 2.23352E-2, -4.59817E-3, -8.24286E-3, -3.90032E-4),
             (6.31446E-2, 1.33830E-2, -8.54735E-2, -6.47349E-3, 3.72266E-2),
             (-2.28987E-3, 7.79323E-3, 7.98332E-3, -3.80435E-3, -4.25035E-3),
             (-8.84093E-2, -2.93629E-2, 1.02966E-1, 1.38957E-2, -4.22093E-2),
             (4.45570E-3, -2.80083E-3, -5.68093E-3, 1.10618E-3, 2.33625E-3),
             (3.46210E-2, 1.23727E-2, -4.04801E-2, -5.68689E-3, 1.66733E-2)])
    # column j=6-10, line 0-10

    # Region II (1 keV<=k_B<=300 keV)
    abII = np.array([(-3.7369800E+1, -9.3647000E+0, 9.2170000E-1, -1.1628100E+1, -8.6991000E+0),
              (3.8036590E+2, 9.5918600E+1, -1.3498800E+1, 1.2560660E+2, 6.3383000E+1),
              (-1.4898014E+3, -3.9701720E+2, 7.6453900E+1, -5.3274890E+2, -1.2889390E+2),
              (2.8614150E+3, 8.4293760E+2, -2.1783010E+2, 1.1423873E+3, -1.3503120E+2),
              (-2.3263704E+3, -9.0730760E+2, 3.2097530E+2, -1.1568545E+3, 9.7758380E+2),
              (-6.9161180E+2, 3.0688020E+2, -1.8806670E+2, 7.5010200E+1, -1.6499529E+3),
              (2.8537893E+3, 2.9129830E+2, -8.2416100E+1, 9.9681140E+2, 1.2586812E+3),
              (-2.0407952E+3, -2.9902530E+2, 1.6371910E+2, -8.8818950E+2, -4.0474610E+2),
              (4.9259810E+2, 7.6346100E+1, -6.0024800E+1, 2.5013860E+2, 2.7335400E+1)])
    # column a_0j-a2j, b0j,b1j; line j=0-8

    cII = np.array([(-5.7752000E+0, 3.0558600E+1, -5.4327200E+1, 3.6262500E+1, -8.4082000E+0),
             (4.6209700E+1, -2.4821770E+2, 4.5096760E+2, -3.1009720E+2, 7.4792500E+1),
             (-1.6072800E+2, 8.7419640E+2, -1.6165987E+3, 1.1380531E+3, -2.8295400E+2),
             (3.0500700E+2, -1.6769028E+3, 3.1481061E+3, -2.2608347E+3, 5.7639300E+2),
             (-3.2954200E+2, 1.8288677E+3, -3.4783930E+3, 2.5419361E+3, -6.6193900E+2),
             (1.9107700E+2, -1.0689366E+3, 2.0556693E+3, -1.5252058E+3, 4.0429300E+2),
             (-4.6271800E+1, 2.6056560E+2, -5.0567890E+2, 3.8008520E+2, -1.0223300E+2)])
    # column CII_2j-cII_6j; line j=0-6

    # Region III (300 keV<= k_BT<=7 MeV)
    abIII = np.array([(5.2163300E+1, 4.9713900E+1, 6.4751200E+1, -8.5862000E+0, 3.7643220E+2),
               (-2.5703130E+2, -1.8977460E+2, -2.1389560E+2, 3.4134800E+1, -1.2233635E+3),
               (4.4681610E+2, 2.7102980E+2, 1.7414320E+2, -1.1632870E+2, 6.2867870E+2),
               (-2.9305850E+2, -2.6978070E+2, 1.3650880E+2, 2.9654510E+2, 2.2373946E+3),
               (0.0000000E+0, 4.2048120E+2, -2.7148990E+2, -3.9342070E+2, -3.8288387E+3),
               (7.7047400E+1, -5.7662470E+2, 8.9321000E+1, 2.3754970E+2, 2.1217933E+3),
               (-2.3871800E+1, 4.3277900E+2, 5.8258400E+1, -3.0600000E+1, -5.5166700E+1),
               (0.0000000E+0, -1.6053650E+2, -4.6080700E+1, -2.7617000E+1, -3.4943210E+2),
               (1.9970000E-1, 2.3392500E+1, 8.7301000E+0, 8.8453000E+0, 9.2205900E+1)])
    # column aII_0j-aIII_2j, bIII_0j, bIII_1j; line j=0-8
    #print "E:",E
    aI = np.hstack((aI1,aI2))

    inum1 = np.arange(11)
    jnum = np.resize(inum1,(11,11))
    inum = np.transpose(jnum)
    aII = np.zeros((9,3))
    bII = np.zeros((9,2))
    [aII, bII] = np.hsplit(abII,np.array([3]))

    numII = np.arange(9)
    aIIj = np.transpose(np.resize(numII,(3,9)))
    aIIi = np.resize(np.arange(3),(9,3))
    bIIj = np.transpose(np.resize(np.arange(9),(2,9)))
    bIIi = np.resize(np.arange(2),(9,2))
    cIIj = np.transpose(np.resize(np.arange(7),(5,7)))
    cIIi = np.resize(np.arange(2,7),(7,5))

    aIII = np.zeros((9,3))
    bIII = np.zeros((9,2))

    kTIII = 1.e3 #1 MeV, in unit of keV
    [aIII, bIII] = np.hsplit(abIII,np.array([3]))


    tao = T/510.0

    #Earray, Eisvec =  util.make_vec(E)
    Earray = np.array(E)
    
    x = Earray/T
    numx=len(Earray)

    GI = np.zeros((numx,))
    AIIr = np.zeros((numx,))
    BIIr = np.zeros((numx,))
    GpwII = np.zeros((numx,))
    GII = np.zeros((numx,))
    FCCII = np.zeros((numx,))
    Ei0 = np.zeros((numx,))
    GpwIII = np.zeros((numx,))
    
    if T<0.05:
        ret = np.zeros(len(x), dtype=float)
    
    elif 0.05<=T<70.:
        # hmm
        GI=np.zeros(len(x))
        theta = (1/1.35) * ( np.log10(tao) + 2.65)
        bigx = (1/2.5) * (np.log10(x) + 1.50)
        for i in range(11):
            for j in range(11):
                GI += aI[i,j]*(theta**i)*(bigx**j)
        GI *= np.sqrt(8/(3*np.pi))
        ret = 1.455e-16*N**2*np.exp(-x)/(x*np.sqrt(tao))*GI
    elif 70.<=T<300.:
        taoII = tao
        for k in range(numx):
            def integrand(t):
                return np.exp(-1.0*t)/t
            [Ei0[k,],error] = scipy.integrate.quad(integrand,x[k,],\
                                             np.Inf,args=())
            AIIr[k,] = np.sum(aII*taoII**(aIIj/8.)*x[k,]**(aIIi))
            BIIr[k,] = np.sum(bII*taoII**(bIIj/8.)*x[k,]**(bIIi))
            FCCII[k,] = 1.+np.sum(cII*taoII**(cIIj/6.)*x[k,]**(cIIi/8.))
            GpwII[k,] = np.sum(aII*taoII**(aIIj/8.)*x[k,]**(aIIi))-\
                          np.exp(x[k,])*(-1.0)*Ei0[k,]*\
                          np.sum(bII*taoII**(bIIj/8.)*x[k,]**(bIIi))
            GII[k,] = GpwII[k,]*FCCII[k,]
        ret = 1.455e-16*N**2*np.exp(-x)/(x*np.sqrt(taoII))*GII
    elif 300.<=T<7000.:
        taoIII = tao
        for k in range(numx):
            GpwIII[k,] = np.sum(aIII*taoIII**(aIIj/8.)*x[k,]**(aIIi))-\
                       np.exp(x[k,])*(-1.0)*Ei0[k,]*\
                       np.sum(bIII*taoIII**(bIIj/8.)*x[k,]**(bIIi))
        ret = 1.455e-16*N**2*np.exp(-x)/\
               (x*np.sqrt(taoIII))*GpwIII
    else:
        taoIV = tao
        GIV = 3./(4.*np.pi*np.sqrt(taoIV))*\
              (28./3.+2.*x+x**2/2.+2.*(8./3.+4.*x/3.+x**2)*\
              (np.log(2.*taoIV)-0.57721)-np.exp(x) \
              *(-1.0*Ei0)*(8./3.-4.*x/3.+x**2))

        ret = 1.455e-16*N**2*np.exp(-x)/(x*np.sqrt(taoIV))*GIV
    ph_keV = ret/150
    #convert to ergs/cm^3/s/Hz
    E_Hz = ph_keV * E*keV * (h/keV)
    return E_Hz