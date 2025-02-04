from pylab import *
#import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
import pickle
from . units import *
from . write_radmc_files import *



class out: #Class that handles and stores the outputs and data for radmc3d after you've made the initial model structure
    """ 
    Attributes
    --------------
    m : wedgeford object model class see: models/make_model.py
    
    r : 1-d array, radial coordinate  in au
    
    theta : 1-d array, polar coordinate
    
    phi : 1-d array, azimuthal coordinate

    rho   : dictionary
            densities in g/cm^3 
            keys: 'dust'+nspec, 'gas'
    T : dictionary
            temperatures in K
            keys: 'dust', 'gas', 'shock', 'UV'
            
    J: nested dictionary
          field = 'xray' or 'uv'
          stores dictionary:
              'J_phot' : ndarray of integrated mean intensity in nphotons/cm^2/s
              'e_phot' : mean energy of photons in integrated wavelength range
              
    wav, freq: wavelengths in microns, frequencies in hertz in local wavelength_micron.inp file
                
    """
    def __init__(self,model):
        if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
        self.m = model #model object is encoded here
        self.r = model.r
        self.theta= model.theta
        self.phi = model.phi
        # dictionaries for storing component quantities
        self.rho = {}
        self.T = {} 
        self.J = {}
        self.wav, self.freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    
    def _fieldReader(self, fname='', ndim=3, nheader = None):
        """Reads a scalar field from file.

        Parameters
        ----------

        fname  : str
                Name of the file containing a scalar variable

        ndim   : int
                Number of dimension of the data field (3 for gas variables, 4 for dust)
                
        nheader : optional
                specify the header length in lines for files with alternate header lengths
                
        Returns
        -------

        Returns a numpy Ndarray with the scalar field
        """

        nx,ny,nz = self.m.grid['N']
        data = None
        with open(fname, 'r') as rfile:
            if ndim == 3:
                if nheader is None:
                    nheader = 2
                hdr = np.fromfile(rfile, count=nheader, sep=" ", dtype=np.int64)
            else:
                if nheader is None:
                    nheader = 3
                hdr = np.fromfile(rfile, count=nheader, sep=" ", dtype=np.int64)

            if hdr[1] != (nx * ny * nz):
                msg = 'Number of grid cells in ' + fname + ' is different from that in amr_grid.inp '\
                              + ' nr cells in ' + fname + ' : ' + ("%d" % hdr[1]) + '\n '\
                              + ' nr of cells in amr_grid.inp : '\
                              + ("%d" % (nx * ny * nz))
                raise ValueError(msg)

            data = np.fromfile(rfile, count=-1, sep=" ", dtype=np.float64)
        if ndim == 3:
            if data.shape[0] == hdr[1]:
                data = np.reshape(data, [1, nz, ny, nx])
            else:
                msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                raise ValueError(msg)
        else:
            if data.shape[0] == hdr[1] * hdr[2]:
                data = np.reshape(data, [hdr[2], nz, ny, nx])
            elif data.shape[0] == hdr[1] * hdr[2] + hdr[2]:
                data = np.reshape(data[hdr[2]:], [hdr[2], nz, ny, nx]) #for mean intensity files
            else:
                msg = 'Internal inconsistency in data file, number of cell entries is different from ' \
                                  'indicated in the file header'
                raise ValueError(msg)
                        
        data = np.swapaxes(data, 0, 3)
        data = np.swapaxes(data, 1, 2)
        return data
    
    
    def read_rho(self,fname_dust=None, fname_gas=None):
        """Reads the density of components, stores the 3D quantity in self.rho keyed by the component name

        Parameters
        ----------

        fname_dust : str, optional
                  Name of the file that contains the dust density. If omitted 'dust_density.inp' is used
        fname_gas : str, optional
                  Name of the file that contains the gas density. If omitted 'gas_density.inp' is used
                  
        """
        model = self.m
        if fname_dust is None:
            fname_dust = model.outdir+'dust_density.inp'
        if fname_gas is None:
            fname_gas =  model.outdir+'gas_density.inp'
        
        if os.path.exists(fname_dust):
            print('Reading '+fname_dust)
            rhodust = self._fieldReader(fname=fname_dust,ndim=4)
            self.ndust = np.shape(rhodust[0,0,0,:])[-1]
            for n in range(self.ndust):
                self.rho['dust' + str(n+1)] = rhodust[:,:,:,n]
        else:
            print('No dust_density.inp found, loading from model parameters')
            self.rho['dust1'] = model.rho_embedded(fluid=1).swapaxes(0,1)
            self.rho['dust2'] = model.rho_embedded(fluid=2).swapaxes(0,1)
        
        if os.path.exists(fname_gas):
            print('Reading '+fname_gas)
            rhogas = self._fieldReader(fname=fname_gas,ndim=3)
            self.rho['gas'] = rhogas.squeeze()
        else:
            self.rho['gas'] = model.rho_embedded(fluid=0).swapaxes(0,1)
            
        return True
    
    def calc_rho2D(self,fluid='dust'):
        """ calculates 2D azimuthal average of the requested component
        
        Parameters
        -----------
        
        fluid: str, model component requested
               'dust' : total dust density
               'dust' + str(number) :  dust density of dust population number 1, 2, ...
               'gas' : the gas density
        
        Returns: 2D numpy array averaged along the azimuth
        
        """
        
        # read in the data if it hasn't been already
        if self.rho == {}:
            self.read_rho()

        if fluid == 'dust':
            rho_tot = np.zeros_like(self.rho['gas'])
            for n in range(self.ndust):
                rho_tot += self.rho['dust' + str(n+1)]
        else:
            rho_tot = self.rho[fluid]
            
        return np.average(rho_tot,axis=-1)
    

            
    def read_Tdust(self,fname=None):
        """Reads the temperature of components, stores the 3D quantity in self.T keyed by 'dust'
        
        Parameters
        ----------
        fname : str, optional
                  Name of the file that contains the dust temperature. If omitted 'dust_temperature.dat' is used.
        
        """
        model = self.m
        if fname is None:
            fname = model.outdir + 'dust_temperature.dat'
        
        if os.path.exists(fname):
            print('Reading '+fname)
            Tdust = self._fieldReader(fname=fname, ndim=4)
        
            #averages dust temperature across all dust populations
            self.T['dust'] = np.average(Tdust,axis=-1)

            return np.average(Tdust,axis=-1)
        else:
            self.T['dust'] = None
            print('ERROR: no dust temperature file found at:', fname)
            return None
    
    
    def calc_T2D(self,fluid='dust'):
        """ calculates 2D azimuthal average of the requested component
        
        Parameters
        -----------
        
        fluid: str, model component requested
               'dust' : total dust density
               'gas' : the gas density
        
        Returns: 2D numpy array averaged along the azimuth
        
        """
        # read in the data if it hasn't been already
        
        if fluid not in self.T.keys():
            if fluid == 'dust':
                self.read_Tdust()
            if fluid == 'gas' or fluid =='uv':
                self.calc_Tgas(ndim=3)
            if fluid == 'shock':
                self.model_Tshock()
        Tvar = self.T[fluid]
        if Tvar is not None:
            return np.average(self.T[fluid],axis=-1)
        else:
            return None
    
        
    def read_Jnu(self,fname=None):
        """Reads the mean intensity.out file 
        Returns: arrays of wavelengths in microns, frequencies in Hz, and 4 dimensional ndarray of the mean intensity at each frequency
        
        Parameters:
        ------------
        fname: str, optional. If none specified mcmono_wavelength_micron.inp is used.
        """
        model = self.m
        
        if fname is None:
            fname = model.outdir + 'mcmono_wavelength_micron.inp'
        wav,nu = read_wavelength(fname=fname)
        #print('Reading '+fname)
        
        fname = model.outdir + 'mean_intensity.out'
        if os.path.exists(fname):
            Jnu = self._fieldReader(fname=fname, ndim=4)
            return wav, nu, Jnu
        else:
            print('ERROR: no radiation field found at :', fname)
            return wav, nu, None
        
    
    def calc_Jint(self,field='uv'):
        """ calculates the integrated mean intensity of the frequency dependent mean intensity field
        output in two fields: 
        'J_phot': the number of photons/cm^2/s
        'e_phot': the average photon energy
        
        Multiplying the two gives you the mean intensity in ergs/cm^2/s
        Parameters
        ------------
        
        field: str, 'uv' or 'xray' for whatever field one wants to integrate 
       
        
        """
        #read in the data
        wav, nu, Jnu = self.read_Jnu()
        if Jnu is None:
            self.J[field] = {}
            self.J[field]['J_phot'] = None
            self.J[field]['e_phot'] = None
            return None
        else:
            wav_hires, nu_hires = read_wavelength(self.m.outdir+'wavelength_micron.inp')
            model = self.m
            if field == 'uv' or field == 'UV':
                nu_index = np.where((wav < uv_max) & (wav > uv_min))
                nu_index_hires = np.where((wav_hires < uv_max) & (wav_hires > uv_min))
            elif field == 'xray' or field == 'Xray':
                nu_index = np.where((wav < xray_max) & (wav > xray_min))
                nu_index_hires = np.where((wav_hires < xray_max) & (wav_hires > xray_min))

            freq = nu[nu_index]
            wav = wav[nu_index]
            efreq = h*freq
            Jfreq = Jnu[:,:,:,nu_index]
            Jphot = Jnu[:,:,:,nu_index]/efreq

            # calculate the relative error you get from summing over a limited number of wavelength points using the original input spectrum as a guide
            fnu0_1d = calc_input_spectrum(model,wav=wav)
            fnu0_1d_hires = calc_input_spectrum(model,wav=wav_hires[nu_index_hires])

            f0_tot = np.trapz(fnu0_1d, x = freq)*-1.
            f0_tot_hires = np.trapz(fnu0_1d_hires, x = nu_hires[nu_index_hires])*-1.
            err_sum = f0_tot/f0_tot_hires - 1.

            J_e = np.trapz(Jfreq, x=np.expand_dims(freq,(0,1,2)),axis=-1)*-1 
            J_n = np.trapz(Jphot, x=np.expand_dims(freq,(0,1,2)),axis=-1)*-1

            self.J[field] = {}
            self.J[field]['J_phot'] = J_n.squeeze()*(1.+ err_sum) #apply the correction factor
            self.J[field]['e_phot'] = (J_e/J_n).squeeze()

            return True

    def model_PDR(self):
        """ calculates the PDR temperature grid for cross-matching to the gas temperatures
        stores results in dictionary T['pdr'] = {'n', 'G', 'Ts'}
       
        'n': log of atomic hydrogen density
        'G': log of FUV field in units of G0
        'Ts': corresponding 2D array of PDF surface temperatures
        
        """
        from scipy import interpolate
        model = self.m
        
        fname = model.models_dir+'templates/tgas.npy'
        if os.path.exists(fname):
            tpdr = np.load(fname)
            log_n = np.cumsum(np.ones(49)*0.125) + 0.75 + 0.125
            log_F = np.linspace(-3.29588079,-3.29588079+ 57*0.125, 57) - np.log10(G0) #first axis
            self.T['pdr'] = {'n':log_n, 'G':log_F, 'Ts': tpdr}
            return True
        else:
            print('ERROR: templates not found')
            self.T['pdr'] = None
            return None
    
    def model_Tshock(self):
        """ calculates the shock temperatures from the infall model
        see: model_.streamline for more info
        """
        model = self.m
        Tshock = model.solve_envelope(prop='Tg')*model.stream_mask
        self.T['shock'] = Tshock.swapaxes(0,1)
        return True
    
    def calc_Tgas(self,ndim=3):
        """ calculates the gas temperature based on the UV radiation field
        (and models shock temperatures if shocks are included)
        Parameters
        ------------
        ndim: int
            2: calculates the gas temperature using two dimensional azimuthal averages (faster)

            3: calculates the gas temperature using the fully 3D arrays

        """
        model = self.m
        shock = model.env['shock']

        if self.rho == {}:
            self.read_rho()

        if self.T == {}:
            self.read_Tdust()

        if 'pdr' not in self.T.keys():
            self.model_PDR()

        if 'uv' not in self.J.keys():
            print('integrating uv field')
            self.calc_Jint(field='uv')
        
        if self.T['dust'] is None:
            self.T['gas'] = None
            return None
        if self.J['uv']['J_phot'] is None:
            self.T['gas'] = None
            return None
        
        if ndim == 3:
            nH = np.log10(2*self.rho['gas']/(mu*mh))
            T_dust = self.T['dust'].copy()
            nG0 = np.log10(4*pi*self.J['uv']['J_phot']*self.J['uv']['e_phot']/G0)  
        else:
            nH = np.log10(2*self.calc_rho2D('gas')/(mu*mh))
            T_dust = self.calc_T2D('dust')
            nG0 = np.log10(np.average(4*pi*self.J['uv']['J_phot']*self.J['uv']['e_phot']/G0,axis=-1))


        T_gas = T_dust
        if self.T['pdr'] is not None:
            G_ = self.T['pdr']['G']
            n_ = self.T['pdr']['n']
            Ts = self.T['pdr']['Ts']

            n_min = np.amin(n_)
            n_max = np.amax(n_)
            n_num = len(n_)

            g_min = np.amin(G_)
            g_max = np.amax(G_)
            g_num = len(G_)

            i = np.clip((n_num*(nH - n_min)/(n_max-n_min) + 0.5).astype(int),a_min = 0, a_max=n_num-1)
            j = np.clip((g_num*(nG0 - g_min)/(g_max-g_min) + 0.5).astype(int),a_min=0, a_max=g_num-1)

            T_UV = Ts[j,i]
            self.T['uv'] = T_UV

            if shock == True and 'shock' not in self.T.keys():
                self.model_Tshock()
            elif shock == True and self.T['shock'] is not None:
                T_shock = self.T['shock']
            else:
                T_shock = np.zeros_like(T_UV)


            Tgas_max = np.maximum(T_UV,T_shock)

            #T_crit = 130*(10**(nH)/1e10)**(0.3)
            #uncoupled = np.where(((Tgas_max/T_crit) >= 1.1))
            #threshold for CO dissociation, below this value molecular self-shielding may kick in
            uncoupled = np.where(nG0-nH > - 6)

            T_gas[uncoupled] = Tgas_max[uncoupled]

            if ndim == 3:
                self.T['gas'] = T_gas
            else:
                self.T['gas']= np.repeat(np.expand_dims(T_gas,axis=-1),len(self.m.phi),axis=-1)
            return
        else:
            self.T['gas'] = None
            return None

        
    def make_rz(self):
        """ generates 2d arrays of r-z coordinate grid (cylindrical)
        Returns: X,Z 
        """
        R,TH = np.meshgrid(self.r,self.theta)
        X = R*np.sin(TH)
        Z = R*np.cos(TH)
        return X,Z
        
def update_model(output,**updated_params):
    """ updates the model parameters to save to the current directory as the .pkl file
    """
    model = output.m
    if os.path.exists(model.outdir+'pars.pkl'):
        params = pickle.load( open( model.outdir+'pars.pkl', "rb" ))
    else:
        print('did not find current param file: {}, using current model params'.format(model.outdir+'pars.pkl'))
        params = model.print_params()
    for key in updated_params.keys():
        if key in params.keys():
            params[key] = updated_params[key]
    pickle.dump(params,open(model.outdir + 'pars.pkl','wb'))
    
def overwrite_model(output,outdir=None):
    """ writes new .inp files for the current model parameters
    !! Note: this option will overwrite current .inp files in the directory !!
    
    Files that you've been using (but not generating through the model) will be lost
    (except for opacities)
    
    Parameters:
    -----------
    output: object,  output class see outputs.py
    
    outdir: optional, default None
            if not None, changes model.outdir and writes .inp files for thermal transfer in new directory
            
    """
    model = output.m
    
    set_params = model.print_params()
    
    if outdir is not None:
        model.outdir = outdir
    
    update_model(output,**set_params)
    
    write_grid(model)
    write_wavelength(model,wav=output.wav)
    write_star(model,wav=output.wav)
    
    write_dust_density(model)
    write_gas_density(model)
    
    if model.rad['viscous_heating'] != False:
        write_viscous_heatsource(model)
        
    if model.rad['G0'] > 0:
        write_external_radfield(model,wav=output.wav)
    
    if 'gas' in output.T.keys():
        write_gas_temperature(model,Tgas=output.T['gas'])
        
    return

def check_spatial_files(output,fname_dens=None,ndim=4):
    """ checks that spatial grid and densities have the same number of cells, 
    returns None if cell amounts do not match.
    
    """
    model = output.m
    fname_grid = model.outdir + 'amr_grid.inp'
    if fname_dens is None:
        fname_dens = model.outdir + 'dust_density.inp'
    data = np.fromfile(fname_grid, count=-1, sep=" ", dtype=np.float64)
    hdr = np.array(data[:10], dtype= int)
    data = data[10:]

    # Check the file format
    if hdr[0] != 1:
        msg = 'Unknown format number in amr_grid.inp'
        raise RuntimeError(msg)

    # Check the coordinate system
    if (hdr[2] >= 100) & (hdr[2] < 200):
        crd_sys = 'sph'
    else:
        raise ValueError('non-spherical coordinate system in' + fname_grid + ' file.')

    # Get the number of cells in each dimension of spatial grid
    nx = hdr[7]
    ny = hdr[8]
    nz = hdr[9]
    
    if os.path.exists(fname_dens):
        rhodust = output._fieldReader(fname=fname_dens,ndim=ndim)
        nx_,ny_,nz_ = np.shape(rhodust[:,:,:,0])
        if nx*ny*nz != nx_*ny_*nz_:
            print('density file has different number of cells than grid file!')
            return None
        if nx*ny*nz != model.grid['N'][0]*model.grid['N'][1]*model.grid['N'][2]:
            print('grid file does not match input model grid')
            return None
        else:
            return 0
    else:
        print('file not found:' + fname_dens)
        return None
    


def prep_thermal_transfer(output,nphot=500000,mrw=1,maxtau=5):
    """ checks for files based on physics specified in the input files 
        and generates the necessary files for a radmc3d mctherm run
    Note: presence of some main input files in the directory means they won't be re-written 

    Parameters
    -----------
    output: output class object to generate files for

    nphot: radmc3d parameter, number of photons for main mc calculation
    
    mrw: radmc3d parameter, turning modified random walk on/off 

    maxtau: radmc3d parameter, maximum scattering optical depth before absorption
    
    Returns:
    --------
    number of errors flagged in input files
    """
    
    model = output.m
    
    err = 0
    fname = model.outdir+'wavelength_micron.inp'
    wav,freq = read_wavelength(fname)
    output.wav = wav
    output.freq = freq
    
    ## if xray is True, the model wavelength grid is recalculated to be set on the same grid of xray   wavelengths
    if model.rad['xray'] == True:
        wav_xray = np.logspace(np.log10(xray_min),np.log10(xray_max), 20)
        wav_update = np.append(wav_xray,wav[wav>xray_max])
        #recalculate spectrum on new wavelength range
        write_wavelength(model,wav=wav_update)
        write_star(model,wav=wav_update)
        output.wav = wav_update
        output.freq = c / wav_update * 1e4
    
    file_list = ['amr_grid.inp', 'dust_density.inp','stars.inp']
    func_list = [write_grid, write_dust_density, write_star]
    
    for file, func in zip(file_list,func_list):
        if os.path.exists(model.outdir+ file) != True:
            func(model)
    
    val = check_spatial_files(output)
    if val is None:
        err += 1
    
    if os.path.exists('heatsource.inp') != True and model.rad['viscous_heating'] != False:
        file_list.append('heatsource.inp')
        write_viscous_heatsource(model)
        
        
    if os.path.exists('dustopac.inp') != True:
        file_list.append('dustopac.inp')
        write_opacities(model,update=False)
        
        
    if os.path.exists('external_source.inp') != True and model.rad['G0'] > 0:
        file_list.append('external_source.inp')
        write_external_radfield(model,wav=wav)
        
    file_list.append('radmc3d.inp')    
    write_main(model,scat=2,nphot=nphot,mrw=mrw,maxtau=maxtau)
    
    return err

def do_thermal_transfer(output,nt=4,prep=False,**prepkw):
    """ starts thermal transfer from python script
    
    Parameters
    -----------
    output: output class object to generate files for

    nt: int, number of threads to run on

    prep: boolean, if True prep_thermal_transfer will be run with **prepkw
    
    """
    model = output.m
    if prep == True:
        err = prep_thermal_transfer(output,**prepkw)
        if err > 0:
            print('errors in preparation of thermal transfer!')
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    os.system('radmc3d mctherm setthreads {}'.format(nt))
        
def prep_he_transfer(output,nx=5,nu=5):
    """ generates the necessary files for a radmc3d mcmono run for high energy radiative transfer
    writes mcmono_wavelength_micron.inp
    updates xray opacities using tables from BB2011
    
    
    Parameters
    -----------
    output: output class object to generate files for
    
    nx: number of xray wavelengths for transfer, no effect if model.rad['xray'] = False
    
    nu: number of uv wavelengths for transfer

    """
    model = output.m
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
        
    x_lam = np.logspace(np.log10(xray_min),np.log10(xray_max),nx)
    u_lam = np.append(uv_min,np.linspace(lam_lya,uv_max,nu-1)) #includes Lya line
    
    if model.rad['xray'] == True:
        he_lam = np.append(x_lam,u_lam)
        write_opacities(model,update=True)
    else:
        he_lam = u_lam
        
    write_wavelength(model,wav=he_lam,fname=model.outdir+'mcmono_wavelength_micron.inp')
    write_gas_density(model)
    write_main(model, scat=2, mrw=1, maxtau=20)
    
def do_he_transfer(output,nphot=100000,nt=4,prep=False,**prepkw):
    """ starts high energy radiative transfer from python script
    
    Parameters
    -----------
    output: output class object to generate files for
    
    nphot: number of photons for mcmono calculation

    nt: int, number of threads to run on

    prep: boolean, if True prep_he_transfer will be run with prepkw
    
    """
    model = output.m
    if prep == True:
        prep_he_transfer(output,**prepkw)
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    os.system('radmc3d mcmono nphot_mono {} setthreads {}'.format(nphot,nt))
        
        
def prep_line_transfer(output,molecules={'names':['co'],'abundances':[1e-4],'lines':[2]},teq=0):
    """ generates the necessary files for radmc3d line transfer
    downloads linelist files
    writes molecular number density files based on abundances
    writes velocity and temperature files
    
    Parameters
    -----------
    output: output class object to generate files for
    
    molecules: dict, entries of 'names', 'abundances, 'lines' for every molecule being computed
    
    teq: 0 or 1, if 0, temperature of gas != temperature of dust
    
    """
    model = output.m
    get_molecule_info(model,names=molecules['names'])
    write_molecule_density(model, names=molecules['names'],abundances = molecules['abundances'])
    write_lines(model)
    write_velocities(model)
    if teq == 0:
        if 'gas' not in output.T.keys():
            output.calc_Tgas(ndim=3)
        write_gas_temperature(model,output.T['gas'])
    write_main(model, scat=1, mrw=1, maxtau=20,dust=0, lines=1, teq=teq)
    
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
