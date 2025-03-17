from pylab import *
#import radmc3dPy as rpy
import numpy as np
import os

from .. units import *
from .. outputs import *
from scipy.ndimage import gaussian_filter


class chemdisk:
    
    """ class object to contain information for the chemistry model
    Attributes:
    -------------
    input:  object, output class containing data structures from radmc3d
    
    chemdir: the environ directory for the chemical code
    
    rundir: the runs directory for the outputs of the chemical code
    
    data: dict, keyed by properties to be written as inputs to the chemical code
    
    inpfiles: dict, contents of the main input files to be written. 
                filled by set_X() functions where X is the number of the input file 
                

    Parameters:
    -----------
    output: object, output class see outputs.py 
    
    chemdir: the name of the directory for the model, folders will be made in 
            chemistry/environ/[chemdir] and chemistry/runs/[chemdir] if one does not exist
            
    rgrid: optional, None or array with the cell faces of the radial grid dimensions in au
            if None, default grid will be logarithmically spaced with half the number of points as radmc3d grid
    
    zgrid: optional, None or array with the cell faces of the vertical grid dimensions in z/r
         if None, default grid will be logarithmically spaced bounded by envelope cavity and have 50 cells.

    iodict: dict, optional
            can use this to change the default .inp files to use as templates for running the chemistry code and the values within them
            if you run set_X() commands with new parameter dictionaries, they will automatically update this dictionary
            keys are the integer of the setup file, while values are dictionaries with keys matching the new flags
    """
    
    def __init__(self,output,chemdir = '/test1/',rgrid = None, zgrid = None, iodict= {}):   

        self.input = output
        model = output.m
        chemdir = chemdir.strip('/')
        self.parentchemdir = '/'.join([self.input.m.parent_dir ,'chemistry'])
        self.rundir = '/'.join([self.parentchemdir, 'runs' ,chemdir]) + '/'
        self.chemdir = self.rundir + 'environ/'
        if os.getcwd() != self.input.m.outdir:
            os.chdir(self.input.m.outdir)
            
        self.iodict = {}
        self.iodict['0'] = {'reactions':'herb0308gr','species':'herb0308gr','name':'', 'abund0':None}
        self.iodict['2'] = {'first':None,'last':None,'total':None}
        self.iodict['3'] = {'name':'','newvals': {}}
        self.iodict['4'] = {'rel':None,'abs':None}
        self.iodict['5'] = {}
        self.iodict['6'] = {'name':''}
        
        for key,valdict in iodict.items():
            if key in self.iodict.keys():
                for param in valdict.keys():
                    if param in self.iodict[key].keys():
                        self.iodict[key][param] = valdict[param]
            
        if rgrid is None:
            rmin = model.grid['min'][0]
            rmax = model.grid['max'][0]
            Nr = model.grid['N'][0]
            self.nr = int(Nr/2)
            self.rgrid = np.logspace(np.log10(rmin),np.log10(rmax), self.nr+1)
        else:
            self.rgrid = rgrid
            self.nr = len(rgrid)
            
        if zgrid is None:
            env_zrmax = 1./np.tan(np.radians(model.env['theta_min']*1.25))
            zrmin = np.tan(np.diff(model.theta)[-1]) #smallest z/r at midplane
            self.nz = 50
            self.zgrid = np.logspace(np.log10(zrmin),np.log10(env_zrmax),self.nz+1)
        else:
            self.zgrid = zgrid
            self.nz = len(zgrid)
           
        if 'gas' not in output.T.keys():
            output.calc_Tgas()
            
        self.data = {}
        self.inpfiles = {}
        
        self.set_1()
        self.set_2()
        self.set_3()
        self.set_4()
        self.set_5()
        self.set_6()
        
        try:
            os.mkdir(self.rundir)
        except:
            print('run directory exists - will overwrite current model if you write to it!')
        
        try:
            os.mkdir(self.chemdir)
        except:
            print('run directory exists !')
    
    
    def make_rz_H(self): 
        """ generates the grid for the chemistry
        2D, cylindrical, logarithmic in r,z
        
        Returns: ndarrays for R, Z cell-centers in two dimensions
        
        """
        model = self.input.m
        
        #turn faces into centers
        zc_norm = 0.5*(self.zgrid[1:] + self.zgrid[:-1])
        rc_pts = (self.rgrid[1:] + self.rgrid[:-1])*0.5
        
        R,Z = np.meshgrid(rc_pts,zc_norm)
        Z *= R #multiply by z/R to get z in au
        
        return R,Z
  
     #def make_quadrant(self,quantity_2d,fill_value=0,order='F',smooth=True): 
#         """ interpolates a 2D ndarray of data on the radmc3d grid onto the 2D grid for chemistry
#         Parameters:
#         ------------
#         quantity_2d: ndarray, 2d, data to be interpolated onto the chemical grid
        
#         fill_value: optional, default 0, value to put where chemistry grid extends past model grid
        
#         order: str, 'F' or 'C' for input arrays that are fortran indexed or c order indexed
        
#         smooth: boolean, if True, gaussian filter will smooth the data before interpolation (can be useful for noisy data)
        
#         Returns: ndarray of interpolated quantity on the chemical grid
#         """
        
#         model = self.input.m
#         r_cyl,z_cyl = model.make_rz()
#         r_new, z_new = self.make_rz_H()
        
#         if smooth == True:
#             quantity_2d = gaussian_filter(quantity_2d, sigma=[1,4])
#         quantity_2d_interp = griddata((r_cyl[:,:,0].flatten(),z_cyl[:,:,0].flatten()), quantity_2d.ravel(order=order), (r_new,z_new),fill_value=fill_value,method='linear',rescale=True)
        
#         return quantity_2d_interp
    
    def make_quadrant(self,quantity):
        from scipy import stats
        """ bins data on the radmc3d grid onto the 2D grid for chemistry

        interpolates any nan values for each zone vertically 
        * should be cautious about regions outside the original domain*

        Parameters:
        ------------
        quantity: ndarray, data to be binnedd onto the chemical grid, can be 3d or 2d

        Returns: 2d ndarray of interpolated quantity on the chemical grid
        """

        model = self.input.m
        r_cyl,z_cyl = model.make_rz()
        if len(np.shape(r_cyl)) > len(np.shape(quantity)):
            r_cyl = r_cyl[:,:,0]
            z_cyl = z_cyl[:,:,0]
        r_new, z_new = self.make_rz_H()

        if np.shape(r_cyl) != np.shape(quantity):
            quantity = quantity.swapaxes(0,1)

        r_old = r_cyl.flatten()
        zr_old = (z_cyl/r_cyl).flatten()
        q_old = quantity.flatten()

        ## bins in z/r ##
        quantity_interp = stats.binned_statistic_2d(r_old, zr_old, q_old, bins=[self.rgrid,self.zgrid])[0].T
        for i in range(self.nr):
            zone = quantity_interp[:,i]
            if len(np.isnan(zone)) >= 1:
                zcoord = z_new[:,i]
                new_zone = np.interp(zcoord,zcoord[~np.isnan(zone)],zone[~np.isnan(zone)])
                quantity_interp[:,i] = new_zone

        return quantity_interp

    
    def regrid_disk(self):
        """ regrids output data from the radmc3d grid onto a new grid for chemistry
        
        saves re-gridded values to self.data 
        
        Returns: 2-d ndarrays for density (gas, two dust populations) and temperature (gas and dust)
        """
        model = self.input.m
        
#         rho_g = self.input.calc_rho2D('gas')
#         rho_d1 = self.input.calc_rho2D('dust1')
#         rho_d2 = self.input.calc_rho2D('dust1')
        
#         T_d = self.input.calc_T2D('dust')
#         T_g = self.input.calc_T2D('gas')
        rho_g = self.input.rho['gas']
        rho_d1 = self.input.rho['dust1']
        rho_d2 = self.input.rho['dust2']
        
        T_d = self.input.T['dust']
        T_g = self.input.T['gas']
        
#         rhod1_2d = self.make_quadrant(rho_d1, fill_value = model.env['rho_amb']*model.env['d2g'],order='F',smooth=False)
#         rhod2_2d = self.make_quadrant(rho_d2, fill_value = 0,order='F',smooth=False)
#         rhog_2d = self.make_quadrant(rho_g, fill_value = model.env['rho_amb'], order='F',smooth=False)
        rhod1_2d = self.make_quadrant(rho_d1)
        rhod2_2d = self.make_quadrant(rho_d2)
        rhog_2d = self.make_quadrant(rho_g)
        
        Td_2d = self.make_quadrant(T_d)
        Tg_2d = self.make_quadrant(T_g)
        r,z = self.make_rz_H()
#         T_amb = model.T(np.sqrt(r**2 + z**2)) #ambient temperature from power-law temp profile
        
#         #temperature gets smoothed before it goes in there as noise from radmc can make interpolation wonky
#         # parts outside model get filled with ambient temperature 
#         Td_2d = self.make_quadrant(T_d,fill_value = 0,order='F')
#         Td_2d[Td_2d < 1] = T_amb[Td_2d < 1]
#         Td_2d = Td_2d
        
#         Tg_2d = self.make_quadrant(T_g,fill_value = 0,order='F')
#         Tg_2d[Tg_2d < 1] = T_amb[Tg_2d < 1]
#         Tg_2d = Tg_2d
        
        
        self.data['r'] = r
        self.data['z'] = z
        
        self.data['Tg'] = Tg_2d
        self.data['Td'] = Td_2d
        
        self.data['rhog'] = rhog_2d
        self.data['rhod'] = rhod1_2d + rhod2_2d
        self.data['rhod1'] = rhod1_2d
        self.data['rhod2'] = rhod2_2d
        
        return rhog_2d, rhod1_2d, rhod2_2d, Td_2d, Tg_2d 
        
        
    def write_out(self,outname='model',ndust=2):
        """ writes the model.out file 
        Parameters:
        -----------
        outname: the name of the .out file, to be saved in the chem environs directory
        
        ndust: int, default 2, number of dust populations to include
        
        Returns: str, name of the .out file 
        """
        model = self.input.m
        if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
        r,z = self.make_rz_H()
        header1 = 'R(AU)    z(AU)   rhogas(g/cm3)  rhodust(g/cm^3)    T(K)   [    4 types of dust]      f_H\n'
        header2 = '-ForChem+Plots--------------------------------------------------------\n'

        rhog_2d, rhod1_2d, rhod2_2d, Td_2d, Tg_2d  = self.regrid_disk()
        rhod_2d = rhod1_2d + rhod2_2d
        
        arrays_to_write = [r,z,rhog_2d, rhod_2d, Td_2d]
        for i in range(len(arrays_to_write)):
            array = arrays_to_write[i]
            arrays_to_write[i] = array.ravel(order='F')
            
        with open(self.chemdir + outname+'.out',"w") as outfile:
            outfile.write(header1)
            outfile.write(header2)
            Nr = len(r.flatten())
            for n in np.arange(Nr):
                outfile.write('%10.7f %10.7f   %5.2e   %5.2e %7.1f %12.7f %12.7f %12.7f %12.7f   %7.2e\n'%
                              ( arrays_to_write[0][n], arrays_to_write[1][n],arrays_to_write[2][n], 
                               arrays_to_write[3][n] ,arrays_to_write[4][n],1.0,0.0,0.0,0.0,1.0))
        return outname + '.out'
        
        
    def read_out(self,outname='model'):
        """ reads the selected .out file
        
        fills the self.data structure with read-in values
        
        """
        if os.getcwd() != self.chemdir:
            os.chdir(self.chemdir)
        if os.path.exists(outname+'.out') == False:
            self.write_out(outname)
        data = np.loadtxt(outname+'.out',skiprows=2,usecols=(0,1,2,3,4))
        r = data[:,0]
        z = data[:,1]
        shape = (50,len(np.unique(r)))
        rhog = data[:,2].reshape(shape)
        rhod = data[:,3].reshape(shape)
        Td = data[:,4].reshape(shape)
        self.data['r'] = r.reshape(shape)
        self.data['z'] = z.reshape(shape)
        self.data['rhog'] = rhog
        self.data['rhod'] = rhod
        self.data['Td'] = Td
        
        
    def calc_zeta(self,cr_model=''):
        """ calculates the cosmic ray ionization rate using cosmic-ray model rates
        fills the 'zeta', 'zcm' columns of self.data structure

        """
        if self.data == {}:
            self.read_out()
        if cr_model == '':
            cr_model = self.input.m.rad['cr_model']

        def zetaeval(ncol,model=cr_model):
            labels = ['m02','w98','ssm','ssx','ttm','ttx']
            #Rates from Cleeves et al. 2014
            #zeta_powerlaw alpha zeta_exponential Sigma
            #allrts = np.array([[6.8e-16,3.7e-18,0.423,210.0],
                #[2.0e-17,9.4e-19,0.021,260.0],
                #[1.1e-18,3.0e-19,-0.00,260.0],
                #[1.6e-19,7.5e-20,-0.01,250.0],
                #[7.0e-21,4e-21,-0.01,290.0],
                #[1.1e-22,3e-23,-0.02,490.0]])

            allrts = np.array([[6.8e-16,3.7e-18,0.423,210.0],
                [2.0e-17,9.4e-19,0.021,260.0],
                [1.3e-18,3.0e-18,-0.00,190.0],
                [2.0e-19,8.0e-19,-0.01,230.0],
                [1.0e-20,2.0e-19,-0.03,270.0],
                [3.0e-22,2.0e-19,-0.03,270.0]])

            modin = labels.index(model)
            zp = allrts[modin,0]
            ze = allrts[modin,1]
            al = allrts[modin,2]
            co = allrts[modin,3]
            mumult = 2.36
            piv = 1e20
    #             F = zp*ze*mumult/(ze*mumult*(ncol/piv)**al+zp*(np.exp(ncol*2.0*mh/co)-1)) #C13 eq5
            F = (zp*ze)/((ze*(ncol/piv)**al) + zp*(np.exp(ncol*mh*mumult/co) - 1.))

            return F

        Ncol_external = 1e10 #column before disk surface

        rhog = self.data['rhog']/(mu*mh)
        axis = 0 #zaxis
        zc = self.data['z']
        dz = np.abs(np.gradient(zc,axis=axis)*AU) #delta z/z 

        #sandwich both disk halves 
        rhog_los = np.concatenate((np.flip(rhog,axis=axis),rhog),axis=axis)
        dz_los = np.concatenate((np.flip(dz,axis=axis),dz),axis=axis)

        #integrate column from above and below
        Ncol_zneg = np.cumsum(rhog_los*dz_los,axis=axis) + Ncol_external
        Ncol_zpos = np.flip(np.cumsum(np.flip(rhog_los*dz_los,axis=axis),axis=axis),axis=axis) + Ncol_external

        zeta = 0.5*zetaeval(Ncol_zpos,model=cr_model) + 0.5*zetaeval(Ncol_zneg,model=cr_model)

        #take the top half 
        self.data['zeta'] = zeta[self.nz:,:]
        self.data['zcm'] = np.flip(np.cumsum(dz_los[:self.nz,:],axis=axis),axis=axis)
        #self.data['Nrz'] = np.ones_like(self.data['z'])*Ncol_external

        
    def calc_dustfrac(self):
        """ calculates the dustfrac and ngr column of self.data structure
        
        to go with updated to chemical code to adjust for multiple size distributions
        and an accurate local dust density
       
        """
        if self.data == {}:
            self.read_out()
        model = self.input.m
        #cross-section assumed by chemistry
        a_chem = 0.1/1e4 #microns to cm
        sa_chem = 4.*pi*a_chem**2 
        
        ndust = len(model.dust['amin'])
        sa_tot = np.zeros_like(self.data['rhod'])
        n_tot = np.zeros_like(self.data['rhod'])
        # Add up the total area and total number across populations
        for d in range(ndust):
            key = 'rhod'+str(d+1)
            p = model.dust['apow'][d]
            amin = model.dust['amin'][d]*1e-4
            amax = model.dust['amax'][d]*1e-4
            rhod = self.data[key]
            # n0: normalization factor --> dn/da = n0 * a**(p) 
            # solved for by equating total mass with density in the cell
            n0 = (rhod/model.dust['rho_si'])*(3/(4*pi))*(4-p)/(amax**(4-p) - amin**(4-p))
            sa_tot += 4*pi*n0*(amax**(3-p) - amin**(3-p)) / (3-p)
            n_tot += n0*(amax**(1-p) - amin**(1-p))/ (1-p)
        self.data['dustfrac'] = ((sa_tot)/(n_tot))/sa_chem
        self.data['ngrain'] = n_tot
        
    def regrid_radiation(self,field='uv'):
        """ regrids the radiation field to generate radiation input files to the chemical code
        radiation is interpolated onto a new spatial grid and a new wavelength grid at the same time
        
        Parameters:
        -----------
        field: str, 'uv' or 'xray', which radiation field to regrid.
        
        Returns: field for writing to a radiation file for the chemistry
        """
        from scipy import interpolate
        model = self.input.m

        wav_rt, freq, Jnu = self.input.read_Jnu()
        
        if field == 'uv':
            field_index = (wav_rt >= uv_min) & (wav_rt <= uv_max)
            wav_chem = np.arange(930,2010,10) #in angstroms
            wav_mu = wav_chem * 1e-4 #angstroms to microns
            field_index_comp = (wav_mu >= uv_min) & (wav_mu <= uv_max)
            wav_conv = 1e-8* c/(wav_mu*1e-4)**2 #convert /Hz to /angstrom dnu/dwav
        
        if field == 'xray':
            wav_chem = np.arange(1,21,1)[::-1] #in keV
            wav_mu = (1e4*h*c)/(wav_chem*keV) #in microns
            field_index = (wav_rt >= xray_min) & (wav_rt <= xray_max)
            field_index_comp = (wav_mu >= xray_min) & (wav_mu <= xray_max)
            wav_conv = (keV/h) #convert /Hz to /keV
            
        freq_field = freq[field_index]
        wav_field = wav_rt[field_index]
        
        #fnu_out_grid = np.average(Jnu[:,:,:,field_index],axis=2)*4.*pi
        fnu_out_grid = Jnu[:,:,:,field_index]*4.*pi

        #input spectra
        fnu0_1d_hires = calc_input_spectrum(model,wav=wav_mu)
        fnu0_1d = calc_input_spectrum(model,wav=wav_field)
        fnu0_1d_interp = 10**np.interp(np.log10(wav_mu),np.log10(wav_field),np.log10(fnu0_1d))
        
        f0_tot = np.trapz(fnu0_1d_interp, x = c/wav_mu * 1e4)*-1.
        f0_tot_hires = np.trapz(fnu0_1d_hires, x = c/wav_mu * 1e4)*-1.

        # error in the total flux due to lower resolution of mcmono wavelengths
        err_sum = f0_tot_hires/f0_tot - 1.

        #chemical grid spatial information
        R2_chem = np.expand_dims((self.data['r']*AU)**2 + (self.data['z']*AU)**2,axis=-1)

        #input spectrum on chemical spatial grid
        fnu0_chem_grid = np.expand_dims(fnu0_1d,axis=(0,1))*(pc**2 / R2_chem)
        fnu0_chem_grid_interp = np.expand_dims(fnu0_1d_interp,axis=(0,1))*(pc**2 / R2_chem)
        fnu0_chem_grid_hires = np.expand_dims(fnu0_1d_hires,axis=(0,1))* (pc**2 / R2_chem)

        #adding in the external ISRF to the initial field
        if model.rad['G0'] > 0 and field=='uv':
            dum, isrf_1d = calc_ISRF(model,wav=wav_field)
            dum, isrf_1d_hires = calc_ISRF(model,wav=wav_mu)

            Rg = np.amax(np.sqrt((self.data['r']*AU)**2 + (self.data['z']*AU)**2))
            Rgg = np.amax(np.sqrt((self.data['r']*AU)**2))

            isrf0_chem_grid = 4*pi*np.expand_dims(isrf_1d,axis=(0,1))*(Rgg**2)/(Rg - np.sqrt(R2_chem))**2
            isrf0_chem_grid_hires = 4*pi*np.expand_dims(isrf_1d_hires,axis=(0,1))*(Rgg**2)/(Rg - np.sqrt(R2_chem))**2

            fnu0_chem_grid += isrf0_chem_grid
            fnu0_chem_grid_hires += isrf0_chem_grid_hires

        #interpolating fluxes onto chemistry (spatial) grid
        for j in np.arange(len(freq_field)):
            fnu_chem_grid_temp = np.expand_dims(self.make_quadrant(fnu_out_grid[:,:,:,j]),axis=-1)
            fnu_chem_grid_temp = np.minimum(fnu_chem_grid_temp, np.expand_dims(fnu0_chem_grid[:,:,j],axis=-1)) 
            if j < 1 :
                fnu_chem_grid = fnu_chem_grid_temp.copy()
            else:
                fnu_chem_grid = np.dstack((fnu_chem_grid,fnu_chem_grid_temp))

        #function to interpolate between attenuated fluxes for wavelengths
        coeff = interpolate.interp1d(x=wav_field, y=(fnu_chem_grid/fnu0_chem_grid),axis=-1,fill_value='extrapolate',kind='slinear')

        #coeff_chem_grid = coeff(wav_field)
        coeff_chem_grid_hires = coeff(wav_mu)

        #f_tot = np.trapz(coeff_chem_grid*fnu0_chem_grid, x = c/wav_field * 1e4,axis=-1)*-1
        f_tot = np.trapz(coeff_chem_grid_hires*fnu0_chem_grid_interp, x = c/wav_mu * 1e4,axis=-1)*-1
        f_tot_hires = np.trapz((coeff_chem_grid_hires*fnu0_chem_grid_hires), x = c/wav_mu * 1e4,axis=-1)*-1
        
        #correction factor to ensure total flux remains the same for interpolated spectrum
        #usually this is on the order of a few percent 
        norm_fac = np.expand_dims(f_tot*(1.+err_sum)/f_tot_hires, axis=-1)
        norm_fac[np.isnan(norm_fac)] = 1.
        
        photon_e = 1./(h*c/(wav_mu*1e-4))
        
        spec_final = coeff_chem_grid_hires*fnu0_chem_grid_hires*norm_fac*wav_conv*photon_e
        spec_final[np.isnan(spec_final)] = 0.0
        return spec_final
    

    
    def set_1(self):
        """ calculates necessary columns for 1environ.inp
        
        + regrids temperature + density in model
        + calculates scaling factors for dust density populations
        + calculates the cosmic ray attentuation
        
        """
        self.regrid_disk()
        self.calc_dustfrac()
        self.calc_zeta()
    
    def set_2(self,tkwargs={}):
        """ sets values for 2times.inp from template with the option 
        to set new values
        
        Parameters:
        -----------
        tkwargs: optional, dict
        
        keys: 'first','last','total'
        the value of the first, last, and total timesteps for integration
        
        """
        template_path = self.input.m.parent_dir+'chemistry/templates/'
        template_file = template_path + '2times.inp'
        file_contents = []
        
        if tkwargs == {}:
            tkwargs = self.iodict['2']
            print('2times.inp: default values from template file')
        else:
            for key in self.iodict['2'].keys():
                if key in tkwargs.keys():
                    self.iodict['2'][key] = tkwargs[key]
                    print('2times.inp: {} set to {}'.format(key,tkwargs[key]))
                    
        with open(template_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    file_contents.append(line)
                else:
                    val = line.split('#')[0].strip()
                    comment = line.split('#')[-1].strip()
                    if comment.startswith('last') and tkwargs['last'] is not None:
                        val = '{:9.3E}'.format(tkwargs['last']).replace('E','D')
                    if comment.startswith('first') and tkwargs['first'] is not None:
                        val = '{:9.3E}'.format(tkwargs['first']).replace('E','D')
                    if comment.startswith('total') and tkwargs['total'] is not None:
                        val = str(tkwargs['total'])
                    new_line = val + '\t #' + comment + '\n'
                    file_contents.append(new_line)
        self.inpfiles['2times'] = file_contents
        
    def set_3(self,fkwargs={}):
        """
        updates values for 3abunds.inp from a template file
        with new values from an input dictionary
        * note: this will only replace values for existing species in the template file
        (including those commented out) but will not write new species
        
        Parameters:
        ------------
        fkwargs, dict with keys:
        
        name: str, optional: file suffix for specific template file 
        
        new_values: dict, optional
                    keys: grain name, case sensitive
                    values: new initial abundance to set
        """
        
        if fkwargs == {}:
            fkwargs = self.iodict['3']
        else:
            for key in self.iodict['3'].keys():
                if key in fkwargs.keys():
                    self.iodict['3'][key] = fkwargs[key]
                    
                    
        template_path = '/'.join([self.parentchemdir,'templates/'])
        template_file = template_path + '3abunds' + fkwargs['name'] + '.inp'
        print('3abunds.inp: using values from template file' + template_file)
        file_contents = []
        
        new_values = fkwargs['newvals']
        
        with open(template_file, 'r') as f:
            for line in f:
                if line.startswith('#') == True:
                    old_species  = line.strip('#').split(' ')[0].strip()
                    if old_species in new_values.keys():
                        print('3abunds.inp: updating ' + old_species )
                        new_val = '{:9.3E}'.format(new_values[old_species]).replace('E','D')
                        new_line = '{:<11} = {}\n'.format(old_species,new_val)
                    else:
                        pass
                else: 
                    old_species = line.split(' ')[0].strip()
                    if old_species in new_values.keys(): #if resetting, update the value
                        print('3abunds.inp: updating ' + old_species )
                        new_val = '{:9.3E}'.format(new_values[old_species]).replace('E','D')
                        new_line = '{:<11} {}\n'.format(old_species,new_val)
                    else:
                        new_line = line
                    file_contents.append(new_line)
        self.inpfiles['3abunds'+fkwargs['name']] = file_contents
    
    def set_4(self, tolkwargs = {}):
        """ sets the values for the 4toleran.inp file from template with the option to set new ones
        Parameters:
        ----------
        tolkwargs: dict, optional 
        keys: 'rel','abs': relative and absolute tolerance levels
        
        """
        template_path = '/'.join([self.parentchemdir,'templates/'])
        template_file = template_path + '4toleran.inp'
        file_contents = []
        
        
        if tolkwargs == {}:
            tolkwargs = self.iodict['4']
            print('4toleran.inp: using default value from template file')
        else:
            for key in self.iodict['4'].keys():
                if key in tolkwargs.keys():
                    self.iodict['4'][key] = tolkwargs[key]
                    print('4toleran.inp: updating {} to {}'.format(key,tolkwargs[key]))
        
        with open(template_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    file_contents.append(line)
                else:
                    val = line.split('#')[0].strip()
                    comment = line.split('#')[-1].strip()
                    if comment.startswith('rel') and tolkwargs['rel'] is not None:
                        val = '{:9.3E}'.format(tolkwargs['rel']).replace('E','D')
                    if comment.startswith('abs') and tolkwargs['abs'] is not None:
                        val = '{:9.3E}'.format(tolkwargs['abs']).replace('E','D')
                    new_line = val + '\t #' + comment + '\n'
                    file_contents.append(new_line)
        self.inpfiles['4toleran'] = file_contents
    
    def set_5(self,fkwargs={}):
        """ sets the values for '5flags.inp' from template
        with the option to set new flags
        
        Parameters:
        -----------
        fkwargs: dict, optional 
        keys: flags in 5flags.inp
        values: 0 or 1 to turn flag on or off
        
        """
        
        if fkwargs == {}:
            fkwargs = self.iodict['5']
            print('5flags.inp: using default values from template file')
        else:
            self.iodict['5'] = fkwargs
            
        
        
        template_path = '/'.join([self.parentchemdir,'templates/'])
        template_file = template_path + '5flags.inp'
        file_contents = []
        with open(template_file, 'r') as f:
            for line in f:
                if line.startswith('#') == True:
                    file_contents.append(line)
                else: #this is a line that sets parameters
                    param, val = line.split('=')
                    param = param.strip()
                    if param in fkwargs.keys(): #if resetting, update the value
                        print('5flags.inp: updating'.foramt(param, fkwargs[param]))
                        new_line = '{:<15} = {}\n'.format(param,fkwargs[param])
                    else:
                        new_line = line
                    file_contents.append(new_line)
        self.inpfiles['5flags'] = file_contents

    def set_6(self,fkwargs={}):
        """ reads in a template 6grainbe.inp file
        Parameters:
        ------------
        name: str, optional, file suffix for template name
        
        """
        
        if fkwargs == {}:
            fkwargs = self.iodict['6']
        else:
            for key in self.iodict['6'].keys():
                if key in fkwargs.keys():
                    self.iodict['6'][key] = fkwargs[key]
                    
        template_path = '/'.join([self.parentchemdir,'templates/'])
        template_file = template_path + '6grainbe' + fkwargs['name'] + '.inp'
        print('6grainbe.inp: using values from template file' + template_file)
        
        file_contents = []
        with open(template_file, 'r') as f:
            for line in f:
                file_contents.append(line)
        self.inpfiles['6grainbe'] = file_contents
                

def contents_to_file(chem_disk,file_contents,fileprefix,filepath=None):
    """ helper function to write input file content to requisite file
    Parameters:
    ----------
    chem_disk: chem_disk object 
    file_contents: lines to write to fle
    file_prefix: name of file to be written
    """
    if filepath is None:
        filepath = chem_disk.rundir + fileprefix + '.inp'
    #filepath = chem_disk.input.m.parent_dir+'chemistry/' + fileprefix + '.inp'
    with open(filepath,'w') as f:
        for line in file_contents:
            f.write(line)

            

def write_0(chem_disk,fkwargs = {}):
    """ writes 0io input file 
    Parameters:
    ----------
    chem_disk: chem_disk object
    
    reactions: str, suffix of rreacs file to use
    
    species: str, suffix of rspecies file to use
    
    name: str, optional, suffix for 0io.inp file
    
    abund0: str, initial 2D abundance file name default: None
    
    """
    
    if fkwargs == {}:
        fkwargs = chem_disk.iodict['0']
    else:
        for key in chem_disk.iodict['0'].keys():
            if key in fkwargs.keys():
                chem_disk.iodict['0'][key] = fkwargs[key]
                    
    mdir = chem_disk.input.m.outdir.strip('/').split('out/')[-1]
    template_path = '/'.join([chem_disk.parentchemdir,'templates/'])
    template_reacs_file = 'rreacs_' + fkwargs['reactions'] + '.dat'
    template_species_file = 'rspecies_' + fkwargs['reactions'] + '.dat'
    
    file_error = 0
    
    with open(chem_disk.rundir + '0io' + fkwargs['name'] + '.inp', 'w') as f:
        f.write('# input & output files:\n')
        f.write('rspecies_{}.dat  \n'.format(fkwargs['species']))
        f.write('rreacs_{}.dat  \n'.format(fkwargs['reactions']))
        f.write('uv_photons_{}.dat  \n'.format(mdir))
        f.write('xray_photons_{}.dat \n'.format(mdir))
        f.write('None  \n') #ISRF is already in outputs from uv/xray field
        f.write('None  \n' ) #Radionucleide file
        f.write('{}    \n'.format(fkwargs['abund0'])) #default is None
    
    for template_file in [template_reacs_file, template_species_file]:
        print('0inp: using {}'.format(template_reacs_file))
        exit = os.system('cp {} {}'.format(template_path+template_file, chem_disk.rundir))
        if exit != 0:
            print('0inp: ERROR copying {} to run directory'.format(template_file))
            file_error += 1
    
    #if you specified an abund0 file it will look for it in the run directory, otherwise, it will look for it in the templates directory
    if fkwargs['abund0'] is not None:
        if os.path.exists(chem_disk.rundir + fkwargs['abund0']) == True:
            print('0inp: using {} in run directory'.format(fkwargs['abund0']))
        else:
            print('0inp: did not find {} in run directory: {}'.format(fkwargs['abund0'], chem_disk.rundir))
            if os.path.exists(template_path + fkwargs['abund0']) == True:
                exit = os.system('cp {} {}'.format(template_path+fkwargs['abund0'], chem_disk.rundir))
                if exit != 0:
                    print('0inp: ERROR copying {} to run directory'.format(fkwargs['abund0']))
                    file_error += 1
                else:
                    print('0inp: using {} in templates directory'.format(fkwargs['abund0']))
            else:
                print('0inp: ERROR did not find {} in either run or templates directory'.format(fkwargs['abund0']))
                file_error += 1
    return file_error
              
        
def write_1(chem_disk,index):
    """ writes a 1environ.inp file for a single zone
    
    Parameters:
    ----------
    chem_disk: chem_disk object
    
    index: int, zone number of the radius row to write
    
    """
    header = ['R','rho','ngr','Tgas','Tdust','zAU','zcm','ZetaCR','DustFrac']
    keys = ['r', 'rhog','ngrain','Tg','Td','z','zcm','zeta','dustfrac']
    vals = [chem_disk.data[key][:,index] for key in keys]
    towrite = dict(zip(keys,vals))
    file_error = 0
    nz = chem_disk.nz
    name = '.e1.{:.4f}'.format(towrite['r'][0])
    with open(chem_disk.chemdir + '1environ.inp' + name, 'w') as f:
        f.write('{:12} {:12} {:12} {:12} {:12} {:12} {:12} {:12} {:12}\n'.format(*header))
        f.write('1\n')
        f.write('{:}\n'.format(nz))
        for j in np.arange(nz):
            row = [towrite[key][::-1][j] for key in keys]
            if np.any(np.isnan(row)):
                print('1environ.inp{} ERROR: nans found'.format(name))
                file_error += 1
            f.write('{:11.5E} {:11.5E} {:11.5E} {:11.5E} {:11.5E} {:11.5E} {:11.5E} {:11.5E} {:11.5E}\n'.format(*row))
    return file_error
            
def write_environ(chem_disk):
    """ writes all the 1environ.inp files
    Parameters:
    -----------
    chem_disk: chem_disk object
    
    """
    fe = write_0(chem_disk)
    nr = chem_disk.nr
    for i in np.arange(nr):
        fe += write_1(chem_disk,i)
    return fe
            
    
def write_uv(chem_disk,filepath = None):
    """ writes the uv.dat file for the chemical code
    Parameters:
    -----------
    chem_disk: chem_disk object
    
    """
    wav_chem = np.arange(930,2010,10) #in angstroms
    wav_mu = wav_chem*1e-4 #in microns
    spectrum = chem_disk.regrid_radiation(field='uv') #z,r,wavelength
    run_name = chem_disk.input.m.outdir.strip('/').split('out/')[-1]
    header = 'Results from RT code with run: {}'.format(run_name)
    file_error = 0
    
    if filepath is None:
        filepath = chem_disk.rundir
    #with open(chem_disk.input.m.parent_dir+'chemistry/' + 'uv_photons_{}.dat'.format(run_name),'w') as f:
    with open(filepath + 'uv_photons_{}.dat'.format(run_name),'w') as f:
        f.write('{}\n'.format(header))
        for j in range(np.shape(spectrum)[1]):
            rheader = 'Radius(AU) {:9.4f} \n'.format(chem_disk.data['r'][0,j])
            zheader = 'z(AU)' + '{:10.4f}'*np.shape(spectrum)[0]
            zheader_vals = chem_disk.data['z'][:,j]
            zhead = zheader.format(*zheader_vals) + '\n'
            wavhead = 'Wavelength A               Photons/cm2/s/A \n'
            f.write(rheader)
            f.write(zhead)
            f.write(wavhead)
            row = zip(wav_chem, [spectrum[:,j,-i] for i in range(np.shape(spectrum)[2])])
            for r in row:
                if np.any(np.isnan(r[1])):
                    print('uv file ERROR: nans found'.format(name))
                    file_error += 1
                f.write('{:7.2f}'.format(r[0]) + ('{:10.2e}'*np.shape(spectrum)[0]).format(*r[1]) + '\n')
    return file_error

def write_xray(chem_disk,filepath=None):
    """ writes the xray.dat file for the chemical code
    Parameters:
    -----------
    chem_disk: chem_disk object
    
    """
    wav_chem = np.arange(1,21,1) #in keV
    wav_mu = (1e4*h*c)/(wav_chem*keV) #in microns
    spectrum = chem_disk.regrid_radiation(field='xray') #order is z, r, wavelength
    run_name = chem_disk.input.m.outdir.strip('/').split('out/')[-1]
    header = 'Results from RT code with run: {}'.format(run_name)
    file_error = 0
    if filepath is None:
        filepath = chem_disk.rundir
    #with open(chem_disk.input.m.parent_dir+'chemistry/' + 'xray_photons_{}.dat'.format(run_name),'w') as f:
    with open(filepath + 'xray_photons_{}.dat'.format(run_name),'w') as f:
        f.write('{}\n'.format(header))
        for j in range(np.shape(spectrum)[1]):
            rheader = 'Radius(AU) {:9.4f}\n'.format(chem_disk.data['r'][0,j])
            zheader = 'z(AU)' + '{:10.4f}'*np.shape(spectrum)[0]
            zheader_vals = chem_disk.data['z'][:,j]
            zhead = zheader.format(*zheader_vals) + '\n'
            wavhead = 'Energy(keV)             Photons/cm2/s/keV \n'
            f.write(rheader)
            f.write(zhead)
            f.write(wavhead)
            row = zip(wav_chem, [spectrum[:,j,-i] for i in range(np.shape(spectrum)[2])])
            for r in row:
                if np.any(np.isnan(r[1])):
                    print('xray file ERROR: nans found'.format(name))
                    file_error += 1
                f.write('{:7.2f}'.format(r[0]) + ('{:10.2e}'*np.shape(spectrum)[0]).format(*r[1]) + '\n')
    return file_error
                
def write_chem_inputs(chem_disk,filepath=None):
    """ writes all the chem input files at once
    Parameters:
    -----------
    chem_disk: chem_disk object
    
    """
    # set default file content
    model = chem_disk.input.m
    file_error = 0
    file_error += write_environ(chem_disk)
    for filename, contents in chem_disk.inpfiles.items():
        contents_to_file(chem_disk,contents,filename,filepath=filepath)
        
    file_error += write_uv(chem_disk,filepath=filepath)
    if model.rad['xray'] == True:
        file_error += write_xray(chem_disk,filepath=filepath)
    
    if filepath is None:
        filepath = chem_disk.rundir
    if os.path.exists(filepath + 'disk_chemistry') != True:
        exit = os.system('cp {} {}'.format(chem_disk.parentchemdir + '/src/disk_chemistry', filepath))
        if exit != 0:
            print('chemistry: ERROR copying disk_chemistry to:'.format(filepath))
            file_error += 1
    if file_error > 0:
        print('WARNING: There are at least {} errors in the input that will cause problems if you try to run as is'.format(file_error))
    else:
        print('SUCCESS: ready to run the chemistry')
                

def plot_out(chem_disk,prop='rhog',log=True,method=contourf,**pk): 
    """ plots basic .out file data in r,z slice
    
    Parameters:
    -----------
    chem_disk: chem_disk object
    prop: str, key of property to plot
    log: default True, plot the logarithm of the property
    
    method: plotting function to use, default is contourf
    
    pk: dict, optional, plotting keywords
    
    """
    if chem_disk.data == {}:
        chem_disk.regrid_disk()
    data = chem_disk.data
    if log == True:
        method(data['r'],data['z'],np.log10(data[prop]),**pk)
    else:
        method(data['r'],data['z'],(data[prop]),**pk)
        
def plot_prechem(chem_disk,rlim=(400,400)): 
    """ plots model properties on the old and new grid for comparison
    
     Parameters:
    -----------
    chem_disk: chem_disk object
    
    rlim: tuple, optional, xlim and ylim for the plot
    
    """
    model = chem_disk.input.m
    X,Z = model.make_rz()
    X = X[:,:,0]
    Z = Z[:,:,0]
    f,ax = subplots(2,4,constrained_layout=True)
    f.set_size_inches(8,6)
    
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
        
    c = ax[0,0].contourf(X, Z, np.log10(chem_disk.input.calc_rho2D('gas').T),levels=np.arange(model.rhomin,model.rhomax,1),extend='both')
   
    ax[0,0].set_ylabel('z [au]')
    ax[0,0].set_title(r'Original $\rho_g$',fontsize=11)
    cb = colorbar(c,ax=ax[:,0],location='bottom',aspect=10,ticks=np.array(c.levels[::2]).astype(int))
    cb.set_label(r'$\rho$ [$\mathrm{g/cm^{-3}}$]')
    
    sca(ax[1,0])
    plot_out(chem_disk,prop='rhog',cmap='viridis',levels=np.arange(model.rhomin,model.rhomax,1),log=True,extend='both')
    ax[1,0].set_ylabel('z [au]')
    ax[1,0].set_title(r'Regridded $\rho_g$ ',fontsize=11)
    
    
    c = ax[0,1].contourf(X, Z, np.log10(chem_disk.input.calc_rho2D('dust').T),levels=np.arange(model.rhomin,model.rhomax,1),extend='both')
    ax[0,1].set_title(r'Original $\rho_d$',fontsize=11)
    cb = colorbar(c,ax=ax[:,1],location='bottom',aspect=10,ticks=np.array(c.levels[::2]).astype(int))
    cb.set_label(r'$\rho$ [$\mathrm{g/cm^{-3}}$]')
    
    sca(ax[1,1])
    plot_out(chem_disk,prop='rhod',cmap='viridis',levels=np.arange(model.rhomin,model.rhomax,1),log=True,extend='both')
    ax[1,1].set_title(r'Regridded $\rho_d$',fontsize=11)
    
    
    
    #from scipy.ndimage import gaussian_filter
    smooth_T = chem_disk.input.calc_T2D('dust').T
    c2 = ax[0,2].contourf(X, Z, smooth_T,levels=np.linspace(5,125,19),cmap='twilight_shifted',extend='both')


    ax[0,2].set_title(r'Original $T_d$',fontsize=11)
    cb2 = colorbar(c2,ax=ax[:,2],location='bottom',ticks=np.array(c2.levels[::3]).astype(int),aspect=10)
    cb2.set_label(r'$T \ \mathrm{[K]}$')
    
    sca(ax[1,2])
    ax[1,2].set_title(r'Regridded $T_d$',fontsize=11)
    plot_out(chem_disk,prop='Td',cmap='twilight_shifted',levels=np.linspace(5,125,19),log=False,extend='both')
    
    #from scipy.ndimage import gaussian_filter
    smooth_T = chem_disk.input.calc_T2D('gas').T
    c2 = ax[0,3].contourf(X, Z, smooth_T,levels=np.linspace(5,125,19),cmap='twilight_shifted',extend='both')
    ax[0,3].set_title(r'Original $T_g$',fontsize=11)
    cb2 = colorbar(c2,ax=ax[:,3],location='bottom',ticks=np.array(c2.levels[::3]).astype(int),aspect=10)
    cb2.set_label(r'$T \ \mathrm{[K]}$')
    
    sca(ax[1,3])
    ax[1,3].set_title(r'Regridded $T_g$',fontsize=11)
    plot_out(chem_disk,prop='Tg',cmap='twilight_shifted',levels=np.linspace(5,125,19),log=False,extend='both')
    
    ax[1,0].set_xlabel('r [au]')
    ax[1,1].set_xlabel('r [au]')
    ax[1,2].set_xlabel('r [au]')
    ax[1,3].set_xlabel('r [au]')
    
    for a in ax[0,:]:
        a.set_xlim(0,rlim[0])
        a.set_ylim(0,rlim[1])
        a.xaxis.set_major_locator(MultipleLocator(int(rlim[0]/5)))
        a.yaxis.set_major_locator(MultipleLocator(int(rlim[1]/5)))
    
    for a in ax[1,:]:
        a.set_xlim(0,rlim[0])
        a.set_ylim(0,rlim[1])
        a.xaxis.set_major_locator(MultipleLocator(int(rlim[0]/5)))
        a.yaxis.set_major_locator(MultipleLocator(int(rlim[1]/5)))
    f.savefig(chem_disk.chemdir+'regrid.pdf')
        