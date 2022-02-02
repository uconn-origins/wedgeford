from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from models.units import *

parent_dir = os.getcwd()
models_dir = parent_dir+'/models/'

class model:
    def __init__(self,stellar_params, disk_params, envelope_params,grid_params, dust_params, RT_params, outdir = '/m1_test/'):
        self.star = stellar_params
        self.disk = disk_params
        self.env = envelope_params
        self.grid = grid_params
        self.dust = dust_params
        self.outdir = parent_dir+outdir
        self.coords = [0,0,0]
        self.rhomin = np.log10(self.env['rho_amb'])
        self.rhomax = -11
        
        try:
            os.mkdir(self.outdir)
        except:
            print('directory exists - will overwrite current model if you write to it!')
        for dim in [0,1,2]:
            N,gmin,gmax,spacing = [self.grid['N'][dim], self.grid['min'][dim], self.grid['max'][dim], self.grid['spacing'][dim]]
            if spacing == 'log':
                self.coords[dim] = np.logspace(np.log10(gmin),np.log10(gmax),N+1)
            elif spacing == 'lin':
                self.coords[dim] = np.linspace(gmin,gmax,N+1)
            else:
                print('Defaulting to linear spacing')
                self.coords[dim] = np.linspace(gmin,gmax,N+1) # these are cell faces
            dim += 1
        self.r = (self.coords[0][1:] + self.coords[0][:-1]) / 2. # these are cell centers
        self.theta = (self.coords[1][1:] + self.coords[1][:-1]) / 2.
        self.phi = (self.coords[2][1:] + self.coords[2][:-1]) / 2.
        
    def write_star(self):
        params = rpy.params.radmc3dPar()
        if os.path.exists(self.outdir+ 'problem_params.inp') != True:
            params.loadDefaults() #loads default parameter dictionary
        else:
            params.readPar()
        if os.getcwd() != models_dir:
            os.chdir(models_dir)
        rad_data = rpy.radsources.radmc3dRadSources()
        rad_data.readStarsinp('stars_og.inp') #reads in a template stars.inp file
        params.ppar['mstar'] = [self.star['Ms']*Msun] # sets the paramaters to the ones passed to the model
        params.ppar['rstar']= [self.star['Rs']*Rsun]
        params.ppar['tstar'] = [self.star['Ts']*1.0]
        rad_data.pstar = params.ppar['pstar']
        rad_data.tstar = params.ppar['tstar']
        if os.getcwd() != self.outdir:
            os.chdir(self.outdir)
        rad_data.writeStarsinp(ppar=params.ppar) #writes a stars.inp in the new directory
        os.chdir(parent_dir)
        
    def calc_ISRF(self,lam=None,gnorm=1,write=True):
        from scipy.interpolate import interp1d
        if os.path.exists('wavelength_micron.inp') != True:
                self.write_wavelength()
        if os.getcwd() != self.outdir:
            os.chdir(self.outdir)
        grid_data = rpy.reggrid.radmc3dGrid()
        grid_data.readWavelengthGrid()
        #goal is to write the ISRF onto the same wavelength grid as the stellar source
        if lam == None:
            lam = grid_data.wav
            nu = grid_data.freq
        else:
            nu = c/(lam*1e-4)
        os.chdir(models_dir)
        ### works right now with the current format of the ISRF file
        ### if you want a different ISRF as the basis, should rewrite for other scalings
        fname='ISRF.csv'
        ISRF = np.loadtxt('ISRF.csv',skiprows=1,delimiter=',')
        lami = ISRF[:,0]#micron
        flam = ISRF[:,1]/(4*pi) #ergs/cm2/s/micron/str
        fnu_ = flam*(lami*(lami*1e-4))/c #conversion to ergs/cms2/s/Hz/str
        Fnu = interp1d(lami, fnu_,fill_value='extrapolate')
        fnu = np.clip(Fnu(lam),a_min=np.amin(fnu_)*1e-2,a_max=None)
        isrf_index = (lam > 0.0912) & (lam < 2.4) #comparing G0
        Ftot = np.trapz(fnu[isrf_index][::-1],x=nu[isrf_index][::-1])
        norm = gnorm*G0/Ftot
        fnu *= norm #scales ISRF by how many G0s you want.
        if write == True:
            with open(self.outdir+'external_source.inp','w') as f:
                f.write('2 \n')
                f.write(str(len(lam))+ '\n')
                f.write('\n')
                lam.tofile(f, sep='\n', format="%13.6e")
                f.write('\n')
                f.write('\n')
                fnu.tofile(f, sep='\n', format="%13.6e")
        os.chdir(parent_dir)
        return lam,fnu
        
    def T(self,R):
        return self.star['Ts']*np.sqrt((self.star['Rs']*Rsun)/(R*AU))
    
    def Rsub(self): #dust sublimation radius at T = 1500 K
        Tsub = 1500.
        T_R = self.T(R=self.r)
        try:
            rsub = np.amin(self.r[T_R>=Tsub])
        except:
            rsub = 0.0
        return rsub
        
    def cs(self,R):
        return 5862*np.sqrt(self.T(1))*(R**(-0.25))
    
    def vk(self,R):
        return sqrt((Gconv*self.star['Ms']*Msun)/(R*AU))
    
    def H(self,R,fluid=0):
        def H_gas(R): #computes the pressure scale height for temp. distribution.
            return self.cs(R)*R/(self.vk(R))
        H0 = self.disk['H0'][fluid-1]
        R0 = self.disk['R0'][fluid-1]
        FI = self.disk['fi'][fluid-1]
        if fluid == 0:
            return H_gas(R)
        elif fluid == 1 and H0 > H_gas(R0): #if H of small dust > H of gas, go with Hgas
            return H_gas(R)
        else:
            return H0*(R/R0)**(1+FI)
        
    def make_grid(self):
        return np.meshgrid(self.r,self.theta,self.phi)
    
    def make_rz(self): #non-uniform spherical grid mapped onto cylindrical coordinates
        R,THETA,PHI = self.make_grid()
        R_CYL = np.clip(R*np.sin(THETA),a_min=1e-5,a_max=None)
        Z_CYL = R*np.cos(THETA)
        return R_CYL, Z_CYL
        
    def sig_profile(self,R,fluid=0): #fluid = 0 for gas, 1 for small dust (follows gas), 2 for large dust (settled)
        i = max(fluid,1)-1
        Rd = self.disk['Rdisk'][i]
        Rout = self.disk['Rout'][i]
        p = self.disk['p'][i]
        R0 = self.disk['R0'][i]
        frac = [1.-self.disk['Mfrac'][0]-self.disk['Mfrac'][1],self.disk['Mfrac'][0],self.disk['Mfrac'][1]]
        Mtot = self.disk['Mdisk']*frac[fluid]*Msun #assigns mass of disk component
        def sig_r(R):
            sig = (R/R0)**(p)*np.exp(-R/Rd)
            sig[R<R0] = 0.
            return sig
        def sig_int(R):
            if len(np.shape(R)) < 2:
                sig = sig_r(R)
                sig[R>Rout] = 0.
                return np.sum(sig*(R*AU)*(np.gradient(R))*AU)
            else:
                r = np.ravel(np.unique(R))
                sig = sig_r(r)
                sig[r>Rout] = 0.
                return np.sum(sig*(r*AU)*np.gradient(r)*AU)
        sig_0 = (Mtot)/(2.*pi*sig_int(R))
        return sig_0*sig_r(R)
    
    def rho_midplane(self,R,fluid=0):
        rho_mid = self.sig_profile(R,fluid=fluid)/(sqrt(pi*2.)*self.H(R,fluid=fluid)*AU)
        self.rhomax = np.log10(np.amax(rho_mid))
        return rho_mid
    
    def rho_disk(self,fluid=0):
        R_CYL,Z_CYL = self.make_rz()
        rho_mid = self.rho_midplane(R_CYL,fluid=fluid)
        h_mid = self.H(R_CYL,fluid=fluid)
        rho_vol = rho_mid*np.exp(-0.5*(Z_CYL/h_mid)**2)
        return rho_vol
    
    def plot_slice(self,rho,plot_params={'levels':np.arange(-27,-11,1)}):
        R_CYL,Z_CYL = self.make_rz()
        contourf(R_CYL[:,:,0],Z_CYL[:,:,0], np.log10(rho[:,:,0]),**plot_params)
        
    def plot_components(self,fluid=0):
        f,ax= subplots(1,3,constrained_layout=True)
        fluids = ['gas','small dust', 'large dust']
        f.suptitle(fluids[fluid])
        f.set_size_inches(9,3)
        components = {'Disk': self.rho_disk(fluid=fluid),'Envelope': self.rho_env(fluid=fluid), 'Disk + Envelope': self.rho_embedded(fluid=fluid)}
        for a,c in zip(ax, components.keys()):
            sca(a)
            self.plot_slice(components[c],{'levels':np.arange(self.rhomin,self.rhomax,1)})
            a.set_title(c,fontsize=14)
            a.set_xlim(self.grid['min'][0],self.grid['max'][0])
            if c == 'Disk':
                for j in np.arange(1,4):
                    a.plot(self.r, self.H(self.r)*j,color='C0',lw=1,label=str(j)+'H') #plot disk scale heights
                    a.legend()
            elif c == "Envelope":
                for j in np.linspace(np.radians(self.env['theta_min']),pi/2.,8):
                    r,th,rho = self.rho_stream(th0=j)
                    l = a.plot(r*np.sin(th),r*np.cos(th),color='C0',ls='dashed',lw=1) #plot streamlines
                a.legend(l,['streamlines'])
            a.set_ylim(a.get_xlim())
        colorbar(ax=ax)
        ax[1].set_xlabel('R [AU]')
        ax[0].set_ylabel('Z [AU]')
        
    def rho_stream(self, th0=pi/2.):
        th = self.theta[self.theta>=th0]
        Rc = self.env['Rc']
        Min = self.env['Min']
        Rin = Rc*np.sin(np.radians(self.env['theta_min']))**2 #inner landing radius of the stream
        Mfac = np.sqrt(1. - sqrt(Rin/Rc)) #adjusts the total mass to fit within the chosen streams
        if th0 != pi/2. and th0 != 0.:
            r = Rc*(np.sin(th0)**2)/(1. - (np.cos(th)/np.cos(th0)))
            r[0] = self.env['Rmax']
            r[-1] = Rc*(np.sin(th0)**2)
        elif th0 <= 1e-5:
            r = np.ones_like(th)*1
        else:
            r= np.ones_like(th)*Rc
        om = np.sqrt(Gconv*self.star['Ms']*Msun*(r*AU)**3)
        t1 = ((Min/Mfac)*Msun/yr)/(8*pi*om)
        t2 = np.sqrt(1+(np.cos(th)/np.cos(th0)))
        t3 = ((np.cos(th)/(2*np.cos(th0))) + (Rc/r)*(np.cos(th0)**2))
        rho1 = t1/(t2*t3)
        return r, th, rho1
    
    
    def rho_env(self,fluid=0):
        thstart = np.radians(self.env['theta_min'])
        rho0 = self.env['rho_0']
        R = np.array([])
        TH = np.array([])
        rho_tot = np.array([])
        for j in self.theta:
            r, th, rho1 = self.rho_stream(th0=j)
            R = np.append(R,r)
            TH = np.append(TH,th)
            if j <= thstart:
                rho_tot = np.append(rho_tot,rho1*0+rho0*r**(-self.env['exf']))
            else:
                rho_tot = np.append(rho_tot,rho1)
        # define grid.
        R_CYL = R*np.sin(TH)
        Z_CYL = R*np.cos(TH)
        r_cyl, z_cyl = self.make_rz()
        RHO = griddata((R_CYL, Z_CYL), rho_tot, (r_cyl[:,:,0],z_cyl[:,:,0]), method='nearest')
        RHO = np.reshape(RHO,(np.shape(RHO)[0],np.shape(RHO)[1],1))
        rho_vol = np.repeat(RHO,len(self.phi),axis=2)
        if fluid == 0:
            return rho_vol
        elif fluid == 1:
            return rho_vol*self.env['d2g']
        else:
            return rho_vol*0
           
    def rho_embedded(self,fluid=0):
        rho_d = self.rho_disk(fluid=fluid)
        if fluid < 2: # for gas and small dust
            rho_e = self.rho_env(fluid=fluid)
        else: # for large dust only
            rho_e = np.zeros_like(rho_d)
        floor = self.env['rho_amb'] 
        factor = min(fluid,1)*self.env['d2g']
        if factor != 0: 
            floor *= factor
        if fluid < 2: #for gas and small dust
            rho_both = np.clip(rho_d + rho_e,a_min=floor,a_max=None)
        else: # for large dust only
            rho_both = rho_d + rho_e
        if fluid != 0: #no dust inside sublimation radius
            R,THETA,PHI = self.make_grid()
            rsub = self.Rsub()
            rho_both[R <= rsub] = 10**(self.rhomin - 6)
        return rho_both
    
    def write_grid(self):
        iformat = 1
        grid_style = 0 # for "regular grid"
        coordsystem = 101 # between 1 and 200 for a spherical grid
        gridinfo = 0 
        incl_x, incl_y, incl_z = [1,1,1] # for a 3 dimensional grid
        nx,ny,nz = self.grid['N']
        r_edges = self.coords[0]*AU
        th_edges = self.coords[1]
        phi_edges = self.coords[2]
        header = str(iformat) + '\n' + str(grid_style) + '\n' + str(coordsystem) + '\n' + str(gridinfo) + '\n' + str(incl_x) + '\t' + str(incl_y) + '\t' + str(incl_z) + '\n' + str(nx) + '\t' + str(ny) + '\t' + str(nz) + '\n'
        with open(self.outdir+"amr_grid.inp","w") as f:
            f.write(header)
            for x,fmt  in zip([r_edges,th_edges,phi_edges],['%13.6e','%17.10e','%13.6e']):
                x.tofile(f, sep= '\t', format=fmt)
                f.write('\n')
        f.close()

    def write_dust_density(self):
        small_dust = self.rho_embedded(fluid=1)
        large_dust = self.rho_embedded(fluid=2)
        Nr = np.prod(np.array(self.grid['N']))
        with open(self.outdir+'dust_density.inp','w+') as f:
            f.write('1\n')                   # Format number
            f.write('%d\n'%(Nr))             # Number of cells
            f.write('2\n')                   # Number of dust species
            for dust in [small_dust,large_dust]:
                dust = dust.swapaxes(0,1) # radmc assumes 'ij' indexing for some reason
                data = dust.ravel(order='F')         # Create a 1-D view, fortran-style indexing
                data.tofile(f, sep='\n', format="%13.6e")
                f.write('\n')
        f.close()
        
    def write_wavelength(self,fname='',lam=None):
        if os.getcwd() != models_dir:
            os.chdir(models_dir)
        grid_data = rpy.reggrid.radmc3dGrid()
        grid_data.readWavelengthGrid()
        if os.getcwd() != self.outdir:
            os.chdir(self.outdir)
        if lam != None:
            grid_data.wav = lam
            grid_data.nwav = grid_data.wav.shape[0]
            grid_data.freq = c / grid_data.wav * 1e4
            grid_data.nfreq = grid_data.nwav
        grid_data.writeWavelengthGrid(fname=fname)
        os.chdir(parent_dir)
        
        
    def write_opacities(self): # note to self: figure out the number of dust species inputs and parameters
        with open(self.outdir+ 'dustopac.inp','w+') as f:
            f.write('2               Format number of this file\n')
            f.write('3               Nr of dust species\n')
            f.write('============================================================================\n')
            f.write('1               Way in which this dust species is read\n')
            f.write('0               0=Thermal grain\n')
            f.write('lg_maps_std        Extension of name of dustkappa_***.inp file\n')
            f.write('----------------------------------------------------------------------------\n')
            f.write('1               Way in which this dust species is read\n')
            f.write('0               0=Thermal grain\n')
            f.write('sm_maps_std        Extension of name of dustkappa_***.inp file\n')
            f.write('----------------------------------------------------------------------------\n')

    def write_main(self,nphot= 100000):        
        with open(self.outdir+'radmc3d.inp','w+') as f:
            f.write('nphot = %d\n'%(nphot))
            f.write('scattering_mode_max = 1\n')
            f.write('iranfreqmode = 1\n')
            f.write('modified_random_walk = 1 \n')
    
    #def RT_Tdust(self):
        #os.chdir(self.outdir)
        #file_list = ['amr_grid.inp', 'dust_density.inp','radmc3d.inp','wavelength_micron.inp','dustopac.inp','stars.inp']
        #func_list = [self.write_grid, self.write_dust_density, self.write_main, self.write_wavelength,self.write_opacities,self.write_star]
        #for file, func in zip(file_list,func_list):
         #   if os.path.exists(file) != True:
                #func()
        #try:
            #os.system('radmc3d mctherm')
        #except:
            #print('something went wrong')
        #return 0
    
    #def get_Tdust(self,fluid=1):
       # if os.getcwd() != self.outdir:
       #     os.chdir(self.outdir)
       # data = rpy.analyze.readData()
       # return data.dusttemp[:,:,:,fluid-1]
    
    #def get_rhodust(self,fluid=1):
       # if os.getcwd() != self.outdir:
       #     os.chdir(self.outdir)
       # data = rpy.analyze.readData()
       # return data.rhodust[:,:,:,fluid-1]
    
    #def make_rz_H(self,return_faces = False): #make a logarithmically spaced cylindrical grid
        #zmin = 1e-10
        #new_r = np.logspace(np.log10(self.grid['min'][0]), np.log10(self.grid['max'][0]), self.grid['N'][0])
        #zf_norm = np.append(zmin, np.logspace(-4,0,50)) #faces of the z-cells
        #zc_norm = 0.5*(zf_norm[1:] + zf_norm[:-1]) #50 points in the Z direction
        #if return_faces == True:
        #    R,Z = np.meshgrid(new_r,zf_norm) #returns the bin edges at the radius, for summing up columns
        #else:
        #    R,Z = np.meshgrid(new_r,zc_norm)
        #Z *= R/np.tan(np.radians(self.env['theta_min']+15))
        #return R,Z
    
    #def make_quadrant(self,quantity_3d,fill_value=0,order='F'): # to be used with the radmc 3d values
        #quantity_2d = quantity_3d[:,:,0] # phi=0 plane
        #r_cyl,z_cyl = self.make_rz()
        #r_new, z_new = self.make_rz_H()
        #quantity_2d_interp = griddata((r_cyl[:,:,0].flatten(),z_cyl[:,:,0].flatten()), quantity_2d.ravel(order=order), (r_new,z_new),fill_value=fill_value,method='linear')
        #return quantity_2d_interp
    
    
    
    
    
    