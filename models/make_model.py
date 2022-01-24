from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
# units and conversions
Msun = 1.9891e33 # g
Rsun = 69.6e9 #cm
Lsun  = 3.8525e33 #erg/s
AU = 1.49598e13 #cm
yr = 3.14e7 #seconds
cmtokm = 10**(-5) #converts cm to km
R_mu = 36149835 # (R/mu) for 2.4 g/mol in PPD
Gconv = 6.6743e-8 # cgs
S0conv = (Msun/(AU**2)) #S0cgs = S0code*S0conv
sigsb = 5.67e-5 #Stefan Boltzmann constant
c = 2.99792458e10 # speed of light
mu = 2.36
h = 6.6260755e-27 #planck
mh = 1.67e-24 # mass of hydrogen cgs
Habing = 1.83590e+08  # Photons/cm^2/s between 930-2000
G0 = 2.7e-3 # erg/cm^2/s^-1

parent_dir = os.getcwd()

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
        
    def T(self,R):
        return self.star['Ts']*np.sqrt((self.star['Rs']*Rsun)/(R*AU))
    
    def cs(self,R):
        return 5862*np.sqrt(self.T(1))*(R**(-0.25))
    
    def vk(self,R):
        return sqrt((Gconv*self.star['Ms']*Msun)/(R*AU))
    
    def H(self,R,fluid=0):
        def H_gas(R):
            return self.cs(R)*R/(self.vk(R))
        H0 = self.disk['H0'][fluid-1]
        R0 = self.disk['R0'][fluid-1]
        FI = self.disk['fi'][fluid-1]
        if fluid == 0:
            return H_gas(R)
        elif fluid == 1 and H0 > H_gas(R0):
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
        #colorbar(location='top')
        
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
        if fluid < 2:
            rho_e = self.rho_env(fluid=fluid)
        else:
            rho_e = np.zeros_like(rho_d)
        floor = self.env['rho_amb'] 
        factor = min(fluid,1)*self.env['d2g']
        if factor != 0: 
            floor *= factor
        if fluid < 2:
            rho_both = np.clip(rho_d + rho_e,a_min=floor,a_max=None)
        else:
            rho_both = rho_d + rho_e
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
        
        
    def write_wavelength(self):
        lam1     = 0.1e0
        lam2     = 7.0e0
        lam3     = 25.e0
        lam4     = 1.0e4
        n12      = 20
        n23      = 100
        n34      = 30
        lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
        lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
        lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
        lam      = np.concatenate([lam12,lam23,lam34])
        nlam     = lam.size
        with open(self.outdir + 'wavelength_micron.inp','w+') as f:
            f.write('%d\n'%(nlam))
            for value in lam:
                f.write('%13.6e\n'%(value))
        
        
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
    
    def RT_Tdust(self):
        os.chdir(self.outdir)
        file_list = ['amr_grid.inp', 'dust_density.inp','radmc3d.inp','wavelength_micron.inp','dustopac.inp']
        func_list = [self.write_grid, self.write_dust_density, self.write_main, self.write_wavelength,self.write_opacities]
        for file, func in zip(file_list,func_list):
            if os.path.exists(file) != True:
                func()
        try:
            os.system('radmc3d mctherm')
        except:
            print('something went wrong with radmc3d')
        return 0
    
    def plot_RT(self):
        if os.getcwd() != self.outdir:
            os.chdir(self.outdir)
        data = rpy.analyze.readData()
        data_r = data.grid.x/AU
        data_th = data.grid.y
        R,TH = np.meshgrid(data_r,data_th)
        X = R*np.sin(TH)
        Z = R*np.cos(TH)
        f,ax = subplots(1,2,constrained_layout=True)
        f.set_size_inches(9,4.5)
        c = ax[0].contourf(X, Z, np.log10(data.rhodust[:,:,0,0].T),levels=np.arange(self.rhomin,self.rhomax,1))
        ax[0].set_xlabel('r [AU]')
        ax[0].set_ylabel('z [AU]')
        ax[0].set_title('Input Model Density',fontsize=12)
        cb = colorbar(c,ax=ax[0],location='bottom')
        cb.set_label(r'$\rho$ [$g/cm^{-3}$]')
        
        c2 = ax[1].contourf(X, Z, data.dusttemp[:,:,0,0].T,levels=[0,10,25,50,100,200],cmap='plasma')
        ax[1].contour(X, Z, np.log10(data.rhodust[:,:,0,1].T),levels=np.arange(-19,-11,1),colors='white')
        ax[1].set_xlabel('r [AU]')
        ax[1].set_title('Dust Temperature',fontsize=12)
        cb2 = colorbar(c2,ax=ax[1],location='bottom')
        cb2.set_label('T [K]')
    
    def get_Tdust(self,fluid=1):
        if os.getcwd() != self.outdir:
            os.chdir(self.outdir)
        data = rpy.analyze.readData()
        return data.dusttemp[:,:,:,fluid-1]
    
    def get_rhodust(self,fluid=1):
        if os.getcwd() != self.outdir:
            os.chdir(self.outdir)
        data = rpy.analyze.readData()
        return data.rhodust[:,:,:,fluid-1]
    
    def make_rz_H(self): #make a cylindrical grid for the chemical models based on scale heights
        r_uni = np.logspace(np.log10(self.grid['min'][1]), np.log10(self.grid['max'][-1]), self.grid['N'][0])
        H_lim = self.H(r_uni)*8
        new_r = r_uni[H_lim<np.amax(r_uni)] #select new radial limit for 8 scale heights up
        N_th = self.grid['N'][1]
        z_norm = np.append([0],np.logspace(-3,0,49)) #50 points in the Z direction
        R,Z = np.meshgrid(new_r,z_norm)
        Z *= R/np.tan(np.radians(self.env['theta_min']+15))
        return R,Z
    
    def make_quadrant(self,quantity_3d,fill_value=0,order='F'): # to be used with the radmc 3d values
        quantity_2d = quantity_3d[:,:,0] # phi=0 plane
        r_cyl,z_cyl = self.make_rz()
        r_new, z_new = self.make_rz_H()
        quantity_2d_interp = griddata((r_cyl[:,:,0].flatten(),z_cyl[:,:,0].flatten()), quantity_2d.ravel(order=order), (r_new,z_new),fill_value=fill_value,method='linear')
        return quantity_2d_interp
    
    
    
    
    
    