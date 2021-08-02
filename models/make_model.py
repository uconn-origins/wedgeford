from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from models.constants import *

class model:
    def __init__(self,stellar_params, disk_params, envelope_params,grid_params, dust_params, RT_params, outdir = '/m1_test/'):
        self.star = stellar_params
        self.disk = disk_params
        self.env = envelope_params
        self.grid = grid_params
        self.dust = dust_params
        self.outdir = os.getcwd()+outdir
        self.coords = [0,0,0]
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
        c = ax[0].contourf(X, Z, np.log10(data.rhodust[:,:,0,0].T),levels=np.arange(-26,-11,1))
        ax[0].contour(X, Z, np.log10(data.rhodust[:,:,0,1].T),levels=np.arange(-26,-11,1),colors='white')
        ax[0].set_xlabel('r [AU]')
        ax[0].set_ylabel('z [AU]')
        cb = colorbar(c,ax=ax[0],location='top')
        cb.set_label(r'$\rho$ [$g/cm^{-3}$]')
        
        c2 = ax[1].contourf(X, Z, data.dusttemp[:,:,0,0].T,levels=[0,5,10,15,25,50,100,200],cmap='magma')
        ax[1].set_xlabel('r [AU]')
        cb2 = colorbar(c2,ax=ax[1],location='top')
        cb2.set_label('T [K]')
        return f,ax
    
    def make_rz_uniform(self): #make a cylindrical grid for the chemical models
        r_uni = np.logspace(np.log10(self.grid['min'][0]), np.log10(self.grid['max'][0]), self.grid['N'][0])
        N_th = self.grid['N'][1]
        th_cav_wide = np.radians(self.env['theta_min'])
        r_wall,th_wall,rho_wall = self.rho_stream(th_cav_wide)
        th_cav = np.amin(th_wall[r_wall <= np.amax(r_uni)]) # getting edge of the cavity lower down
        th_polar = np.linspace(th_cav/2.,th_cav,2) #thetas in the cavity
        th_eq = np.linspace(th_cav,pi/2.,50-2 ) #thetas not in the cavity
        th_uni = np.unique(np.append(th_polar,th_eq)) # all the thetas together
        R_CYL,TH = np.meshgrid(r_uni,th_uni)
        Z_UNI = (R_CYL/np.sin(TH))*np.cos(TH)
        z_norm = (2.*np.amax(r_uni)/np.amax(Z_UNI,axis=0))*np.ones_like(Z_UNI)
        return R_CYL, Z_UNI*z_norm
    
    
    def make_quadrant(self,fluid=0):
        rho_3d = self.rho_embedded(fluid=fluid)
        rho_2d = rho_3d[:,:,0] # phi=0 plane
        r_cyl,z_cyl = self.make_rz()
        r_new, z_new = self.make_rz_uniform()
        rho_2d_interp = griddata((r_cyl[:,:,0].flatten(),z_cyl[:,:,0].flatten()), rho_2d.ravel(), (r_new[1:,:],z_new[1:,:]),fill_value=self.env['rho_amb'])
        return rho_2d_interp

        
    def sig_profile(self,R,fluid=0): #fluid = 0 for gas, 1 for small dust (follows gas), 2 for large dust (settled)
        i = max(fluid,1)-1
        Rd = self.disk['Rdisk'][i]
        Rout = self.disk['Rout'][i]
        p = self.disk['p'][i]
        R0 = self.disk['R0'][i]
        frac = [1.-self.disk['Mfrac'][0]-self.disk['Mfrac'][1],self.disk['Mfrac'][0],self.disk['Mfrac'][1]]
        Mtot = self.disk['Mdisk']*frac[fluid]*Msun #assigns mass of disk component
        def sig_r(R):
            return (R/R0)**(p)*np.exp(-R/Rd)
        def sig_int(R):
            if len(np.shape(R)) < 2:
                sig = sig_r(R)
                sig[R<R0] = 0.
                sig[R>Rout] = 0.
                return np.sum(sig*(R*AU)*(np.gradient(R))*AU)
            else:
                r = np.ravel(np.unique(R))
                sig = sig_r(r)
                sig[r<R0] = 0.
                sig[r>Rout] = 0.
                return np.sum(sig*(r*AU)*np.gradient(r)*AU)
        sig_0 = (Mtot)/(2.*pi*sig_int(R))
        return sig_0*sig_r(R)
    
    def rho_midplane(self,R,fluid=0):
        return self.sig_profile(R,fluid=fluid)/(sqrt(pi*2.)*self.H(R,fluid=fluid)*AU)
    
    def rho_disk(self,fluid=0):
        R_CYL,Z_CYL = self.make_rz()
        rho_mid = self.rho_midplane(R_CYL,fluid=fluid)
        h_mid = self.H(R_CYL,fluid=fluid)
        rho_vol = rho_mid*np.exp(-0.5*(Z_CYL/h_mid)**2)
        return rho_vol
    
    def plot_slice(self,rho,plot_params={'levels':np.arange(-27,-11,1)}):
        R_CYL,Z_CYL = self.make_rz()
        contourf(R_CYL[:,:,0],Z_CYL[:,:,0], np.log10(rho[:,:,0]),**plot_params)
        colorbar()
        
    def rho_stream(self, th0=pi/2.):
        th = self.theta[self.theta>=th0]
        Rc = self.env['Rc']
        Min = self.env['Min']
        if th0 != pi/2. and th0 != 0.:
            r = Rc*(np.sin(th0)**2)/(1. - (np.cos(th)/np.cos(th0)))
            r[0] = self.env['Rmax']
            r[-1] = Rc*(np.sin(th0)**2)
        elif th0 <= 1e-5:
            r = np.ones_like(th)*1
        else:
            r= np.ones_like(th)*Rc
        om = np.sqrt(Gconv*self.star['Ms']*Msun*(r*AU)**3)
        t1 = (Min*Msun/yr)/(8*pi*om)
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
        else:
            return rho_vol*self.env['d2g']
        
    
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
            f.write('%d\n'%(Nr))           # Number of cells
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
            f.write('1               Way in which this dust species is read\n')
            f.write('0               0=Thermal grain\n')
            f.write('smdsharp        Extension of name of dustkappa_***.inp file\n')
            f.write('----------------------------------------------------------------------------\n')

    def write_main(self,nphot= 100000):        
        with open(self.outdir+'radmc3d.inp','w+') as f:
            f.write('nphot = %d\n'%(nphot))
            f.write('scattering_mode_max = 1\n')
            f.write('iranfreqmode = 1\n')