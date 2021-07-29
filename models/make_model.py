from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from constants import *

class model:
    def __init__(self,stellar_params, disk_params, envelope_params,grid_params, dust_params, outdir = '/m1_test/'):
        self.star = stellar_params
        self.disk = disk_params
        self.env = envelope_params
        self.grid = grid_params
        self.dust = dust_params
        self.outdir = outdir
        self.coords = [0,0,0]
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
        if fluid < 2:
            if fluid == 1 and self.disk['H0'][0] == H_gas(self.disk['R0'][0]):
                return H_gas
            else:
                return self.disk['H0'][fluid-1]*(R/self.disk['R0'][fluid-1])**(1+self.disk['fi'][fluid-1])
        else:
            return self.disk['H0'][fluid-1]*(R/self.disk['R0'][fluid-1])**(1+self.disk['fi'][fluid-1])
        
    def make_grid(self):
        return np.meshgrid(self.r,self.theta,self.phi)
    
    def make_rz(self):
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
            return (R/R0)**(p)*np.exp(-R/Rd)
        def sig_int(R):
            sig = sig_r(R)
            sig[R<R0] = 0.
            sig[R>Rd] = 0.
            return np.sum(sig*(R*AU)*(np.gradient(R))*AU)
        
        sig_0 = (Mtot)/(2*pi*sig_int(R))
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
        with open(outdir+"amr_grid.inp","w") as f:
            f.write(header)
            for x,fmt  in zip([r_edges,th_edges,phi_edges],['%13.6e','%17.10e','%13.6e']):
                x.tofile(f, sep= '\t', format=fmt)
                f.write('\n')
        f.close()

    
    def write_dust_density(self):
        small_dust = self.rho_embedded(fluid=1)
        large_dust = self.rho_embedded(fluid=2)
        Nr = np.prod(np.array(self.grid['N']))
        with open(outdir+'dust_density.inp','w+') as f:
            f.write('1\n')                   # Format number
            f.write('%d\n'%(Nr))           # Number of cells
            f.write('2\n')                   # Number of dust species
            for dust in [small_dust,large_dust]:
                data = dust.ravel(order='F')         # Create a 1-D view, fortran-style indexing
                data.tofile(f, sep='\n', format="%13.6e")
                f.write('\n')
        f.close()