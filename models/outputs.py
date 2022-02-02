from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from models.units import *


class out: #Class that handles and stores the outputs and data for radmc3d after you've made the initial model structure
    def __init__(self,model):
        if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
        self.data = rpy.data.radmc3dData()
        self.rad_data = rpy.radsources.radmc3dRadSources()
        self.r = self.data.grid.x/AU
        self.theta = self.data.grid.y
        self.phi = self.data.grid.z
        self.m = model # original model access
        # dictionaries that load 2D slices in
        self.rho = {}
        self.T = {} 
        self.J = {}
        self.rad = {}
        
    def thermal_RT(self):
        os.chdir(self.m.outdir)
        file_list = ['amr_grid.inp', 'dust_density.inp','radmc3d.inp','wavelength_micron.inp','stars.inp','external_source.inp']
        func_list = [self.m.write_grid, self.m.write_dust_density, self.m.write_main, self.m.write_wavelength,self.m.write_star,self.m.calc_ISRF]
        for file, func in zip(file_list,func_list):
            if os.path.exists(file) != True:
                func()
        try:
            os.system('radmc3d mctherm')
        except:
            print('something went wrong')
        return 0
    
    def HE_RT(self,nphot = 1000000):
        os.chdir(self.m.outdir)
        try:
            os.system('radmc3d mcmono nphot_mono {}'.format(nphot))
        except:
            print('something went wrong with radmc3d')
        return 0
    
    def rz(self):
        R,TH = np.meshgrid(self.r,self.theta)
        X = R*np.sin(TH)
        Z = R*np.cos(TH)
        return X,Z
        
    def rho2D(self):
        self.data.readDustDens() #self.data.rhodust is now loaded in
        d2g = self.m.disk['Mfrac'][0] #dust to gas ratio for small dust
        self.rho['gas'] = self.data.rhodust[:,:,0,0]/d2g 
        ndust = np.shape(self.data.rhodust[0,0,0,:])[-1]
        for n in range(ndust):
            self.rho['dust' + str(n+1)] = self.data.rhodust[:,:,0,n]
    
    def T2D(self):
        self.data.readDustTemp()#self.data.Tdust is now loaded in
        self.T['dust'] = self.data.dusttemp[:,:,0,0] #only load first dust temp
        
    def rad_source(self,field='star',fname=''):
        self.rad_data.readStarsinp(fname=fname)
        self.rad[field] = {'lam': self.rad_data.grid.wav, 'fnu':self.rad_data.fnustar}
        
    def Jnu(self,field = 'UV',fname='mean_intensity.out'):
        shape = (len(self.r),len(self.theta),len(self.phi))
        def read_intensity(fname=fname,freq=None):
            with open(fname,'r') as f:
                header_info = f.readlines()[:4]
            header = {'iformat':None,'nrcells':None,'nfreq':None,'freq':None}
            for line,key in zip(header_info,header.keys()):
                header[key] = np.array([i for i in line.strip('\n').split(' ') if i != '']).astype(float)
                if len(header[key]) == 1:
                    header[key] = int(header[key][0])
            if freq != None and isinstance(freq,list):
                indices = np.array(freq)
                indices = indices[indices < header['nfreq']]
            else:
                indices = np.arange(0,header['nfreq'])
            Jnu = np.zeros_like(indices,dtype='object')
            nu = header['freq'][indices]
            for j,k in zip(indices,range(len(nu))):
                jj = np.genfromtxt(fname,skip_header=4+header['nrcells']*j,max_rows=header['nrcells'])
                Jnu[k] = jj.reshape(shape,order='F')[:,:,0]
            return nu, Jnu
        self.J[field] = read_intensity(fname)
        
        
    def plot_dustRT(self): #convenience function to plot the dust radiative transfer real quick
        if os.getcwd() != self.m.outdir:
            os.chdir(self.m.outdir)
        X,Z = self.rz()
        f,ax = subplots(1,2,constrained_layout=True)
        f.set_size_inches(9,4.5)
        self.rho2D()
        self.T2D()
        c = ax[0].contourf(X, Z, np.log10(self.rho['dust1'].T),levels=np.arange(self.m.rhomin,self.m.rhomax,1))
        ax[0].set_xlabel('r [AU]')
        ax[0].set_ylabel('z [AU]')
        ax[0].set_title('Input Model Density',fontsize=12)
        cb = colorbar(c,ax=ax[0],location='bottom')
        cb.set_label(r'$\rho$ [$g/cm^{-3}$]')
        
        c2 = ax[1].contourf(X, Z, self.T['dust'].T,levels=[0,10,25,50,100,200],cmap='plasma')
        ax[1].contour(X, Z, np.log10(self.rho['dust1'].T),levels=np.arange(-19,-11,1),colors='white')
        ax[1].set_xlabel('r [AU]')
        ax[1].set_title('Dust Temperature',fontsize=12)
        cb2 = colorbar(c2,ax=ax[1],location='bottom')
        cb2.set_label('T [K]')