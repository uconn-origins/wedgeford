from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy import interpolate

from . units import *


def read_par(file="infall.par"):
    disk_par = {'AspectRatio': 0, 'Sigma0': 0 , 'Alpha': 1.0e-7, 'SigmaSlope':0,
               'FlaringIndex':0, 'rc':0 , 'rcd':0 , 'minfall':0,
               'rin': 0, 'rout':0,'Nx':0, 'Ny': 0, 'Xmin':0, 'Xmax':0, 'Ymin':0, 'Ymax':0,
                'Ntot':0, 'Ninterm':0, 'DT':0, 'OutputDir': '0', 'Invstokes1':0, 'Invstokes2':0, 
                'epsilon1': 0, 'epsilon2':0, 'epsilonism':0.01}
    if os.path.isfile(file):
        file1 = open(file, 'r') 
        count = 0

        while True: 
            count += 1

            # Get next line from file 
            line = file1.readline() 

            # if line is empty 
            # end of file is reached 
            if not line: 
                break
            line_keys = line.split()
            for key in disk_par:
                if key in line_keys:
                    if key != 'OutputDir':
                        disk_par[key] = float(line_keys[1])
                    else:
                        disk_par[key] = line_keys[1]
        file1.close()
        return disk_par
    else:
        print('Parfile not found!')

def load_vars(t,OutputDir,var='dens',fluid='gas'): 
    var_list = ['dens','energy','vx','vy']
    if var in var_list:
        filebase = fluid + var
        filename = filebase + str(t) + '.dat'
        #data = np.fromfile(OutputDir+filename)
        data = np.memmap(OutputDir+filename,dtype=float,mode='c')
        return data

def load_scalar(OutputDir,var='mass'):
    file = var+'.dat'
    tvals = np.loadtxt(OutputDir+file)[:,0]
    mvals = np.loadtxt(OutputDir+file)[:,1]
    return tvals, mvals

class fargo_disk:
    def __init__(self,OutputDir,Mstar=1,R0 = 1, nx=512,ny=512):
        S0conv = (Msun*Mstar/((R0*AU)**2)) #S0cgs = S0code*S0conv
        Tconv = 1./(2.*pi) # 2PI = 1 orbit at R0 = 1 yr T(years) = Tcode*Tconv
        Tconv_cgs = (1./(2.*pi))*yr #T(cgs) = Tcode*Tconv_cgs
        Vconv_cgs = sqrt(Gconv*Mstar*Msun/(R0*AU)) #V(cgs) = Vcode*Vconv_cgs
        Vconv_kms = Vconv_cgs/10**5
        Econv_cgs = Msun*Mstar*(R0*AU)**2/(Tconv_cgs)**2
        self.pars = read_par(file=OutputDir+'summary0.dat')
        self.dust = {}
        if self.pars != 0:
            self.nx = int(self.pars['Nx'])
            self.ny = int(self.pars['Ny'])
            self.dtoutput = self.pars['DT']*self.pars['Ninterm']*Tconv
            if self.pars['Invstokes1'] != 0:
                self.dust['dust1'] = 1./self.pars['Invstokes1']
                self.Mdust = {}
                if self.pars['Invstokes2'] != 0:
                    self.dust['dust2'] = 1./self.pars['Invstokes2']
        else:
            self.dtoutput = Tconv*10*10*2.*pi
            self.nx = nx
            self.ny = ny
        self.ri = np.loadtxt(OutputDir+'domain_y.dat')[4:self.ny+4+1]
        self.phii = np.loadtxt(OutputDir+'domain_x.dat')[0:self.nx+1]
        self.r = (self.ri[1:] + self.ri[:-1]) / 2.
        self.phi = (self.phii[1:] + self.phii[:-1]) / 2
        t, Md = load_scalar(OutputDir)
        if self.dust != {}:
            for key in self.dust.keys():
                self.Mdust[key] = load_scalar(OutputDir+'monitor/'+key+'/',var='mass')
        self.Mgas = (t*Tconv,Md)
        self.Mstar = Mstar
        self.var_units = {'sigma': S0conv, 'vphi': Vconv_cgs, 'vr': Vconv_cgs, 'time': self.dtoutput}
        self.var_names = {'sigma': 'dens', 'vphi': 'vx', 'vr': 'vy'}
        self.OutputDir = OutputDir
            
    def get(self,time,prop,fluid='gas',dim=None):
        if prop in self.var_names:
            if os.path.isfile(self.OutputDir+fluid+self.var_names[prop]+str(time)+'.dat'):
                data = load_vars(time,self.OutputDir, self.var_names[prop],fluid=fluid).reshape(self.ny,self.nx)*self.var_units[prop]
                if dim == 1:
                    return np.average(data,axis=1)
                else:
                    return data
            else:
                return {}
        else:
            return {}
        
    def sound_speed(self,dim=None):
        if self.pars != 0:
            aspectratio = self.pars['AspectRatio']
            flare_exponent = self.pars['FlaringIndex']
        else:
            aspectratio = 0.046
            flare_exponent = 0.25
        if dim == 1:
            R = self.r
        else:
            R,PHI = np.meshgrid(self.r,self.phi)
        cs = aspectratio*(R**(flare_exponent)).T*self.vkep(dim=dim)
        return cs
    
    def T(self,dim=None):
        return (self.sound_speed(dim=dim)**2)*(mu*m_H)/kB
    
    def vkep(self,dim=None,r='all'):
        if r == 'all':
            R,PHI = np.meshgrid(self.r, self.phi)
        else:
            R,PHI = np.meshgrid(r, self.phi)
        vkep_2d = np.sqrt(Gconv*self.Mstar*Msun/(R*AU)).T
        if dim == 1:
            return np.average(vkep_2d,axis=1)
        else:
            return vkep_2d

    def make_2d(self, arr_1d):
        R,PHI = np.meshgrid(self.r,self.phi)
        ARR_2D = arr_1d*np.ones_like(R)
        return ARR_2D.T
        
    def dr(self):
        return np.diff(self.ri)
    
    def dphi(self):
        return np.diff(self.phii)
    
    def xy(self):
        R,PHI = np.meshgrid(self.r, self.phi)
        x = R*np.cos(PHI)
        y = R*np.sin(PHI)
        return x,y
    
    def alpha_rey(self,time,dim=None):
        cs = self.sound_speed()
        sig = self.get(time,'sigma')
        num = sig*self.dev_ave(time,'vphi')*self.dev_ave(time,'vr')
        den = sig*cs**2
        alpha = num/den
        if dim == 0:
            return np.abs(np.average(num)/np.average(den))
        elif dim == 1:
            return np.abs(np.average(num, axis=1)/np.average(den,axis=1))
        else:
            return alpha
        
    def eta(self,time,dim=None):
        eta = -1.*self.dev_v_kep(time,dim=None)/disk.vkep(dim=None)
        if dim == 1:
            return np.average(eta,axis=1)
        else:
            return eta
    
    def t_yr(self,time,dim=1):
        return time*self.dtoutput
    
    def t_snap(self,physical_time):
        return int(physical_time/self.dtoutput)
    
    def at_R(self,time,func_R,rval,fluid='gas'):
        if 'kwargs' in func_R:
            kwargs = func_R['kwargs']
        else:
            kwargs = {}
        if ['dim'] in kwargs:
            if kwargs['dim'] == 1:
                dim = 1
                R = self.r
        else:
            dim = None
            R,PHI = np.meshgrid(self.r,self.phi)
        if isinstance(func_R['func'],str):
            field = self.get(time,func_R['func'],dim = dim,fluid=fluid)
        elif callable(func_R):
            field = func_R['func'](time,**kwargs)
        r_diff = np.abs(R-rval)
        if dim == 1:
            closest_index = np.argmin(r_diff)
        else:
            closest_index = np.argmin(r_diff,axis=1)
        return field[closest_index]
    
def params_from_fargo(simdir,par={'Ms':1,'Rs':2,'Mfrac':[0.01,0.01]}): #crib parameters straight from the FARGO model
    disk = fargo_disk(simdir,Mstar=par['Ms'])
    par['Ts'] = disk.T[0]*(disk.r[0]/par['Rs'])**(2*disk.pars['FlaringIndex'])
    par['d2g'] = disk.pars['epsilonism']
    par['Min'] = disk.pars['minfall']*par['Ms']
    par['Rc'] = disk.pars['rc']
    par['theta_min'] = np.arccos(disk.pars['rin']/disk.pars['rc'])
    par['Mdisk'] = disk.Mgas[1][0]*par['Ms']
    if disk.dust != {}:
        par['Mfrac'] = [disk.pars['epsilon1'],disk.pars['epsilon2']]
    par['Rdisk'] = [disk.pars['Rdisk'],disk.pars['Rdisk']]
    par['p'] = [-disk.pars['SigmaSlope'],-disk.pars['SigmaSlope']]
    return par


def disk_from_fargo(output,simdir,snapshot,fluid=0):
    model = output.m
    disk = fargo_disk(simdir,Mstar=model.star['Ms'])
    R,THETA,PHI = model.make_grid()
    X_C = R*np.sin(THETA)*np.cos(PHI)
    Y_C = R*np.sin(THETA)*np.sin(PHI)
    Z_C = R*np.cos(THETA)
    R_CYL = np.sqrt(X_C**2 + Y_C**2)
    sigma_2d = disk.get(time=snapshot,prop='sigma')
    x_disk, y_disk = disk.xy()
    if fluid < 2: #
        sigma_2d_scaled = sigma_2d
        if fluid == 1:
            sigma_2d_scaled *= model.disk['Mfrac'][0]
    else:
        a1 = disk.dust['dust1']*sigma_2d/((pi/2.)*model.dust['rho_si'])
        a2 = disk.dust['dust2']*sigma_2d/((pi/2.)*model.dust['rho_si'])
        sigma_2d_1 = disk.get(time=snapshot,prop='sigma',fluid='dust1')
        sigma_2d_2 = disk.get(time=snapshot,prop='sigma',fluid='dust2')
        sigma_2d_tot = sigma_2d_1+sigma_2d_2
        amin = model.dust['amin'][-1]*1e-4 #micron to cm
        amax = model.dust['amax'][-1]*1e-4
        powerlaw = model.dust['apow']
        amid = (np.log10(a1) + np.log10(a2))/2.
        mmid = np.log10(sigma_2d_2*((10**(amid)/a2)**(-(powerlaw-3))))
        a0 = np.log10(amax)
        m0 = mmid + (powerlaw - 3.)*(amid - a0)
        sigma_2d_scaled = (10**m0)*((10**a0)**(powerlaw-3))*(amax**(-(powerlaw-4)) - amin**(-(powerlaw-4)))/(-(powerlaw-4))
    sig_xy = interpolate.griddata((x_disk.flatten(),y_disk.flatten()), sigma_2d_scaled.T.flatten(), (X_C[-1,:,:],Y_C[-1,:,:]), method='linear',fill_value=0,rescale=True)
    hmid = model.H(R_CYL,fluid=fluid)
    rho_vol = np.clip(np.expand_dims(sig_xy,axis=0)/(sqrt(pi*2.)*hmid*AU),a_min=model.env['rho_amb'],a_max=None)*np.exp(-0.5*(Z_C/hmid)**2)
    return rho_vol