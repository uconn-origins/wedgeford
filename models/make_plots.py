from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from models.units import *
from models.outputs import *

##### plotting convenience functions #####


def plot_slice(self,rho,log=True,plot_params={'levels':np.arange(-27,-11,1)}):
    R_CYL,Z_CYL = self.make_rz()
    if np.shape(rho) != np.shape(R_CYL[:,:,0]):
        rho2d = rho[:,:,0]
    else:
        rho2d = rho
    ax = gca()
    if log == True:
        im=ax.contourf(R_CYL[:,:,0],Z_CYL[:,:,0], np.log10(rho2d),**plot_params)
    else:
        im=ax.contourf(R_CYL[:,:,0],Z_CYL[:,:,0], rho2d,**plot_params)
    return im
     
        
def plot_components(self,fluid=0):
    f,ax= subplots(1,3,constrained_layout=True)
    fluids = ['gas','small dust', 'large dust']
    f.suptitle(fluids[fluid])
    f.set_size_inches(9,3)
    components = {'Disk': self.rho_disk(fluid=fluid),'Envelope': self.rho_env(fluid=fluid), 'Disk + Envelope': self.rho_embedded(fluid=fluid)}
    for a,c in zip(ax, components.keys()):
        sca(a)
        im=plot_slice(self,rho=components[c],log=True,plot_params={'levels':np.arange(self.rhomin,self.rhomax,1)})
        a.set_title(c,fontsize=14)
        a.set_xlim(self.grid['min'][0],self.grid['max'][0])
        if c == 'Disk':
            for j in np.arange(1,4):
                a.plot(self.r, self.H(self.r)*j,color='C0',lw=1,label=str(j)+'H') #plot disk scale heights
                a.legend()
        elif c == "Envelope":
            for j in np.linspace(np.radians(self.env['theta_min']),pi/2.,8):
                streamline = self.stream(th0=j,shock=False)
                r = streamline['path'][0]
                th = streamline['path'][1]
                l = a.plot(r*np.sin(th),r*np.cos(th),color='C0',ls='dashed',lw=1) #plot streamlines
            a.legend(l,['streamlines'])
        a.set_ylim(a.get_xlim())
    colorbar(im,ax=ax)
    ax[1].set_xlabel('R [AU]')
    ax[0].set_ylabel('Z [AU]')
    

def plot_dustRT(self): #convenience function to plot the dust radiative transfer real quick
    if os.getcwd() != self.outdir:
        os.chdir(self.outdir)
    X,Z = self.make_rz()
    f,ax = subplots(1,2,constrained_layout=True)
    f.set_size_inches(9,4.5)
    opt = out(self)
    opt.rho2D()
    opt.T2D()
    c = ax[0].contourf(X[:,:,0], Z[:,:,0], np.log10(opt.rho['dust1'].T),levels=np.arange(self.rhomin,self.rhomax,1))
    ax[0].set_xlabel('r [AU]')
    ax[0].set_ylabel('z [AU]')
    ax[0].set_title('Input Model Density',fontsize=12)
    cb = colorbar(c,ax=ax[0],location='bottom')
    cb.set_label(r'$\rho$ [$g/cm^{-3}$]')
        
    c2 = ax[1].contourf(X[:,:,0], Z[:,:,0], opt.T['dust'].T,levels=[0,10,20,30,40,50,100,150,200],cmap='plasma')
    ax[1].contour(X[:,:,0], Z[:,:,0], np.log10(opt.rho['dust1'].T),levels=np.arange(self.rhomax-4,self.rhomax,1),colors='white')
    ax[1].set_xlabel('r [AU]')
    ax[1].set_title('Dust Temperature',fontsize=12)
    cb2 = colorbar(c2,ax=ax[1],location='bottom')
    cb2.set_label('T [K]')
    
    
def plot_flux_components(self,fname=''):
    if os.getcwd() != self.outdir:
        os.chdir(self.outdir)
    opts = out(self)
    opts.rad_data.nstar = 1
    opts.rad_data.getStarSpectrum(grid=opts.data.grid,ppar=self.radmcpar)
    fnu_s = opts.rad_data.fnustar
    lam = opts.rad_data.grid.wav
    nu = c/(1e-4*lam)
    nufnu_s = np.reshape(nu,np.shape(nu)[0])*np.reshape(fnu_s,np.shape(fnu_s)[0])
    loglog(lam, nufnu_s,label='Star')
    
    opts.rad_data.getSpotSpectrum(grid=opts.data.grid,ppar=self.radmcpar)
    fnu_acc = opts.rad_data.fnuspot
    nufnu_acc = np.reshape(nu,np.shape(nu)[0])*np.reshape(fnu_acc,np.shape(fnu_acc)[0])
    loglog(lam, nufnu_acc,label='UV Accretion')
    
    opts.rad_data.getAccdiskSpectra(grid=opts.data.grid,ppar=self.radmcpar)
    fnu_d = opts.rad_data.fnuaccdisk
    nufnu_d = np.reshape(nu,np.shape(nu)[0])*np.reshape(fnu_d,np.shape(fnu_d)[0])
    loglog(lam, nufnu_d,label='Viscous Disk')
    
    opts.rad_data.getAccdiskTemperature(grid=opts.data.grid,ppar=self.radmcpar)
    
    opts.rad_data.incl_accretion=True
    res=opts.rad_data.getTotalLuminosities(readInput=False)
    print(res)
    loglog(lam, nufnu_s+nufnu_acc+nufnu_d, label='SED (estimated)',color='black')
    xlabel(r'$\log \ \lambda \ \mathrm{[\mu m]}$')
    ylabel(r'$\log \ \nu F_{\nu} \ \mathrm{[ergs \ cm^{-2} \ s^{-1}]}$')
    title(r'$T_{\rm s}$' + ' = {} K,'.format(self.star['Ts']) + r' $\ \dot{M}_{\rm s}$' + '= {} '.format(self.star['accrate']) + r'$\mathrm{M_{\odot} \ yr^{-1}}$',size=12)
    legend()
    ylim(np.amin(nufnu_s)*10000,)
    
def plot_shock_model(self):
    colors=cm.get_cmap('viridis',10)
    thstart = np.radians(self.env['theta_min'])
    f,ax = subplots(3)
    f.set_size_inches(3,6)
    for j,k in zip(np.linspace(thstart+0.01,pi/2-0.01,10),range(10)):
        streams = self.stream(th0=j,shock=True)
        rland = streams['path'][0][-1]
        r = streams['path'][0]
        th = streams['path'][1]
        z = r*np.cos(th)
        rhod = self.rho_disk_profile(r=r*np.sin(th),z=r*np.cos(th))
        Tcrit = 130*((streams['rho']+rhod)/(mu*mh*1e10))**(0.3)
        ax[0].loglog(z,streams['rho'],color=colors(k))
        ax[1].loglog(z,(streams['Tg']+self.T(streams['path'][0])),color=colors(k))
        ax[2].loglog(z,np.sqrt(streams['v'][0]**2+streams['v'][1]**2)/1e5,color=colors(k))
    ax[0].set_ylim(10**self.rhomin,)
    ax[1].set_ylim(5,)
    ax[2].set_ylim(0.2,)

    for a in ax:
        a.set_xlim(1,300)
    ax[-1].set_xlabel('z [au]')
    ax[0].set_ylabel(r'$\rho_{\rm in} \ \mathrm{[g \ cm^{-3}]}$')
    ax[1].set_ylabel(r'$T_{\rm gas}$ [K]')
    ax[2].set_ylabel(r'$v_{\rm in} \ \mathrm{[km \ s^{-1}]}$')
    ax[0].set_title('Shock Model')
    return f