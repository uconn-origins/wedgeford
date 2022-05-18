from pylab import *
import radmc3dPy as rpy
import numpy as np
import os

from .. units import *
from .. outputs import *
from scipy.ndimage import gaussian_filter

class chemdisk:
    def __init__(self,output,chemdir = '/chemistry/environ/test1/'):
        self.input = output
        self.chemdir = self.input.m.parent_dir + chemdir
        if os.getcwd() != self.input.m.outdir:
            os.chdir(self.input.m.outdir)
        if output.rho == {}:
            output.rho2D()
        if output.T == {}:
            output.T2D()
        if 'gas' not in output.T.keys():
            calc_gas_T(output)
        self.data = {}
        
        try:
            os.mkdir(self.chemdir)
        except:
            print('directory exists - will overwrite current model if you write to it!')
    
    def regrid_disk(self):
        model = self.input.m
        rho_g = self.input.rho['gas']
        rho_d1 = self.input.rho['dust1'] 
        rho_d2 = self.input.rho['dust2']
        T_d = self.input.T['dust']
        T_g = self.input.T['gas'].T
        rhod1_2d = make_quadrant(model,rho_d1, fill_value = model.env['rho_amb'],order='F',smooth=False).T
        rhod2_2d = make_quadrant(model,rho_d2, fill_value = model.env['rho_amb'],order='F',smooth=False).T
        rhog_2d = make_quadrant(model,rho_g, fill_value = model.env['rho_amb'], order='F',smooth=False).T
        r,z = self.make_rz_H()
        #temperature gets smoothed before it goes in there as noise from radmc can make interpolation wonky
        Td_2d = make_quadrant(model,T_d,fill_value = 0,order='F')
        #fill in the bits outside the model quadrant with an ambient temperature from the model
        T_amb = model.T(np.sqrt(r**2 + z**2))
        Td_2d[Td_2d < 1] = T_amb[Td_2d < 1]
        Td_2d = Td_2d.T
        
        Tg_2d = make_quadrant(model,T_g,fill_value = 0,order='F')
        #fill in the bits outside the model quadrant with an ambient temperature from the model
        Tg_2d[Tg_2d < 1] = T_amb[Tg_2d < 1]
        Tg_2d = Tg_2d.T
        return rhog_2d, rhod1_2d, rhod2_2d, Td_2d, Tg_2d 
    
    def fill_data(self):
        rhog_2d, rhod1_2d, rhod2_2d, Td_2d, Tg_2d = self.regrid_disk()
        self.data['Tg'] = Tg_2d.flatten().reshape(np.shape(self.data['r']))
        self.data['rhod1'] = rhod1_2d.flatten().reshape(np.shape(self.data['r']))
        self.data['rhod2'] = rhod2_2d.flatten().reshape(np.shape(self.data['r']))
        
         
    def make_rz_H(self): #make a logarithmically spaced cylindrical grid
        model = self.input.m
        zmin = 1e-10
        new_r = np.logspace(np.log10(model.grid['min'][0]), np.log10(model.grid['max'][0]), model.grid['N'][0])
        zf_norm = np.append(zmin, np.logspace(-4,0,50)) #faces of the z-cells
        zc_norm = 0.5*(zf_norm[1:] + zf_norm[:-1]) #50 points in the Z direction
        R,Z = np.meshgrid(new_r,zc_norm)
        Z *= R/np.tan(np.radians(model.env['theta_min']))
        return R,Z
    
    def make_quadrant(self,quantity_2d,fill_value=0,order='F',smooth=True): # to be used with the radmc 3d values
        model = self.input.m
        r_cyl,z_cyl = model.make_rz()
        r_new, z_new = self.make_rz_H()
        if smooth == True:
            quantity_2d = gaussian_filter(quantity_2d, sigma=[1,4])
        quantity_2d_interp = griddata((r_cyl[:,:,0].flatten(),z_cyl[:,:,0].flatten()), quantity_2d.ravel(order=order), (r_new,z_new),fill_value=fill_value,method='linear',rescale=True)
        return quantity_2d_interp
    
    def write_out(self,outname='model',ndust=2):
        model = self.input.m
        if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
        r,z = self.make_rz_H()
        header1 = 'R(AU)    z(AU)   rhogas(g/cm3)  rhodust(g/cm^3)    T(K)   [    4 types of dust]      f_H\n'
        header2 = '-ForChem+Plots--------------------------------------------------------\n'

        rhog_2d, rhod1_2d, rhod2_2d, Td_2d, Tg_2d  = self.regrid_disk()
        rhod_2d = rhod1_2d + rhod2_2d
        with open(self.chemdir + outname+'.out',"w") as outfile:
            outfile.write(header1)
            outfile.write(header2)
            Nr = len(r.flatten())
            for n in np.arange(Nr):
                outfile.write('%10.4f %10.4f   %5.2e   %5.2e %7.1f %12.7f %12.7f %12.7f %12.7f   %7.2e\n'%(r.T.flatten()[n],z.T.flatten()[n],rhog_2d.flatten()[n], rhod_2d.flatten()[n] ,Td_2d.flatten()[n],1.0,0.0,0.0,0.0,1.0))
        return outname + '.out'

    def read_out(self,outname='model.out'): #convenience function for reading the .out file and returning the bits accessible as dictionary
        if os.getcwd() != self.chemdir:
            os.chdir(self.chemdir)
        if os.path.exists(outname) == False:
            self.write_out(outname)
        data = np.loadtxt(outname,skiprows=2,usecols=(0,1,2,3,4))
        r = data[:,0]
        z = data[:,1]
        shape = (len(np.unique(r)),50)
        rhog = data[:,2].reshape(shape)
        rhod = data[:,3].reshape(shape)
        T = data[:,4].reshape(shape)
        self.data = {'r':r.reshape(shape), 'z':z.reshape(shape),'rhog':rhog,'rhod':rhod,'Td':T}
        self.fill_data()

def plot_out(chem_disk,prop='rhog',outname='model.out',log=True,method=contourf,**pk): #plotting convenience function to inspect .out file contents
    if chem_disk.data == {}:
        chem_disk.read_out()
    data = chem_disk.data
    if log == True:
        method(data['r'],data['z'],np.log10(data[prop]),**pk)
    else:
        method(data['r'],data['z'],(data[prop]),**pk)
        
def plot_prechem(chem_disk,rlim=400): #convenience function to plot the dust radiative transfer real quick
    model = chem_disk.input.m
    X,Z = model.make_rz()
    X = X[:,:,0]
    Z = Z[:,:,0]
    f,ax = subplots(2,3,constrained_layout=True)
    f.set_size_inches(6,6)
    
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    if chem_disk.input.rho == {}:
        chem_disk.input.rho2D()
    if chem_disk.input.T == {}:
        chem_disk.input.T2D()
        
    c = ax[0,0].contourf(X, Z, np.log10(chem_disk.input.rho['gas'].T),levels=np.arange(model.rhomin,model.rhomax,1),extend='both')
   
    ax[0,0].set_ylabel('z [au]')
    ax[0,0].set_title(r'Original $\rho_g$',fontsize=11)
    cb = colorbar(c,ax=ax[:,0],location='bottom',aspect=10,ticks=np.array(c.levels[::2]).astype(int))
    cb.set_label(r'$\rho$ [$\mathrm{g/cm^{-3}}$]')
    
    sca(ax[1,0])
    plot_out(chem_disk,prop='rhog',cmap='viridis',levels=np.arange(model.rhomin,model.rhomax,1),log=True,extend='both')
    ax[1,0].set_ylabel('z [au]')
    ax[1,0].set_title(r'Regridded $\rho_g$ ',fontsize=11)
    
    
    c = ax[0,1].contourf(X, Z, np.log10(chem_disk.input.rho['dust1'].T + chem_disk.input.rho['dust2'].T),levels=np.arange(model.rhomin,model.rhomax,1),extend='both')
    ax[0,1].set_title(r'Original $\rho_d$',fontsize=11)
    cb = colorbar(c,ax=ax[:,1],location='bottom',aspect=10,ticks=np.array(c.levels[::2]).astype(int))
    cb.set_label(r'$\rho$ [$\mathrm{g/cm^{-3}}$]')
    
    sca(ax[1,1])
    plot_out(chem_disk,prop='rhod',cmap='viridis',levels=np.arange(model.rhomin,model.rhomax,1),log=True,extend='both')
    ax[1,1].set_title(r'Regridded $\rho_d$',fontsize=11)
    
    
    
    from scipy.ndimage import gaussian_filter
    smooth_T = gaussian_filter(chem_disk.input.T['dust'].T, sigma=[4,1])
    c2 = ax[0,2].contourf(X, Z, smooth_T,levels=np.linspace(5,125,19),cmap='twilight_shifted',extend='both')


    ax[0,2].set_title(r'Original $T_d$',fontsize=11)
    cb2 = colorbar(c2,ax=ax[:,2],location='bottom',ticks=np.array(c2.levels[::3]).astype(int),aspect=10)
    cb2.set_label(r'$T \ \mathrm{[K]}$')
    
    sca(ax[1,2])
    plot_out(chem_disk,prop='Td',cmap='twilight_shifted',levels=np.linspace(5,125,19),log=False,extend='both')
    
    ax[1,0].set_xlabel('r [au]')
    ax[1,1].set_xlabel('r [au]')
    ax[1,2].set_xlabel('r [au]')
    
    ax[1,2].set_title(r'Regridded $T_d$',fontsize=11)
   
    
    for a in ax[0,:]:
        a.set_xlim(0,rlim)
        a.set_ylim(0,rlim)
        a.xaxis.set_major_locator(MultipleLocator(int(rlim/5)))
        a.yaxis.set_major_locator(MultipleLocator(int(rlim/5)))
    
    for a in ax[1,:]:
        a.set_xlim(0,rlim)
        a.set_ylim(0,rlim)
        a.xaxis.set_major_locator(MultipleLocator(int(rlim/5)))
        a.yaxis.set_major_locator(MultipleLocator(int(rlim/5)))
    f.savefig(chem_disk.chemdir+'regrid.pdf')
        