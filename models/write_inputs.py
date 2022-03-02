from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from models.units import *
from models.outputs import *

parent_dir = os.getcwd()
models_dir = parent_dir+'/models/'

def write_star(self):
    if os.getcwd() != self.outdir:
        os.chdir(self.outdir)
    if os.path.exists('wavelength_micron.inp') != True:
        write_wavelength(self)
    grid = rpy.reggrid.radmc3dGrid()
    grid.readWavelengthGrid()
    rad_data = rpy.radsources.radmc3dRadSources(ppar=self.radmcpar,grid=grid)
    rad_data.incl_accretion = True
    if os.getcwd() != self.outdir:
        os.chdir(self.outdir)
    rad_data.writeStarsinp(ppar=self.radmcpar) #writes a stars.inp in the new directory
    os.chdir(parent_dir)
    
def calc_ISRF(self,lam=None,gnorm=1,write=True):
    from scipy.interpolate import interp1d
    if os.getcwd() != self.outdir:
        os.chdir(self.outdir)
    if os.path.exists('wavelength_micron.inp') != True:
        write_wavelength(self)
    grid_data = rpy.reggrid.radmc3dGrid()
    grid_data.readWavelengthGrid()
    gnorm = self.rad['G0']
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


def write_viscous_heatsource(self,write=True):
    if os.getcwd() != self.outdir:
        os.chdir(self.outdir)
    opts = out(self)
    opts.rad_data.nstar = 1
    opts.rad_data.incl_accretion=True
    opts.rad_data.getAccdiskTemperature(grid=opts.data.grid,ppar=self.radmcpar)
    tacc = opts.rad_data.tacc
    facc = sigsb*(tacc**4)
    lacc = facc/(self.H(self.r)*AU)
    lacc[self.r< self.disk['R0'][0]] = np.amin(lacc)
    R_CYL,Z_CYL = self.make_rz()
    lacc = np.reshape(lacc,(1,np.shape(lacc)[0],1))
    D_disk = lacc*np.exp(-0.5*(Z_CYL/self.H(R_CYL))**2)
    #im=plot_slice(self,rho=D_disk,plot_params={'levels':np.arange(-20,-11,1),'cmap':'magma'})
    #colorbar(im)
    Nr = np.prod(np.array(self.grid['N']))
    if write == True:
        with open(self.outdir+'heatsource.inp','w+') as f:
            f.write('1\n')                   # Format number
            f.write('%d\n'%(Nr))             # Number of cells
            heat = D_disk.swapaxes(0,1) # radmc assumes 'ij' indexing for some reason
            heat = heat.ravel(order='F')         # Create a 1-D view, fortran-style indexing
            heat.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')
    os.chdir(parent_dir)


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
        f.write('scattering_mode_max = 2\n')
        f.write('iranfreqmode = 1\n')
        f.write('modified_random_walk = 1 \n')
