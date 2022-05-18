from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata

from . units import *
from . write_radmc_files import *



class out: #Class that handles and stores the outputs and data for radmc3d after you've made the initial model structure
    def __init__(self,model):
        if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
        grid = rpy.reggrid.radmc3dGrid()
        grid.readGrid()
        self.data = rpy.data.radmc3dData(grid=grid)
        self.r = self.data.grid.x/AU
        self.theta = self.data.grid.y
        self.phi = self.data.grid.z
        self.m = model #model object is encoded here
        # dictionaries that load 2D slices in
        self.rho = {}
        self.T = {} 
        self.J = {}
        
    def rz(self):
        R,TH = np.meshgrid(self.r,self.theta)
        X = R*np.sin(TH)
        Z = R*np.cos(TH)
        return X,Z
        
    def rho2D(self):
        model = self.m
        if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
        self.data.readDustDens() #self.data.rhodust is now loaded in
        d2g = self.m.disk['Mfrac'][0] #dust to gas ratio for small dust
        self.rho['gas'] = model.rho_embedded(fluid=0).swapaxes(0,1)[:,:,0]
        ndust = np.shape(self.data.rhodust[0,0,0,:])[-1]
        for n in range(ndust):
            self.rho['dust' + str(n+1)] = self.data.rhodust[:,:,0,n]
    
    def T2D(self):
        if os.getcwd() != self.m.outdir:
            os.chdir(self.m.outdir)
        self.data.readDustTemp()#self.data.Tdust is now loaded in
        self.T['dust'] = self.data.dusttemp[:,:,0,0] #only load first dust temp
        
    def Jnu(self,field = 'UV',fname='mean_intensity.out'):
        shape = (len(self.r),len(self.theta),len(self.phi))
        if os.getcwd() != self.m.outdir:
            os.chdir(self.m.outdir)
        with open(fname,'r') as f:
            header_info = f.readlines()[:4]
        header = {'iformat':None,'nrcells':None,'nfreq':None,'freq':None}
        for line,key in zip(header_info,header.keys()):
            header[key] = np.array([i for i in line.strip('\n').split(' ') if i != '']).astype(float)
            if len(header[key]) == 1:
                header[key] = int(header[key][0])
        lam = (c/header['freq'])*1e4
        if field.lower() == 'uv':
            indices = np.where((lam >= uv_min) & (lam <= uv_max))[0]
        elif field.lower() == 'xray':
            indices = np.where((lam >= xray_min) & (lam <= xray_max))[0]
        else:
            indices = np.arange(header['nfreq'])
        Jnu = np.zeros_like(indices,dtype='object')
        nu = header['freq'][indices]
        for j,k in zip(indices,range(len(nu))):
            jj = np.genfromtxt(fname,skip_header=4+header['nrcells']*j,max_rows=header['nrcells'])
            Jnu[k] = jj.reshape(shape,order='F')[:,:,0]
        self.J[field] = nu,Jnu
        
    def integrate_intensity(self,field='uv',photons=False):
        freq = self.J[field][0][::-1]
        integrated_J = self.J[field][1][-1]
        efactor = 1
        for j in np.arange(1,len(freq)):
            if photons == True:
                efactor = h*freq[j]
            integrated_J = np.dstack((integrated_J,self.J[field][1][-j]/efactor))
        integrated_J = np.trapz(integrated_J, x=np.expand_dims(freq,(0,1)),axis=-1)
        return integrated_J


def prep_thermal_transfer(output,infallheat=False,nphot=500000,mrw=1,maxtau=5):
    model = output.m
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    file_list = ['amr_grid.inp', 'dust_density.inp','wavelength_micron.inp','stars.inp']
    func_list = [write_grid, write_dust_density,write_wavelength,write_star]
    for file, func in zip(file_list,func_list):
        #if os.path.exists(file) != True:
        func(model)
    model.disk_transport(setnew=infallheat)
    xlam,xfnu = insert_xray_radiation(model)
    write_opacities(model,update=False)
    ilam,ifnu = calc_ISRF(model,write=model.rad['G0'])
    write_viscous_heatsource(model,write=model.rad['viscous_heating'])
    write_main(model, scat=2,nphot=nphot,mrw=mrw,maxtau=maxtau)

def do_thermal_transfer(output,nt=4,prep=False,**prepkw):
    model = output.m
    if prep == True:
        prep_thermal_transfer(output,**prepkw)
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    os.system('radmc3d mctherm setthreads {}'.format(nt))
        
def prep_he_transfer(output):
    model = output.m
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    write_opacities(model,update=True)
    write_main(model, scat=2, mrw=1, maxtau=20)
    
def do_he_transfer(output,nphot=100000,nt=4,prep=False):
    model = output.m
    if prep == True:
        prep_he_transfer(output)
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    os.system('radmc3d mcmono nphot_mono {} setthreads {}'.format(nphot,nt))
        
def calc_gas_T(output):
    model = output.m
    shock = model.env['shock']
    if output.rho == {}:
        output.rho2D()
    if output.T == {}:
        output.T2D()
    if 'uv' not in output.J.keys():
        output.Jnu(field='uv')
    nG0 = np.log10(np.clip(output.integrate_intensity(field='uv')/(2.*pi*G0),a_min=10**-0.5,a_max=10**6.5))
    nH = np.log10(2.*output.rho['gas']/(mu*mh))
    nH_ = np.clip(nH,a_min=1,a_max=7)
    
    T_UV = [model.T_heat(i,j)[0] for i,j in zip(nH_.flatten(), nG0.flatten())]
    T_UV = np.reshape(np.array(T_UV),np.shape(nH)).T
    
    if shock == True:
        T_shock = model.solve_envelope(prop='Tg')[:,:,0]
    else:
        T_shock = np.zeros_like(T_UV)
        
    Tcrit = 130*(10**nH.T/1e10)**(0.3)
    Tgas = np.maximum(T_UV,T_shock)
    
    #if gas temp is below the critical temperature, its coupled to the dust
    coupled = np.where((Tgas-Tcrit)<= 1.0) 
    Tgas[coupled] = 0.0
    T_shock[coupled] = 0.0
    T_UV[coupled] = 0.0
    #assuming that the floor of the gas temperature in uncoupled regions is the dust temperature
    Tgas = np.maximum(Tgas,output.T['dust'].T)
        
    output.T['uv'] = T_UV
    output.T['shock'] = T_shock
    output.T['crit'] = Tcrit
    output.T['gas'] = Tgas
    
    
def prep_line_transfer(output,molecules={'names':['co'],'abundances':[1e-4],'lines':[2]},heated=False):
    model = output.m
    get_molecule_info(names=molecules['names'],outdir=model.outdir)
    write_molecule_density(model, names=molecules['names'],abundances = molecules['abundances'])
    write_lines(model)
    write_velocities(model)
    if heated == True:
        calc_gas_T(output)
        write_gas_temperature(model,output.T['gas'])
        
def insert_hydrodust(output,small_hydro,large_hydro):
    model = output.m
    if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
    output.data.readDustDens() #self.data.rhodust is now loaded in
    small_dust_all0 = output.data.rhodust[:,:,:,0].swapaxes(0,1)
    large_dust_all0 = output.data.rhodust[:,:,:,1].swapaxes(0,1)
    small_dust_disk0 = model.rho_disk(fluid=1)
    large_dust_disk0 = model.rho_disk(fluid=2)
    small_dust = small_dust_all0 - small_dust_disk0 + small_hydro
    large_dust = large_dust_all0 - large_dust_disk0 + large_hydro
    Nr = np.prod(np.array(model.grid['N']))
    with open(model.outdir+'dust_density.inp','w+') as f:
        f.write('1\n')                   # Format number
        f.write('%d\n'%(Nr))             # Number of cells
        f.write('2\n')                   # Number of dust species
        for dust in [small_dust,large_dust]:
            data = dust.swapaxes(0,1).ravel(order='F') # radmc assumes 'ij' indexing for some reason, Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')
    f.close()