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
        self.T0 = {}
        self.H0 = {}
        
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
        self.rho['gas'] = np.average(model.rho_embedded(fluid=0).swapaxes(0,1),axis=-1)
        ndust = np.shape(self.data.rhodust[0,0,0,:])[-1]
        for n in range(ndust):
            self.rho['dust' + str(n+1)] = np.average(self.data.rhodust[:,:,:,n],axis=-1)
    
    def T2D(self):
        if os.getcwd() != self.m.outdir:
            os.chdir(self.m.outdir)
        self.data.readDustTemp()#self.data.Tdust is now loaded in
        self.T['dust'] = np.average(self.data.dusttemp[:,:,:,0],axis=-1)
        #only load first dust temp
        
    def T3D(self):
        if os.getcwd() != self.m.outdir:
            os.chdir(self.m.outdir)
        self.data.readDustTemp()#self.data.Tdust is now loaded in
        return self.data.dusttemp[:,:,:,0] #only load first dust temp
    
    def rho3D(self,fluid=0):
        model = self.m
        if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
        if fluid == 0:
            return model.rho_embedded(fluid=0).swapaxes(0,1)
        self.data.readDustDens() #self.data.rhodust is now loaded in
        if fluid > 0:
            return self.data.rhodust[:,:,:,fluid-1]
        
    def Jnu(self,field = 'UV',fname='mean_intensity.out',average=True,iphi=0):
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
            if average == True:
                Jnu[k] = np.average(jj.reshape(shape,order='F'),axis=-1)
            else:
                Jnu[k] = jj.reshape(shape,order='F')#[:,:,iphi]
        self.J[field] = nu,Jnu
        
    def integrate_intensity(self,field='uv',photons=False):
        freq = self.J[field][0][::-1]
        integrated_J = self.J[field][1][-1]
        efactor = 1
        for j in np.arange(1,len(freq)):
            if photons == True:
                efactor = h*freq[j]
            integrated_J = np.dstack((integrated_J,self.J[field][1][-j]/efactor))
        if (self.J[field][1][0]).squeeze().ndim > 2:
            integrated_J = np.reshape(integrated_J, (len(self.r),len(self.theta),len(self.phi),len(freq)))
            integrated_J = np.trapz(integrated_J, x=np.expand_dims(freq,(0,1,2)),axis=-1)
        else:
            integrated_J = np.trapz(integrated_J, x=np.expand_dims(freq,(0,1)),axis=-1)
        return integrated_J


def prep_thermal_transfer(output,infallheat=False,nphot=500000,mrw=1,maxtau=5):
    model = output.m
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    file_list = ['amr_grid.inp', 'dust_density.inp','wavelength_micron.inp','stars.inp','stars_thermal.inp']
    func_list = [write_grid, write_dust_density,write_wavelength,write_star,insert_xray_radiation]
    for file, func in zip(file_list,func_list):
        if os.path.exists(file) != True:
            print('writing new file:' + file)
            func(model)
    if os.path.exists('heatsource.inp') != True and model.rad['viscous_heating'] != False:
        model.disk_transport(setnew=infallheat)
        write_viscous_heatsource(model,write=model.rad['viscous_heating'])
    if os.path.exists('dustopac.inp') != True:
        #xlam,xfnu = insert_xray_radiation(model)
        write_opacities(model,update=False)
    if os.path.exists('external_source.inp') != True:
        ilam,ifnu = calc_ISRF(model,write=model.rad['G0'])
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
    Tgas[coupled] = output.T['dust'].T[coupled]
    T_shock[coupled] = 0.0
    T_UV[coupled] = 0.0
    #assuming that the floor of the gas temperature in uncoupled regions is the dust temperature
    #Tgas = np.maximum(Tgas,output.T['dust'].T)
        
    output.T['uv'] = T_UV
    output.T['shock'] = T_shock
    output.T['crit'] = Tcrit
    output.T['gas'] = Tgas
    

def calc_gas_T3D(output):
    model = output.m
    shock = model.env['shock']
    nH = np.log10(2.*output.rho3D()/(mu*mh))
    
    T_dust = output.T3D()
    T_gas = T_dust
    
    n_min = 1
    n_max = 7
    n_num = 60
    
    g_min = -0.5
    g_max = 6.5
    g_num = 60
    
    Ts = model.T_heat(np.linspace(n_min,n_max,n_num),np.linspace(g_min,g_max,g_num))
    
    if output.J['uv'][1].ndim != 3:
        output.Jnu(field='uv', average=False)
    nG0 = np.log10(output.integrate_intensity(field='uv')/(2.*pi*G0))
    # should be:
    # nG0 = np.log10(4*pi*output.integrate_intensity(field='uv')/G0)
    
    
    i = np.clip((n_num*(nH - n_min)/(n_max-n_min)).astype(int),a_min = 0, a_max=n_num-1)
    j = np.clip((g_num*(nG0 - g_min)/(g_max-g_min)).astype(int),a_min=0, a_max=g_num-1)
    
    T_UV = Ts[j,i]
    
    if shock == True:
        T_shock = model.solve_envelope(prop='Tg')*model.stream_mask
    else:
        T_shock = np.zeros_like(T_UV)
    
    T_crit = 130*(10**(nH)/1e10)**(0.3)
    
    Tgas = np.maximum(T_UV,T_shock)
    
    uncoupled = np.where((Tgas - T_crit) >= 1.0)
    
    T_gas[uncoupled] = Tgas[uncoupled]
    
    return T_gas.swapaxes(0,1)
        
        
    
def prep_line_transfer(output,molecules={'names':['co'],'abundances':[1e-4],'lines':[2]},heated=False):
    model = output.m
    get_molecule_info(names=molecules['names'],outdir=model.outdir)
    write_molecule_density(model, names=molecules['names'],abundances = molecules['abundances'])
    write_lines(model)
    write_velocities(model)
    if heated == True:
        calc_gas_T(output)
        write_gas_temperature(model,output.T['gas'])
        
