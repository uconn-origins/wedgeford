from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata

from . units  import *
from . radiation import *
from . dust import *



####### basic spatial grid/density info for thermal MC
def write_grid(model):
    iformat = 1
    grid_style = 0 # for "regular grid"
    coordsystem = 101 # between 1 and 200 for a spherical grid
    gridinfo = 0 
    incl_x, incl_y, incl_z = [1,1,1] # for a 3 dimensional grid
    nx,ny,nz = model.grid['N']
    r_edges = model.coords[0]*AU
    th_edges = model.coords[1]
    phi_edges = model.coords[2]
    header = str(iformat) + '\n' + str(grid_style) + '\n' + str(coordsystem) + '\n' + str(gridinfo) + '\n' + str(incl_x) + '\t' + str(incl_y) + '\t' + str(incl_z) + '\n' + str(nx) + '\t' + str(ny) + '\t' + str(nz) + '\n'
    with open(model.outdir+"amr_grid.inp","w") as f:
        f.write(header)
        for x,fmt  in zip([r_edges,th_edges,phi_edges],['%13.6e','%17.10e','%13.6e']):
            x.tofile(f, sep= '\t', format=fmt)
            f.write('\n')
            
def write_dust_density(model,envelope=True):
    if envelope == True:
        small_dust = model.rho_embedded(fluid=1)
        large_dust = model.rho_embedded(fluid=2)
    else:
        small_dust = model.rho_disk(fluid=1)
        large_dust = model.rho_disk(fluid=2)
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
    
####### basic wavelength, radiation source information for thermal MC
                     
def write_wavelength(model,fname='',lam=[]):
    templates_dir = model.models_dir+'templates/'
    if os.getcwd() != templates_dir:
        os.chdir(templates_dir)
    grid_data = rpy.reggrid.radmc3dGrid()
    grid_data.readWavelengthGrid()
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    xlam = np.logspace(np.log10(xray_min),np.log10(xray_max),20)
    if size(lam)==0:
        new_lam = np.append(xlam,grid_data.wav) 
    else:
        new_lam = lam
    grid_data.wav = new_lam
    grid_data.nwav = grid_data.wav.shape[0]
    grid_data.freq = c / grid_data.wav * 1e4
    grid_data.nfreq = grid_data.nwav
    grid_data.writeWavelengthGrid(fname=fname)
    os.chdir(model.parent_dir)
    
def write_star(model):
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    if os.path.exists('wavelength_micron.inp') != True:
        write_wavelength(model)
    grid = rpy.reggrid.radmc3dGrid()
    grid.readWavelengthGrid()
    rad_data = rpy.radsources.radmc3dRadSources(ppar=model.radmcpar,grid=grid)
    rad_data.incl_accretion = True
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    rad_data.writeStarsinp(ppar=model.radmcpar) #writes a stars.inp in the new directory
    os.chdir(model.parent_dir)
    
def calc_ISRF(model,lam=[],write=True):
    from scipy.interpolate import interp1d
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    if os.path.exists('wavelength_micron.inp') != True:
        write_wavelength(model)
    grid_data = rpy.reggrid.radmc3dGrid()
    grid_data.readWavelengthGrid()
    gnorm = model.rad['G0']
    #goal is to write the ISRF onto the same wavelength grid as the stellar source
    if size(lam) == 0:
        lam = grid_data.wav
        nu = grid_data.freq
    else:
        nu = c/(lam*1e-4)
    templates_dir = model.models_dir+'templates/'
    os.chdir(templates_dir)
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
        with open(model.outdir+'external_source.inp','w') as f:
            f.write('2 \n')
            f.write(str(len(lam))+ '\n')
            f.write('\n')
            lam.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')
            f.write('\n')
            fnu.tofile(f, sep='\n', format="%13.6e")
    os.chdir(model.parent_dir)
    return lam,fnu

def write_viscous_heatsource(model,write=True):
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    grid = rpy.reggrid.radmc3dGrid()
    grid.readWavelengthGrid()
    grid.readSpatialGrid()
    rad_data = rpy.radsources.radmc3dRadSources(ppar=model.radmcpar,grid=grid)
    rad_data.nstar = 1
    rad_data.incl_accretion = True
    rad_data.getAccdiskTemperature(grid=grid,ppar=model.radmcpar)
    tacc = rad_data.tacc
    facc = sigsb*(tacc**4)
    lacc = facc/(np.sqrt(2.*pi)*model.H(model.r)*AU)
    lacc[model.r > model.env['Rc']] = np.amin(lacc)
    lacc[model.r< model.disk['R0'][0]] = np.amin(lacc)
    R_CYL,Z_CYL = model.make_rz()
    lacc = np.reshape(lacc,(1,np.shape(lacc)[0],1))
    D_disk = lacc*np.exp(-0.5*(Z_CYL/model.H(R_CYL))**2)
    Nr = np.prod(np.array(model.grid['N']))
    if write == True:
        with open(model.outdir+'heatsource.inp','w') as f:
            f.write('1\n')                   # Format number
            f.write('%d\n'%(Nr))             # Number of cells
            heat = D_disk.swapaxes(0,1).ravel(order='F')# radmc assumes 'ij' indexing for some reason Create a 1-D view, fortran-style indexing
            heat.tofile(f, sep='\n', format="%13.6e")
    os.chdir(model.parent_dir)

###### inserts the x-ray accretion spectrum into the stars file, should do this before thermal MC
    
def insert_xray_radiation(model,nx=5,nu=5,write=True): #writes the wavelength and spectrum info for the new fields
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    grid = rpy.reggrid.radmc3dGrid()
    grid.readWavelengthGrid()
    rad_data = rpy.radsources.radmc3dRadSources(ppar=model.radmcpar,grid=grid)
    rad_data.nstar = 1
    rad_data.incl_accretion = True
    rad_data.getStarSpectrum(grid=grid,ppar=model.radmcpar)
    rad_data.getSpotSpectrum(grid=grid,ppar=model.radmcpar)
    
    
    lam_ = np.squeeze(rad_data.grid.wav)
    fnu_ = np.squeeze(rad_data.fnustar) + np.squeeze(rad_data.fnuspot)
       
    x_lam = np.logspace(np.log10(xray_min),np.log10(xray_max),nx)
    u_lam = np.append(uv_min,np.linspace(lam_lya,uv_max,nu-1)) #includes Lya line
    he_lam = np.append(x_lam,u_lam)
    write_wavelength(model,fname='mcmono_wavelength_micron.inp',lam=he_lam)
    
    #add lya flux
    fpeak = Lya_line(model,fLya=model.rad['fLya'])
    fnu_[np.argmin(np.abs(lam_-lya_lam))] += fpeak
    
    #find xray wavelengths
    x_lam = lam_[lam_ <= xray_max]
    nx = size(x_lam)
    
    #compute spectrum at those wavelengths
    x_lam, x_fnu = Xray_accretion(model,fX=None,nx=nx)
    fnu_[lam_<= xray_max] += x_fnu
    
    if write == True:
        therm_star_file = model.outdir+'stars.inp'
        if os.path.exists(therm_star_file):
            os.rename(therm_star_file, model.outdir+'stars_thermal.inp')

        with open(model.outdir+'stars_thermal.inp',"r") as f:
            f.readline()
            f.readline()
            line = f.readline()
        f.close()

        with open(model.outdir+'stars.inp','w') as f:
            f.write('2 \n')
            f.write('1 \t')
            f.write(str(len(lam_))+ '\n')
            f.write(line)
            f.write('\n')
            lam_.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')
            f.write('\n')
            fnu_.tofile(f, sep='\n', format="%13.6e")
    return x_lam, x_fnu

###### writes the opacity information and updates the dust opacities to include the x-ray photoelectric cross-section
###### run with update = False before thermal MC
###### run with update = True before high energy mcmono

def write_opacities(model,ndust=2,filenames=['',''],update=True): 
    iformat = 3
    exts = {}
    for fluid,fname in zip(range(1,ndust+1),filenames):
        if fname == '':
            kappa_file = "{}dustkappa_dust-{}.inp".format(model.outdir,fluid)
            exts[str(fluid)] = 'dust-{}'.format(fluid)
        else:
            kappa_file = '{}dustkappa_{}.inp'.format(model.outdir,fname)
            exts[str(fluid)] = fname
        if os.path.exists(kappa_file) != True and fname != '':
            print('{}: Dust opacity file not found in model directory'.format(kappa_file))
        elif os.path.exists(kappa_file) != True and fname == '':
            print('Running optool to generate new opacities')
            run_optool(model,fluid=fluid,na=50)
            exts[str(fluid)] = 'dust-{}'.format(fluid)
        if update == True:
            kappa_file = "{}dustkappa_dust-{}x.inp".format(model.outdir,fluid)
            if os.path.exists(kappa_file) != True:
                lam, kabs, kscat, g  = amend_kappa(model,fluid=fluid,filename=fname)
                arrays = np.stack((lam,kabs,kscat,g),axis=-1)
                nlam = len(lam)
                header = '{} \n{} \n'.format(iformat,nlam)
                np.savetxt(kappa_file,arrays,header=header,fmt='%12.6e',comments='',delimiter='\t')
            exts[str(fluid)] = 'dust-{}x'.format(fluid)
    with open(model.outdir+ 'dustopac.inp','w') as f:
        f.write('2               Format number of this file\n')
        f.write('{}               Nr of dust species\n'.format(ndust))
        f.write('============================================================================\n')
        for fluid in range(1,ndust+1):
            f.write('1               Way in which this dust species is read\n')
            f.write('0               0=Thermal grain\n')
            f.write('{}    Extension of name of dustkappa_***.inp file\n'.format(exts[str(fluid)]))
            f.write('----------------------------------------------------------------------------\n')
            
##################################################################################            
            
def write_main(model,nphot= 100000,scat=2,mrw =1,maxtau=15,noheat=1,dust=1,lines=0):        
    with open(model.outdir+'radmc3d.inp','w') as f:
        f.write('nphot = {}\n'.format(nphot))
        f.write('scattering_mode_max = {}\n'.format(scat))
        f.write('iranfreqmode = 1\n')
        f.write('modified_random_walk = {}\n'.format(mrw))
        f.write('mc_scat_maxtauabs = {}\n'.format(maxtau))
        f.write('tgas_eq_tdust = {}\n'.format(noheat))
        f.write('incl_dust = {}\n'.format(dust))
        f.write('incl_lines = {}\n'.format(lines))
            
##### for line transfer ##########################################################
def write_lines(model,names=['co']):
    with open(model.outdir+'lines.inp','w+') as f:
        f.write('2\n')                   # Format number
        f.write('{}\n'.format(len(names)))             # Number of molecules
        for name in names:
            f.write('{} leiden 0 0 0\n'.format(name))

def write_velocities(model,envelope = True):
    if envelope == True:
        velocities = model.v_embedded()
    else:
        velocities = model.v_disk()
    Nr = np.prod(np.array(self.grid['N']))
    with open(model.outdir+'gas_velocity.inp','w+') as f:
        f.write('1\n')                   # Format number
        f.write('%d\n'%(Nr))             # Number of cells
        vr = velocities['r'].swapaxes(0,1).ravel(order='F')
        vth = velocities['theta'].swapaxes(0,1).ravel(order='F')
        vphi = velocities['phi'].swapaxes(0,1).ravel(order='F')
        for i,j,k in zip(vr,vth,vphi):
            f.write("%9e %9e %9e\n" % (i, j, k))
            
def write_gas_temperature(model,Tgas):
    if Tgas.ndim != 3:
        np.repeat(np.expand_dims(Tgas,axis=-1),len(self.phi),axis=-1)
    Nr = np.prod(np.array(self.grid['N']))
    with open(model.outdir+'gas_temperature.inp','w+') as f:
        f.write('1\n')                   # Format number
        f.write('%d\n'%(Nr))             # Number of cells
        data = Tgas.swapaxes(0,1).ravel(order='F') # radmc assumes 'ij' indexing for some reason
        data.tofile(f, sep='\n', format="%13.6e")
        

def get_molecule_info(names=['co'],outdir=''):
    leiden_url = 'https://home.strw.leidenuniv.nl/~moldata/datafiles/'
    co_molecules = dict(zip(['co','13co','c17o','c18o'],['co','13co','c17o','c18o']))
    shock_tracers = dict(zip(['sio-hot','sio','so'],['sio-h2-highT','sio-h2','so@lique']))
    if os.getcwd() != outdir:
        os.chdir(outdir)
    for name in names:
        if name in co_molecules.keys():
            molecule = co_molecules[name]
        elif name in shock_tracers.keys():
            molecule = shock_tracers[name]
        else:
            molecule = name
        molecule_url = leiden_url + molecule + '.dat'
        molecule_file = 'molecule_'+molecule + '.inp'
        if os.path.exists(molecule_file) != True:
            if os.uname().sysname.lower().startswith('darwin'):
                os.system("curl {} -o {}".format(molecule_url, molecule_file))
            else:
                os.system("wget -0 {} {}".format(molecule_file, molecule_url))
        else:
            print('molecule file: {} exists in {}'.format(molecule_file,outdir))
            
def write_molecule_density(model,names=['co'],abundances=[1e-4]):
    co_molecules = dict(zip(['co','13co','c17o','c18o'],['co','13co','c17o','c18o']))
    shock_tracers = dict(zip(['sio-hot','sio','so'],['sio-h2-highT','sio-h2','so@lique']))
    for name, X in zip(names,abundances):
        if name in co_molecules.keys():
            molecule = co_molecules[name]
        elif name in shock_tracers.keys():
            molecule = shock_tracers[name]
        else:
            molecule = name
        if isinstance(X,str):
            print('This is where putting in an abundance file should go')
            # X = np.load(.... )
        else:
            molecule_field = (model.rho_embedded(fluid=0)/(mu*mh)) * X
        Nr = np.prod(np.array(model.grid['N']))
        with open(model.outdir+'numberdens_' + molecule + '.inp','w+') as f:
            f.write('1\n')                   # Format number
            f.write('%d\n'%(Nr))             # Number of cells
            data = molecule_field.swapaxes(0,1).ravel(order='F') # radmc assumes 'ij' indexing for some reason
           # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')