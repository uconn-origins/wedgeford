from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata

from . units  import *
from . radiation import *
from . dust import *

"""
Functions meant for writing the input files for radmc3d

"""

def _scalarfieldWriter(model, data=None, fname=''):
    """Writes a scalar field to a file.

        Parameters
        ----------

        data   : ndarray
                Scalar variable to be written

        fname  : str
                Name of the file containing a scalar variable

        """
    nx,ny,nz = model.grid['N']
    with open(fname, 'w') as wfile:
        if len(data.shape) == 3:
            hdr = np.array([1, nx * ny * nz], dtype=int)
            hdr.tofile(wfile, sep=" ", format="%d\n")
            #need fortran indexing
            data = data.ravel(order='F')
            data.tofile(wfile, sep=" ", format="%.9e\n")

        elif len(data.shape) == 4:
            hdr = np.array([1, nx * ny * nz, data.shape[3]], dtype=int)
            hdr.tofile(wfile, sep=" ", format="%d\n")
            #need fortran indexing
            data = data.ravel(order='F')
            data.tofile(wfile, sep=" ", format="%.9e\n")
        else:
            raise ValueError('Incorrect shape. Data stored in a regular grid should have 3 or 4 dimensions'
                                         + ' with the fourth dimension being the dust species. '
                                         + ' The data to be written has a dimension of ' + ("%d" % len(data.shape))
                                         + '\n No data has been written')
            
            
def write_wavelength(model,wav = None, fname=''):
    """Writes the wavelength grid to a file (e.g. wavelength_micron.inp).

    Parameters
        ----------
    wav: optional
        wavelengths to write to file, if omitted, a template wavelength will be read in from the templates directory
        
    fname  : str, optional
        File name into which the wavelength grid should be written. If omitted 'wavelength_micron.inp'
                 will be used
    """
    if wav is None:
        templates_dir = model.models_dir+'templates/'
        wav,freq = read_wavelength(fname=templates_dir+'wavelength_micron.inp')
        if model.rad['xray'] == True:
            xwav = np.logspace(np.log10(xray_min),np.log10(xray_max),20)
            wav = np.append(xwav,wav[wav>xray_max])
    
    if fname == '':
        fname = model.outdir+ 'wavelength_micron.inp'
    nwav = len(wav)
    print('Writing ' + fname)
    
    with open(fname, 'w') as wfile:
        wfile.write('%d\n' % nwav)
        for ilam in range(nwav):
            wfile.write('%.9e\n' % wav[ilam])
    return 
    

def write_grid(model):
    """ writes the spatial grid file for radmc3d, a 3d spherical grid with mirror symmetry
    Parameters
    ------------
    model: model class object to write grid information from 

    """
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
    print('writing amr_grid.inp')
    with open(model.outdir+"amr_grid.inp","w") as f:
        f.write(header)
        for x,fmt  in zip([r_edges,th_edges,phi_edges],['%13.6e','%17.10e','%13.6e']):
            x.tofile(f, sep= '\t', format=fmt)
            f.write('\n')
    return
            
            
def write_dust_density(model,envelope=True,rhodust = None):
    """ writes the dust density of small and large dust components in a dust_density.inp file
    for radmc3d

    Parameters
    ----------

    model: model class object to write densities from

    envelope: boolean, if True then write the full density distribution of envelope + disk
            if False, writes the densities of the disk only
            
    rhodust: if None, densities from model, otherwise, custom array will be written
    """
    if rhodust is None:
        if envelope == True:
            small_dust = model.rho_embedded(fluid=1)
            large_dust = model.rho_embedded(fluid=2)
        else:
            small_dust = model.rho_disk(fluid=1)
            large_dust = model.rho_disk(fluid=2)
        rhodust = np.stack((small_dust,large_dust),axis=-1)

    print('writing dust_density.inp')
    fname = model.outdir + 'dust_density.inp'
    _scalarfieldWriter(model,data=rhodust.swapaxes(0,1), fname=fname)
    return
    
def write_gas_density(model,envelope=True,rhogas=None):
    """ writes the gas density in a gas_density.inp file g/cm^3
    for radmc3d

    Parameters
    ----------

    model: model class object to write densities from

    envelope: boolean, if True then write the full density distribution of envelope + disk
            if False, writes the densities of the disk only
            
    rhogas: if None, densities from model, otherwise, custom array will be written
    """
    
    if rhogas is None:
        if envelope == True:
            rhogas = model.rho_embedded(fluid=0).swapaxes(0,1)
        else:
            rhogas = model.rho_disk(fluid=0).swapaxes(0,1)
            
    print('writing gas_density.inp')
    fname = model.outdir + 'gas_density.inp'
    _scalarfieldWriter(model,data=rhogas, fname=fname)
    return
    
def write_star(model,wav=None,accrate=None,fnu=None):
    """ writes the stars.inp file for radmc3d
    summing up three components: blackbody emission, uv continuum from accretion onto central star,
    and (optional) xray free-free emission from accretion onto central star
    
    Parameters
    ----------
    model: model class object, uses model.star parameters to write file
    
    wav: wavelength in microns, if None, reads from local wavelength_micron.inp file
    
    accrate: accretion rate in Msun/yr, if None, uses model accretion rate
    
    fnu: flux in ergs/cm^2/s/Hz, if None, adds up contributions from component models, otherwise custom flux can be entered
    
    """
    Rstar = model.star['Rs']*Rsun
    Mstar = model.star['Ms']*Msun
    
    if wav is None:
        wav, freq = read_wavelength(fname=model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
    
    if accrate is None:
        accrate = model.star['accrate']
        
    nwav = len(wav)
    
    if fnu is None:
        fnu_tot = calc_input_spectrum(model,wav=wav)
    else:
        if len(fnu) == nwav:
            fnut_tot = fnu
        else:
            print('Error: input spectrum does not match wavelength grid')
            return None
            
    print('Writing stars.inp')
    with open(model.outdir + 'stars.inp', 'w') as wfile:
        wfile.write('%d\n' % 2)
        wfile.write('%d %d\n' % (1, nwav))
        wfile.write('%.9e %.9e %.9e %.9e %.9e\n' % (Rstar, Mstar,0,0,0))
        wfile.write('%s\n' % ' ')
        for ilam in range(nwav):
            wfile.write('%.9e\n' % wav[ilam])
        wfile.write('%s\n' % ' ')
        
        for ilam in range(nwav):
            wfile.write('%.9e\n' % (fnu_tot[ilam]))
    return 0
    

def write_external_radfield(model,wav=None,fnu=None):
    """ writes the external_source.inp file for radmc3d
    Parameters
    ----------
    model: model class object, uses model.star parameters to write file
    
    wav: wavelength in microns, if None, reads from local wavelength_micron.inp file
    
    fnu: flux in ergs/cm^2/s/Hz/str, if None, calculates ISRF from template, otherwise custom intensity can be entered
    
    """
    if wav is None:
        wav, freq = read_wavelength(model.outdir+'wavelength_micron.inp')
    else:
        freq = c / wav * 1e4
    
    nwav = len(wav)
    
    if fnu is None:
        wav, fnu_field = calc_ISRF(model,wav=wav)
    else:
        if len(fnu) == nwav:
            fnut_field = fnu
        else:
            print('Error: input spectrum does not match wavelength grid')
            return None
        
    print('writing external_source.inp')
    with open(model.outdir+'external_source.inp','w') as f:
        f.write('2 \n')
        f.write(str(len(wav))+ '\n')
        f.write('\n')
        wav.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')
        f.write('\n')
        fnu_field.tofile(f, sep='\n', format="%13.6e")
    return 0
        

def write_viscous_heatsource(model,accrate=None):
    """ writes the heatsource.inp file for radmc3d
    Parameters:
    -----------
    model: model class object, uses model.star parameters to write file
    
    accrate: accretion rate in Msun/yr, if None, uses model accretion rate
    """
    if accrate is None:
        accrate = model.star['accrate']*Msun/yr
    else:
        accrate = accrate*Msun/yr
        
    D_disk = model_viscous_dissipation(model,accrate=accrate)
    Nr = np.prod(np.array(model.grid['N']))
    if write == True:
        with open(model.outdir+'heatsource.inp','w') as f:
            f.write('1\n')                   # Format number
            f.write('%d\n'%(Nr))             # Number of cells
            heat = D_disk.swapaxes(0,1).ravel(order='F')# radmc assumes 'ij' indexing for some reason Create a 1-D view, fortran-style indexing
            heat.tofile(f, sep='\n', format="%13.6e")
    return

def write_opacities(model,ndust=2,filenames=['',''],update=True):
    """ checks for the dust_kappa..inp files, runs optool calculation if not present to write new ones
    updates the dust_opac.inp file with the names of the dust_kappa.inp files
    Parameters:
    -----------
    model: model_ class object
    
    ndust: number of dust populations to include in dust_opac.inp
    
    filenames: list of str, optional, names of the dust opacity files 
    
    update: boolean, if True, opacities are updated with the x-ray opacities from BB11
    
    """
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
            print('{}: user dust opacity file not found in model directory'.format(kappa_file))
        elif os.path.exists(kappa_file) != True and fname == '':
            print('Running optool to generate new opacities')
            run_optool(model,fluid=fluid,na=60)
            exts[str(fluid)] = 'dust-{}'.format(fluid)
    
        if update == True:
            print('Updating x-ray opacities')
            kappa_file = "{}dustkappa_dust-{}x.inp".format(model.outdir,fluid)
            #if os.path.exists(kappa_file) != True:
            wav, kabs, kscat, g  = read_kappa(model,fluid=fluid,filename=fname)
            wav, kabs_new = model_kappa_xray(model,wav=wav,fluid=fluid)
            kabs[kabs_new > 0] = kabs_new[kabs_new > 0]
            arrays = np.stack((wav,kabs,kscat,g),axis=-1)
            nwav = len(wav)
            header = '{} \n{} \n'.format(iformat,nwav)
            np.savetxt(kappa_file,arrays,header=header,fmt='%12.6e',comments='',delimiter='\t')
            exts[str(fluid)] = 'dust-{}x'.format(fluid)
    print('updating dustopac.inp')      
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
            
def write_main(model,nphot= 100000,scat=2,mrw =1,maxtau=15,teq=1,dust=1,lines=0,coup=0):   
    """ writes radmc3d.inp file with desired parameters
    Parameters
    ----------
    model: model class object, model.outdir stores where to write the file
    
    all others:
    see radm3d.inp input documentation
    """
    print('updating radmc3d.inp')
    with open(model.outdir+'radmc3d.inp','w') as f:
        f.write('nphot = {}\n'.format(nphot))
        f.write('scattering_mode_max = {}\n'.format(scat))
        f.write('iranfreqmode = 1\n')
        f.write('modified_random_walk = {}\n'.format(mrw))
        f.write('mc_scat_maxtauabs = {}\n'.format(maxtau))
        f.write('tgas_eq_tdust = {}\n'.format(teq))
        f.write('incl_dust = {}\n'.format(dust))
        f.write('incl_lines = {}\n'.format(lines))
        f.write('itempdecoup = {}\n'.format(coup))
    return True
            
##### for line transfer ##########################################################
def write_lines(model,names=['co']):
    """ writes the lines.inp file for radmc3d in leiden format
    Parameters
    ---------
    model: model class object, model.outdir stores where to write the file
    
    names: list of str, names of molecules to include in line transfer
    """
    with open(model.outdir+'lines.inp','w+') as f:
        f.write('2\n')                   # Format number
        f.write('{}\n'.format(len(names)))             # Number of molecules
        for name in names:
            f.write('{} leiden 0 0 0\n'.format(name))
    return True

def write_velocities(model,envelope = True):
    """ writes the gas_velocity.inp file for radmc3d
    Parameters
    ----------
    model: model class object, model.outdir stores where to write the file
    
    envlope: boolean, if True, envelope is included in the calculation of the velocity field, 
                      if False, disk only
    """
    if envelope == True:
        velocities = model.v_embedded()
    else:
        velocities = model.v_disk()
    Nr = np.prod(np.array(model.grid['N']))
    print('writing gas_velocity.inp')
    with open(model.outdir+'gas_velocity.inp','w+') as f:
        f.write('1\n')                   # Format number
        f.write('%d\n'%(Nr))             # Number of cells
        vr = velocities['r'].swapaxes(0,1).ravel(order='F')
        vth = velocities['theta'].swapaxes(0,1).ravel(order='F')
        vphi = velocities['phi'].swapaxes(0,1).ravel(order='F')
        for i,j,k in zip(vr,vth,vphi):
            f.write("%9e %9e %9e\n" % (i, j, k))
    return
            
def write_gas_temperature(model,Tgas=None):
    """ writes the gas temperature to a file 
    Parameters
    ----------
    model: model class object, model.outdir stores where to write the file
    
    Tgas: 3d ndarray of the gas temperature 
    """
    print('writing gas_temperature.inp')
    fname = model.outdir+'gas_temperature.inp'
    _scalarfieldWriter(model,Tgas,fname=fname)
    return

def get_molecule_info(model,names=['co']):
    """ downloads molecule files from leiden database, unless files already exist locally
    Parameters:
    ----------
    model: model class object, model.outdir stores where to write the file
    
    names: list of str, molecules to download from leiden database, names are from data file names
    
    https://home.strw.leidenuniv.nl/~moldata/datafiles/
    
    """
    leiden_url = 'https://home.strw.leidenuniv.nl/~moldata/datafiles/'
    co_molecules = dict(zip(['co','13co','c17o','c18o'],['co','13co','c17o','c18o']))
    shock_tracers = dict(zip(['sio-hot','sio','so'],['sio-h2-highT','sio-h2','so@lique']))
    for name in names:
        if name in co_molecules.keys():
            molecule = co_molecules[name]
        elif name in shock_tracers.keys():
            molecule = shock_tracers[name]
        else:
            molecule = name
        molecule_url = leiden_url + molecule + '.dat'
        molecule_file = 'molecule_'+molecule + '.inp'
        if os.path.exists(model.outdir+molecule_file) != True:
            if os.uname().sysname.lower().startswith('darwin'):
                os.system("curl {} -o {}".format(molecule_url, molecule_file))
            else:
                os.system("wget {} {}".format(molecule_file, molecule_url))
        else:
            print('molecule file: {} exists in {}'.format(molecule_file,model.outdir))
    return 
            
def write_molecule_density(model,names=['co'],abundances=[1e-4]):
    """ writes the number density of molecules into a numberdens_X.inp file for radmc3d line transfer
    where X is the molecule name
    
    Parameters:
    ----------
    model: model class object, model.outdir stores where to write the file
    
    names: list of str, molecules from leiden database, file extension name
    
    abundances: list of scalar or list of ndarray
                if scalar, then the factor to multiple the number density of H2 by 
                if ndarray, then an array of the number density of the molecule species
    
    
    """
    co_molecules = dict(zip(['co','13co','c17o','c18o'],['co','13co','c17o','c18o']))
    shock_tracers = dict(zip(['sio-hot','sio','so'],['sio-h2-highT','sio-h2','so@lique']))
    
    for name, X in zip(names,abundances):
        if name in co_molecules.keys():
            molecule = co_molecules[name]
        elif name in shock_tracers.keys():
            molecule = shock_tracers[name]
        else:
            molecule = name
        if X.ndim > 1:
            molecule_field = X #custom molecule field instead
        else:
            molecule_field = (model.rho_embedded(fluid=0)/(mu*mh)) * X
            molecule_field = molecule_field.swapaxes(0,1)
        fname = model.outdir+'numberdens_' + molecule + '.inp'
        print('writing',fname)
        _scalarfieldWriter(model,molecule_field, fname=fname)
    return
