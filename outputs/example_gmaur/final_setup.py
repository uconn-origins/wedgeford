#radmc3d mctherm
#radmc3d sed incl 53.21 phi 57.17
#
# Import NumPy for array handling
#
import numpy as np
import sys
import astropy.io.fits as pyfits
from scipy import interpolate
#
# Import plotting libraries (start Python with ipython --matplotlib)
#
#from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt

def get_Sig(mass,r,rc,gamma):
    sigc = mass/(2.*np.pi*np.trapz(r*(r/rc)**gamma*np.exp(-1*(r/rc)**(2+gamma)),r))
    return sigc

def get_mass(Sig,r):
    mass = 2.*np.pi*np.trapz(Sig*r,r)
    return mass

#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
mj  = 1.899e30       # Jupiter mass             [g]
ts  = 3.9e3        # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity, not used        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]

lgm18 = 0.4164262351428959*mj
#sgm18 = 0.6749021502820188*mj
#
# Monte Carlo parameters
#
nphot    = 100000
#
# Grid parameters
#
nr       = 600
ntheta   = 100
nphi     = 1
rin      = 1.0*au
rin_out = 1.80*au
rim = 1.0*au
rout = 650.0*au#300.0*au
thetaup  = np.pi*0.5 - 0.7e0
#
# Disk parameters
rcrit       = 111.0*au
rref = 100*au
massdsm =1.03e-4*ms#
print("total mass small grains:",massdsm/ms)
plsig    = -1.            # Powerlaw of the surface density
plsig_d  = -1.
hr0      = 0.12#0.08#17.4/176.              # H_p/r at rref AU
plhc      = 1.35 # Powerlaw of flaring
plh      = 1.35 
chi      =3.75/7.5
plh_lg =1.35#1.0 

#
# Star parameters
#
mstar    = ms*1.1
rstar    = rs*1.9
tstar    = 4350.0 
pstar    = np.array([0.,0.,0.])
#
# Make the coordinates
#
ri       = np.logspace(np.log10(rin),np.log10(rout),nr+1)
thetai   = np.linspace(thetaup,0.5e0*np.pi,ntheta+1)
phii     = np.linspace(0.e0,np.pi*2.e0,nphi+1)
rc       = 0.5 * ( ri[0:nr] + ri[1:nr+1] )
thetac   = 0.5 * ( thetai[0:ntheta] + thetai[1:ntheta+1] )
phic     = 0.5 * ( phii[0:nphi] + phii[1:nphi+1] )

#get sigma
#sigmad0 = get_Sig(massd,np.logspace(np.log10(rin),np.log10(rout_mm),nr+1),rcrit,plsig_d)
od = np.where(ri>=rim)
cav = np.where(ri<=rin_out)
masscav = 4.0e-9*ms
print("total mass cavity dust:",masscav/ms)
plcav =plsig
sigmasm0 = get_Sig(massdsm,ri[od],rcrit,plsig)
sigcav0 = get_Sig(masscav,ri[cav],rcrit,plcav)

#
#
# Make the grid
#
qq       = np.meshgrid(rc,thetac,phic,indexing='ij')
rr       = qq[0]
tt       = qq[1]
zr       = np.pi/2.e0 - qq[1]
#
# Make the dust density model

#new
sigin = np.loadtxt('/home/kschwarz/rac-2d-master/inp/GM_Aur_sm_designed_1.txt',skiprows=1)
f = interpolate.interp1d(sigin[:,0]*au,sigin[:,1],kind='linear',bounds_error=False,fill_value=0.0)
sigmasm = f(rr)

hhr = hr0*(rr/rref)**plh
hh = hhr*rr
rhosm = ( sigmasm / (np.sqrt(2.e0*np.pi)*hh) ) * np.exp(-(zr**2/hhr**2)/2.e0)

hr0cav = hr0#0,15
hhr = hr0cav*(rr/rref)**plhc
hh = hhr*rr
sigcav = sigcav0*(rr/rcrit)**plcav*np.exp(-1.*(rr/rcrit)**(2.+plcav))
sigcav[np.where(rr>rin_out)] = 0.
rhocav = ( sigcav / (np.sqrt(2.e0*np.pi)*hh) ) * np.exp(-(zr**2/hhr**2)/2.e0)

sigin = np.loadtxt('/home/kschwarz/rac-2d-master/inp/GMAur_sigma_lg.txt',skiprows=1)
f = interpolate.interp1d(sigin[:,0]*au,sigin[:,1],kind='linear',bounds_error=False,fill_value=0.0)
sigmalg = f(rr)
hhr = chi*hr0*(rr/rref)**plh_lg
hh = hhr*rr
rholg = ( sigmalg / (np.sqrt(2.e0*np.pi)*hh) ) * np.exp(-(zr**2/hhr**2)/2.e0)

#
# Write the wavelength_micron.inp file
#
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
nlam     = lam.size
#
# Write the wavelength file
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    for value in lam:
        f.write('%13.6e\n'%(value))
#
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    for value in lam:
        f.write('%13.6e\n'%(value))
    f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system: spherical
    f.write('0\n')                       # gridinfo
    f.write('1 1 0\n')                   # Include r,theta coordinates
    f.write('%d %d %d\n'%(nr,ntheta,1))  # Size of grid
    for value in ri:
        f.write('%13.6e\n'%(value))      # X coordinates (cell walls)
    for value in thetai:
        f.write('%13.6e\n'%(value))      # Y coordinates (cell walls)
    for value in phii:
        f.write('%13.6e\n'%(value))      # Z coordinates (cell walls)
#radmc3d mctherm
#radmc3d sed incl 53.21 phi 57.17
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nr*ntheta*nphi))     # Nr of cells
    f.write('3\n')                       # Nr of dust species
    data = rholg.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
    data = rhosm.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
    data = rhocav.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    data.tofile(f, sep='\n', format="%13.6e")
    f.write('\n')
#
# Dust opacity control file
#
with open('dustopac.inp','w+') as f:
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
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('smdsharp        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 1\n')
    f.write('iranfreqmode = 1\n')
