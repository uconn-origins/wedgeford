#Converts radmc3d model input/output into format needed for UV/X-ray radiative transfer and chemistry
import numpy as np
import os
from pylab import *

Msun = 1.9891e33 # g
Rsun = 69.6e9 #cm
Lsun  = 3.8525e33 #erg/s
AU = 1.49598e13 #cm
yr = 3.14e7 #seconds
cmtokm = 10**(-5) #converts cm to km
R_mu = 36149835 # (R/mu) for 2.4 g/mol in PPD
Gconv = 6.6743e-8 # cgs
S0conv = (Msun/(AU**2)) #S0cgs = S0code*S0conv
sigsb = 5.67e-5 #Stefan Boltzmann constant
c = 2.99792458e10 # speed of light
mu = 2.36
h = 6.6260755e-27 #planck
mh = 1.67e-24 # mass of hydrogen cgs
Habing = 1.83590e+08  # Photons/cm^2/s between 930-2000
G0 = 2.7e-3 # erg/cm^2/s^-1



def write_out(model,outname='model',ndust=2): #write the .out file for RT and chemistry, pass the model class into this
    m = model
    if os.getcwd() != m.outdir:
        os.chdir(m.outdir)
    r,z = m.make_rz_H()
    T = m.get_Tdust()
    rho_g = m.rho_embedded(fluid=0)
    rho_d = np.zeros_like(T)
    for fluid in np.arange(ndust):
        rho_d += m.get_rhodust(fluid=fluid+1)
    rhod_2d = m.make_quadrant(rho_d, fill_value = np.amin(rho_d),order='F').T.flatten()
    T_2d = m.make_quadrant(T,fill_value = np.amin(T[T>1]),order='F').T.flatten()
    rhog_2d = m.make_quadrant(rho_g, fill_value = np.amin(rho_g), order='C').T.flatten()
    
    header1 = 'R(AU)    z(AU)   rhogas(g/cm3)  rhodust(g/cm^3)    T(K)   [    4 types of dust]      f_H\n'
    header2 = '-ForChem+Plots--------------------------------------------------------\n'
    outfile = open(outname+'.out',"w")
    outfile.write(header1)
    outfile.write(header2)
    for n in np.arange(len(r.flatten())):
        outfile.write('%10.4f %10.4f   %5.2e   %5.2e %7.1f %12.7f %12.7f %12.7f %12.7f   %7.2e\n'%(r.T.flatten()[n],z.T.flatten()[n],rhog_2d[n], rhod_2d[n] ,T_2d[n],1.0,0.0,0.0,0.0,1.0))
    outfile.close()
    return outname + '.out'

def read_out(outname='model.out'): #convenience function for reading the .out file and returning the bits accessible as dictionary
    data = np.loadtxt(outname,skiprows=2,usecols=(0,1,2,3,4))
    r = data[:,0]
    z = data[:,1]
    shape = (len(np.unique(r)),50)
    rhog = data[:,2].reshape(shape)
    rhod = data[:,3].reshape(shape)
    T = data[:,4].reshape(shape)
    out_dict = {'r':r.reshape(shape), 'z':z.reshape(shape),'rhog':rhog,'rhod':rhod,'Td':T}
    return out_dict

def plot_out(prop='rhog',outname='model.out',log=True,method=contourf,**pk): #plotting convenience function to inspect .out file contents
    data = read_out(outname)
    if log == True:
        method(data['r'],data['z'],np.log10(data[prop]),**pk)
    else:
        method(data['r'],data['z'],(data[prop]),**pk) 








