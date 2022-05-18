#Converts radmc3d model input/output into format needed for UV/X-ray radiative transfer and chemistry
import numpy as np
import os
from pylab import *
from models.units import *

#need to rewrite to use radmc3d output class****

def make_rz_H(model): #make a logarithmically spaced cylindrical grid
        zmin = 1e-10
        new_r = np.logspace(np.log10(model.grid['min'][0]), np.log10(model.grid['max'][0]), model.grid['N'][0])
        zf_norm = np.append(zmin, np.logspace(-4,0,50)) #faces of the z-cells
        zc_norm = 0.5*(zf_norm[1:] + zf_norm[:-1]) #50 points in the Z direction
        if return_faces == True:
            R,Z = np.meshgrid(new_r,zf_norm) #returns the bin edges at the radius, for summing up columns
        else:
            R,Z = np.meshgrid(new_r,zc_norm)
        Z *= R/np.tan(np.radians(self.m.env['theta_min']+15))
        return R,Z
    
def make_quadrant(model,quantity_3d,fill_value=0,order='F'): # to be used with the radmc 3d values
        quantity_2d = quantity_3d[:,:,0] # phi=0 plane
        r_cyl,z_cyl = model.make_rz()
        r_new, z_new = model.make_rz_H()
        quantity_2d_interp = griddata((r_cyl[:,:,0].flatten(),z_cyl[:,:,0].flatten()), quantity_2d.ravel(order=order), (r_new,z_new),fill_value=fill_value,method='linear')
        return quantity_2d_interp


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








