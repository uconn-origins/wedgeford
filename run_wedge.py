from pylab import *
import numpy as np
import os
from scipy import interpolate
import sys
import models as wedge

#basic script for a setup and run of wedgeford with thermal, high energy transfer, and chemical computations

outname = 'test' #your name here
nthreads = 4 #number of threads to do computations with
nphot_therm = 500000 # number of photons for thermal computation
nphot_he = 100000 # number of photons per high energy band

p0 = 'default'
outdir = '/out/'+ outname + '/'
pdef = wedge.new_model(p0)

pdef['Ms'] = 1 #stellar mass in Msun
pdef['Rs'] = 2.0 #stellar radius in Rsun
pdef['Ts'] =  4000 #stellar temp in K
pdef['xray'] = True #compute Xray-accretion model
pdef['accrate'] = 1e-07 #stellar accretion rate in Msun/yr
pdef['f'] = 0.01 # filling factor for stellar accretion 
pdef['xmodel'] =  None #default Xray model
pdef['Mdisk'] = 0.1 #disk mass total in Msun
pdef['Mfrac'] = [0.01, 0.02] # fraction of mass in dust for 1) small dust 2) large dust
pdef ['R0'] = [5, 5] #inner radius of hole in the dust for 1) small dust 2) large dust
pdef['H0'] = [1, 0.2] #scaling factor for the dust scale height compared to gas scale height 
                     # for 1) small dust 2) large dust. H0 = 1  dust H = gas H
pdef['p']= [-1, -1] # slope of surface density for 1) gas+small dust 2) large dust
pdef['Rdisk'] = [125, 125] # outer disk radius in AU (critical radius for LBP viscous disk)
pdef['Tfac'] = 1 # factor to multiply the disk temperature by compared to equilibrium temp from stellar flux
pdef['q'] = 0.5 # temperature power-law
pdef['Min'] = 1e-06 #infall rate in Msun/yr
pdef['Rc'] = 125 #centrifugal radius of infall in AU
pdef['rho_amb'] = 1e-25 #background ambient density in g/cm^3, minimum density in the box
pdef['rho_0'] = 3e-22 # base density at the bottom of the envelope cavity in g/cm^3
pdef['theta_min'] = 25 # minimum streamline angle for envelope in degrees, 0 = no cavity, 45 = half cavity, etc.
pdef['exf']= 0.25 #density drop-off power-law in the cavity
pdef['Rmax'] =  15000.0 # outer radius for Ulrich cavity parameters in AU (should be outside the grid)
pdef['d2g']= 0.01 #dust to gas ratio of envelope
pdef['shock'] = False #whether to compute with an accretion shock
pdef['nstreams'] =  1 #how many individual streams to use in azimuth, nstreams=1 is the whole envelope 
pdef['stream_frac'] =  1 #filling factor of streamlines, stream_frac = 1 is the whole envelope
pdef['N'] =  [180, 90, 48] # number of grid cells in radius, inclination, azimuth
pdef['min'] = [0.1, 0.19634954084936207, 0] # minimum grid radius, inclination, azimuth angle
pdef['max']= [400, 1.5707963267948966, 6.283185307179586] #maximum grid radius, inclination, azimuth angle 
pdef['spacing']= ['log', 'lin', 'lin'] #grid spacing for r, theta, phi log= log, lin=linear
pdef['rho_si']= 3.1518 #silicate density for dust
pdef['amin'] = [0.005, 0.005] #amin of power law dust size distribution in microns for 1) small dust 2) large dust
pdef['amax']= [1, 1000.0] #amax of power law dust size distribution in microns for 1) small dust and 2) and large dust
pdef['apow']= [3.5, 3.5] # power law of dust size distribution for 1) small dust and 2) large dust
pdef['cr_model']= 'ssx' #cosmic ray ionization model based on Cleeves 2013
pdef['G0'] =  1 #background external UV radiation in G0 
pdef['viscous_heating'] = False #whether to use internal viscous heating for radmc temperature computation
pdef['fLya'] = 0.0001 # fraction of luminosity in the lyman alpha line

m0 = wedge.initialize_model(pdef,outdir=outdir)
     
wedge.prep_thermal_transfer(m0,nphot=nphot_therm)
wedge.do_thermal_transfer(m0,nt=nthreads)

nu = 5 #number of UV pointsfor mcmono computation
nx = 8 #number of Xray points for mcmono computation
wedge.prep_he_transfer(m0,nu=nu,nx=nx)
wedge.do_he_transfer(m0,nphot=nphot_he,nt=nthreads,prep=False)

m0.calc_Tgas(ndim=3)

import models.prepchem as wc
c0 = wc.chemdisk(m0,chemdir=outname,rgrid=None,zgrid=None)
wc.plot_prechem(c0)
wc.write_chem_inputs(c0)
