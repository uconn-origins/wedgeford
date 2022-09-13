#!/usr/bin/env python

import os
import matplotlib
import matplotlib.pyplot as plt
import optparse, time
import numpy as N
import pdb as pdb
from scipy import *
from matplotlib.mlab import griddata
from matplotlib import rc
from matplotlib import rcParams
from pylab import *
from subprocess import call
import fnmatch
import shutil

##################################################################################
# Incorporate additional gas temperature and CR information into environment files.
# This code will expect that the UV file is in the main disk chemistry directory and 
# that the environment file contains a copy of the relevant gas_Model.out file. 
#
# The script can also produce 0io files and copy into all relevant directories.

# Location of disk chemical model:
diskchem_dir = '/Nirgal1/kamberrs/disk_chemistry/MasterChemistry'
#cp -r /Nirgal1/kamberrs/torus/environ/Ms2.2Md0.0* .

# Name(s) of torus environment output:
mass = ['0.03Ms']
lg = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]#,0.9,0.9,0.9,0.9,0.0,0.9,0.9,0.99,0.99,0.99]
Routmm = [200.0,200.0,200.0,200.0,200.0,200.0,200.0,200.0]#158.0,123.0,60.0,200.0,200.0,158.0,158.0,123.0,60.0]
in2d = [0,0,0,0,0,0,0,0]#,0,1,1,1,0,1,1,1,1,1]
srt_time_vals = ['1.000D-00','1.000D-00','1.000D-00','1.000D-00','1.000D-00','1.000D-00','1.000D-00','1.000D-00']#,'1.000D-00','1.000D-00']#,'1.000D-00','1.321D+06','2.151D+06','2.971D+06','1.000D-00','5.001D+05','1.321D+06','1.501D+06','2.151D+06','2.971D+06']#,'None']
end_time_vals = ['6.000D+06','6.000D+06','6.000D+06','6.000D+06','6.000D+06','6.000D+06','6.000D+06','6.000D+06']#,'6.000D+06','6.000D+06']#,'1.320D+06','2.150D+06','2.970D+06','6.000D+06','5.001D+05','1.320D+06','1.500D+06','2.146D+06','2.970D+06','6.000D+06']
flags = ['g','g','g','g','g','g','g','g']#,'d','d','d','d','dg','dg','dg','dg','dg','dg']

# lg = [0.0,0.9,0.9,0.99,0.99,0.99]
# Routmm = [200.0,200.0,158.0,158.0,123.0,60.0]
# in2d = [0,1,1,1,1,1]
# srt_time_vals = ['1.000D-00','5.001D+05','1.321D+06','1.501D+06','2.151D+06','2.971D+06']#,'None']
# end_time_vals = ['5.001D+05','1.320D+06','1.500D+06','2.146D+06','2.970D+06','6.000D+06']
# flags = ['dg','dg','dg','dg','dg','dg']

for m in range(len(mass)):
	for l in range(len(lg)):
		dest = diskchem_dir+"/store_inp/"+mass[m]+"_"+str(lg[l])+"_"+str(Routmm[l])+"_"+flags[l]+"/"
		oldpath = os.path.exists(dest)
		if oldpath == False:
			os.makedirs(dest)
		f = open(dest+"/2times.inp","w+")
		f.write("# Time steps:\n")
		f.write(end_time_vals[l]+"              #last time step (yr),\n")
		f.write(srt_time_vals[l]+"              #first time step (yr),\n")
		f.write("180                     #total amount of time steps\n")
		f.close()
		f2 = open(dest+"/5flags.inp","w+")
		f2.write("###############################################\n")
		f2.write("# Parameters\n")
		f2.write("#################################################\n")
                f2.write("# fg (#)      = dust / gas ratio in disk\n")
		f2.write("# freezeeffic = Computed by freezeeffic = (r / 0.1 mu)**-1.5.  Simulates net area loss due to grain growth.\n")
		f2.write("fg            = 1.00D+2\n")
		f2.write("freezeeffic   = 1.0\n")
		f2.write("epsilon       = 1.0\n")
		f2.write("#################################################\n")
		f2.write("# Flags\n")
		f2.write("#################################################\n")
		f2.write("# 0 = off, anything else = on\n")
		f2.write("# Crdesorp       = CR desorption\n")
		f2.write("# CRionization   = CR ionization\n")
		f2.write("# photodesorp    = Photodesorption\n")
		f2.write("# LyAphotodesorp = Lyman Alpha photodesorption\n")
		f2.write("# thermaldesorp  = Thermal Desorption / Classical Evaporation\n")
		f2.write("# include_lya    = if false (0), remove the lyman alpha flux from the UV field\n")
		f2.write("# xraydust       = dust-dependent x-ray opacity\n")
		f2.write("# incl_radionuc  = turn on/off radionuclide ionization\n")
		f2.write("CRdesorp        = 1\n")
		f2.write("CRionization    = 1\n")
		f2.write("photodesorp     = 1\n")
		f2.write("LyAphotodesorp  = 1\n")
		f2.write("thermaldesorp   = 1\n")
		f2.write("include_lya     = 1\n")
		f2.write("xraydust        = 1\n")
		f2.write("incl_radionuc   = 0\n")
		f2.write("incl_isrf       = 1\n")
		f2.write("incl_2dabun     = "+str(in2d[l])+"\n")
		f2.write("write_2dabun    = 1\n")
		f2.write("spatial_dust    = 1\n")
		f2.write("#################################################\n")
		f2.write("# Testing\n")
		f2.write("#################################################\n")
		f2.write("# ratetest        = print out file with formation/destruction rates.  0 = off\n")
		f2.write("# testspec (name) = species to print out formation/destruction rate of\n")
		f2.write("#       (max=25, 160 characters)\n")
		f2.write("# shieldtest      = print out file with self shielding test data.  0 = off\n")
		f2.write("ratetest       = 1\n")
		f2.write("testspec       = CO,C2H,H2O,N2,N2H+,C3H2,HCO+,CN,HCN,NH3,H2CO,NO,He+,H3+,C+,O+\n")
		f2.write("shieldtest    = 1\n")
		f2.close()

		


