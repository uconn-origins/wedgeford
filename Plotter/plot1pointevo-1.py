#!/usr/bin/env python
# to run: python plot_disk.py &

import os
import matplotlib
import matplotlib.pyplot as plt
import optparse, time
import numpy as np
from scipy import *
from matplotlib.mlab import griddata
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
import pdb as pdb
import newcolormaps
import sys

rc('text',usetex=True)

############################################################
def readin_rho(photfile, colnum=5):

	radii = [] #empty array?
	radcount = -1

	row = 0     # row = zones
	col = 0     # col = radii
	gaparray = 0.0
	gaparray = genfromtxt(photfile) #load data from txt file

	nel = size(gaparray[:,0])
	ind = nonzero(gaparray[1:nel,0]-gaparray[0:nel-1,0])#ind is what is nonzero of radius[n]-radius[n-1]
	ind = ind[0]
	Nzones = ind[0] + 1 #number of z levels?
	Nradii = nel/Nzones #number of radius levels?
	gapnew = gaparray
	r = gapnew[:,0] #read in radii
	z = gapnew[:,1] #read in heights
	rho = gapnew[:,colnum] #read in value we are making contours of
	heights = np.zeros((Nzones,Nradii)) #heights is matrix of deminsions Nzones by Nradii filled with zeros
	abunsto = np.zeros((Nzones,Nradii))
	rhosto = np.zeros((Nzones,Nradii))
	radcount = -1
	col = 0

	for t in range(0,Nradii*Nzones): #Nradii*Nzones is number of elements in the matrix
		if np.mod(t,Nzones) == 0: #if Nzones goes into t evenly
			radval = gapnew[t,0] # 0th row if tth column of datafile
			radii.append(radval) #add radval to bottom of radii
			radcount += 1
			colv = gapnew[t:t+Nzones,:]
			heights[:,col] = colv[:,1]
			abunsto[:,col] = colv[:,colnum]
			rhosto[:,col] = colv[:,2]
			col += 1

	return radii, heights, abunsto, rhosto, Nradii, Nzones

############################################################

def main():

	fig = plt.figure(1)
	fig.clf()

	xloc=np.array([0,100,200,300,400])
	yloc=np.array([0,50,100,150,200,250,300])
	steps = 180

	col2plot = 5  # if 0=RA, 1=DEC, 2=surface density
	abun = np.zeros(steps)
	abun0 = np.zeros(steps)
	tstep = np.zeros(steps)
	abun2 = np.zeros(steps)
	mAU = 1.496e11 #1 AU = 1.496e11
	hnum = 44#40-1
	rnum = 66#63

	t1 = np.logspace(np.log10(1),np.log10(6e6),num=steps,endpoint=True,base=10.0)
	t2 = np.logspace(np.log10(1),np.log10(6e6),num=steps,endpoint=True,base=10.0)

	for t in range(0,steps):
    		nmf = '/data/disk_chemistry/MasterChemistry/runs/new_0.03_0.5mslgm_d200.0_Tgas_drl/e1/lime_CO/CO_time'+str(t)+'.dat'
		radii1, heights1, Sigma1, rho1, Nradii, Nzones = readin_rho(nmf,colnum=5)
		radii1 = np.vstack((radii1, radii1)) #stack arrays, so double size of radii1
		for i in range(Nzones-2):
			radii1 = np.vstack((radii1, radii1[0]))
		heights1 = np.array(heights1)
		radii1 = radii1/mAU
		heights1 = heights1/mAU
		abun[t]=Sigma1[hnum,rnum]

		tstep[t]=t




	plt.loglog(t2,abun,label="0.03")

    	plt.xlabel("time (Myr)")
    	plt.ylabel("n(CO)/n(H2)")
    	plt.legend(loc=3)
    	plt.show()




main()
