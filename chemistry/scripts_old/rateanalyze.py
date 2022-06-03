# Code to compare two chemical models
# Column-wise, or single wise

import os
import fnmatch
import matplotlib
import matplotlib.pyplot as plt
import optparse, time
import numpy as N
import pdb as pdb
from scipy import *
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
from matplotlib import rc
from matplotlib import rcParams
from pylab import *
from subprocess import call
import re
import sys

def reacanalysis(routfile,reacfile,time,zone,mnrt):

    f = open(routfile)
    contents = f.readlines()
    f.close()
    nl = len(contents)
    timesv = N.zeros(500)-100.0
    inds = list(range(500))

    ct = 0
    for i in range(nl-3):
        line = contents[i+3].split()
                #print 'line:'
                #print line
        if N.float(line[0]) == zone:
            th = contents[i+3].split()
            timesv[ct] = N.float(th[1])
            inds[ct] = i+3
            ct = ct + 1
        elif N.float(line[0]) > zone:
            break
    timesv = timesv[0:ct]+0.0
    inds = inds[0:ct]

    gm = N.argmin(N.abs(timesv-time))
    reacdata = contents[inds[gm]].split()
    reacdata = reacdata[2:]
    lndata = len(reacdata)
    #print reacdata
    rshp = N.reshape(reacdata,(int(lndata/2),2))

    id = list(range(int(lndata/2)))
    rat = N.zeros(int(lndata/2))

    for i in range(int(lndata/2)):
        id[i] = N.int(rshp[i,0])
        rat[i] = N.float(rshp[i,1])
        #print id,rat

    f = open(reacfile)
    reacall = f.readlines()
    f.close()
    nl = len(reacall)

    ct = 0
    idfind = id[ct]
    idfindold = ''

    reactants = []
    for i in range(nl):
        line = reacall[i].split()
        idl = line[0]
        if idl[0] == '#':
            continue

        if idfind == idfindold:  # catch doubles
            for b in range(6):
                # check if string or float
                if re.match("^[A-Z]", line[b]):  # Checks if the first letter is a capital letter == Chemical.
                    #print "Counted:", idfind,line[b]
                    reactants.append(line[b])
                else:
                    #print 'didnt count:'+line[b]
                    reactants.append('')
            ct = ct+1 # Find next specie
            if ct < lndata/2:
                idfind = id[ct]
            else:
                break

        if N.int(idl) == idfind:
            for b in range(6):
                # check if string or float
                if re.match("^[A-Z]", line[b+1]):  # Checks if the first letter is a capital letter == Chemical.
                    #print "Counted:", idfind,line[b+1]
                    reactants.append(line[b+1])
                else:
                    #print 'didnt count:'+line[b+1]
                    reactants.append('')

                #print reactants[nl*6:nl+6]
            ct = ct+1 # Find next specie
            #print ct
            idfindold = idfind
            #print idfindold
            if ct < lndata/2:
                idfind = id[ct]
            else:
                break


    strt = N.argsort(N.abs(rat))
    strt = strt[::-1]

    print( 'Printing rate results for t = {0} Myr and zone = {1}.'.format(time/1e6,zone))
    ntrt = 0.0
    for i in range(int(lndata/2)):
        print( reactants[strt[i]*6:(1+strt[i])*6])
        if N.abs(rat[strt[i]]) > mnrt:
            print( "Rate: {0:10.3e}".format(rat[strt[i]]))
            #print "Rate: {0:10.3e}   R1,R2: {1:10} + {2:10} -->    {3:10} + {4:10} + {5:10} + {6:10} ".format(rat[strt[i]],reactants[strt[i]*6],reactants[(strt[i])*6+1],reactants[(strt[i])*6+2],reactants[(strt[i])*6+3],reactants[(strt[i])*6+4],reactants[(strt[i])*6+5])
        ntrt = rat[strt[i]] + ntrt
    print( "Net Rate: {0:10.3e}".format(ntrt))

    return


def main():
#0.1Mslg_0.6_w98_x30_0.03cr
#6Myr_new_0.03Mslg_0.6m_d200.0_Tgas_w98_1g0_x30
    # First, tell me what output file (i.e. molecule) you want to analyze.  These have extensions '*.rout'.
    routfile = '/data/disk_chemistry/MasterChemistry/runs/rr_0.003Mslg_0.5m_d200.0_Tgas_drl/e1/CO2(gr)_e1_15.5866.rout'
    # Second, point me to the right reaction file -- i.e. the one that you used to run the chemistry (just in case you have different ones).
    reacfile = '/home/kschwarz/MasterChemistry/rreacs_herb0308_CDR_FullGr.dat'

    # Now tell me specifically where you want to look, in time and vertical position:
    time =  1e4
#1.048E+03  # Time step in years.
    zone = 44   # Z-value of zone/height you want to analyze.
    minrt = 1e-30   # Set the minimum rate you care to print out here.  You can play with this.

    # Magic happens!
    reacanalysis(routfile,reacfile,time,zone,minrt)
    return

main()
