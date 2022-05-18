#!/usr/bin/env python

from math import *
import os
import pdb as pdb
import re
import numpy as N
# renumber rreacs file (in case reactions added, subtracted, commented out or moved)

def all_indices(value, qlist):
    indices = []
    idx = -1
    while True:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices

def main(name_reac):

    # Read in entire new rate file with : delimited.
    inp = '/Shire1/cleeves/Research/PhotoDissUMIST2012/RATE12_binding_energies.dist.txt'
    
    fumist = open(inp, 'r')
    dat = fumist.readlines()
    fumist.close()

    spec    = []
    BE      = []
    
    ctt = 0
    
    # Read in lines from file data, dat:
    for line in dat:
        spln = line.split()
        if spln[0] != '#':
            spec.append(spln[0])
            BE.append(spln[1])            

    fin = open(name_reac+'.dat', 'r')
    fout = open(name_reac+'_NewBE.dat', 'w')

    nreac = 1
    for line in fin:
        if line[0] == "#":
            fout.write(line)
            continue
        outline = " %4d %s" % (nreac, line[6:])
        if outline.split()[-1] == '21' or outline.split()[-1] == '22':
            r1 = outline.split()[1].strip()
            rsp = r1[:-4]
            if rsp in spec:  # Check if r1 exists in the UMIST list.
                umist_BE = BE[spec.index(rsp)]  
                print "%13s:  Old BE = %s -> New BE = %s" % (rsp,outline.split()[-4],umist_BE)
                newline = line[0:86]+"%8.2e" % (N.float(umist_BE))+line[94:]
                fout.write(newline)
                #ctt += 1

            else:
                print "UMIST doesn't have %s." % r1
                fout.write(outline)
        else:
            fout.write(outline)        
        nreac += 1

    fout.close()
    fin.close()
    print 'Successfully matched %i reactions.' % ctt
    #cmd = 'mv -f '+name_reac+'.dat '+name_reac+'_oldpd.dat'
    #os.system(cmd)
    
    #cmd = 'mv -f '+name_reac+'_NewPD.dat '+name_reac+'.dat'
    #os.system(cmd)

nlist = ['rreacs_herb0308_CDR_CDtl_HDO_CO2_v2','rreacs_herb0308_CDR_CDtl_HDO_CO2_ice'] #,'_CR_Deut_Carbon_Dummy']
nlist = ['rreacs_herb0308'] #,'_CR_Deut_Carbon_Dummy']

for j in range(len(nlist)):
    main(nlist[j]) 
