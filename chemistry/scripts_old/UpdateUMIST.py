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
    inp = '/Shire1/cleeves/Research/PhotoDissUMIST2012/RATE12.dist.txt'
    fumist = open(inp, 'r')
    dat = fumist.readlines()
    fumist.close()
    
    line_ph = []
    reac    = []
    p1      = []
    p2      = []
    p3      = []
    p4      = []
    rate    = []
    
    ctt = 0
    
    # Read in lines from file data, dat:
    for line in dat:
        spln = re.split(':',line)
        if spln[1] == 'PH' and line != '\n' and line != '':
            line_ph.append(line.strip())
            reac.append(spln[2])
            
            if 'e-' in spln[4:8]:
                spln[spln.index('e-')] = 'E'
            
            p1.append(spln[4])
            p2.append(spln[5])
            p3.append(spln[6])
            p4.append(spln[7])
            rate.append(N.float(spln[9]))
    dir = '/Shire1/cleeves/Research/PhotoDissUMIST2012/'
    fin = open(name_reac+'.dat', 'r')
    fout = open(name_reac+'_NewPD.dat', 'w')

    nreac = 1
    for line in fin:
        if line[0] == "#":
            fout.write(line)
            continue
        outline = " %4d %s" % (nreac, line[6:])
        if outline.split()[-1] == '13':
            r1 = outline.split()[1]
            if r1 in reac:  # Check if r1 exists in the UMIST list.
                a_i = all_indices(r1, reac)  # Return all instances of r1 photodissociation.
                in_prod = [outline[32:45].strip(),outline[45:58].strip(),outline[58:71].strip(),outline[71:84].strip()]
                flag = 0
                a_found = -1
                
                for p in range(len(a_i)):
                    umistp = [p1[a_i[p]],p2[a_i[p]],p3[a_i[p]],p4[a_i[p]]]
                    if (p1[a_i[p]] in in_prod) and (p2[a_i[p]] in in_prod) and (p3[a_i[p]] in in_prod) and (p4[a_i[p]]) in in_prod:
#                        print 'Matched reaction! Found: ',umistp
                        a_found = a_i[p]
                        flag = 1
                        break
                if flag == 0:
                    print 'Was unsuccessful at matching.  ',r1,'goes to: ',in_prod
                    fout.write('\\_'+outline)
                else:    
#                    print "%13s:  Old rate = %s -> New rate = %13.2e" % (r1,outline.split()[-4],rate[a_found])
                    newline = line[0:86]+"%8.2e" % rate[a_found]+line[94:]
                    fout.write(newline)
                    ctt += 1
            else:
                print "UMIST doesn't have %s." % r1
                fout.write(outline)
        else:
            fout.write(outline)        
        nreac += 1

    fout.close()
    fin.close()
    print 'Successfully matched %i reactions.' % ctt
    cmd = 'mv -f '+name_reac+'.dat '+name_reac+'_oldpd.dat'
    os.system(cmd)
    
    cmd = 'mv -f '+name_reac+'_NewPD.dat '+name_reac+'.dat'
    os.system(cmd)

nlist = ['rreacs_herb0308_CDR_CDtl_HDO_CO2_v2','rreacs_herb0308_CDR_CDtl_HDO_CO2_ice'] #,'_CR_Deut_Carbon_Dummy']
for j in range(len(nlist)):
    main(nlist[j]) 
