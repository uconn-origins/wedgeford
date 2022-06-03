#!/usr/bin/env python

from math import *
import os
import pdb as pdb
import re
import numpy as N
import sys 
import shutil

#if i == len(m) or m[i+1] is lower:
#this_species = m[i:i+1]
#else
#this_species = m[i]

#search this_species in Mnm

def calcreacmass(reacs):
    Mnm = ['H','O','C','N','He','Si','Mg','Fe','S','D','Na','P','Cl','F','Z','Y']
    Mnu = N.zeros(len(Mnm))
    chg = 0
    for m in reacs:
        i = 0
        il = -1
        if m == 'PHOTON' or m=='GRAIN0' or m=='GRAIN' or m=='CRP' or m=='XRAY' or m=='TEMP' or m=='C-RAY' or m=='LYAPHOTON' or m=='RADNUC':
            continue
        if m[-4:] == '(gr)':
            m = m[:-4]
        if m == 'GRAIN-' or m == 'E':
            chg-=1
            i+=1
            continue            
        while i < len(m):
            if m[i] == '+':
                chg+=1
                i+=1
                continue
            if m[i] == '-': 
                chg-=1
                i+=1
                continue
            if m[i] == 'H':
                try:
                    if m[i+1] == 'e':
                        Mnu[4]+=1
                        i+=1
                        il = 4
                    else: 
                        Mnu[0]+=1
                        il = 0
                except:        
                    Mnu[0]+=1
                    il = 0
            elif m[i] == 'O':        
                Mnu[1]+=1
                il = 1
            elif m[i] == 'P':        
                Mnu[11]+=1
                il = 11
            elif m[i] == 'C':        
                try:
                    if m[i+1] == 'l':
                        Mnu[12] += 1
                        i+=1
                        il = 12
                    else: 
                        Mnu[2]+=1
                        il = 2
                except:        
                    Mnu[2]+=1
                    il = 2               
            elif m[i] == 'N':        
                try:
                    if m[i+1] == 'a':
#                        print 'Na.'
                        Mnu[10]+=1
                        i+=1
                        il = 10
                    else: 
#                        print 'Nitrogen.'
                        Mnu[3]+=1
                        il = 3
                except:        
#                    print 'Nitrogen.' 
                    Mnu[3]+=1
                    il = 3
            elif m[i] == 'D':        
#                print 'Deut.'
                Mnu[9]+=1
                il = 9
            elif m[i] == 'S':
                try:
                    if m[i+1] == 'i':
#                        print 'Silicon.'
                        Mnu[5]+=1
                        il = 5
                        i+=1
                    else: 
#                        print 'Sulfur.'
                        Mnu[8]+=1
                        il = 8
                except:        
#                    print 'Sulfur.' 
                    Mnu[8]+=1
                    il = 8
            elif m[i] == 'M':
                try:
                    if m[i+1] == 'g':
#                        print 'Mg.'
                        Mnu[6]+=1
                        il = 6
                        i+=1
                except:        
                    print "Can't find Mg."
                    sys.exit()   
            elif m[i] == 'F':
                try:
                    if m[i+1] == 'e':
                        Mnu[7]+=1
                        il = 7
                        i+=1
                    else: 
                        Mnu[13]+=1
                        il = 13
                except:        
                    Mnu[13]+=1
                    il = 13
            elif m[i] in '23456789':
                Mnu[il]+=float(m[i])-1
                #pdb.set_trace()
            elif m[i] == '1':
                Mnu[il]+=float(m[i:i+2])-1
                i+=1
                #pdb.set_trace()
            else:
                print "Can't find ",m[i],' of ',m
                sys.exit()            
            i+=1             
        
    return Mnu,Mnm,chg

def calcprodmass(prods,Mnm):
    Mnu = N.zeros(len(Mnm))
    return Mnu

######################################################################

def main(name_reac):

    # Read in entire new rate file with : delimited.
    spec    = []
    BE      = []
    
    ctt = 0

    fin = open(name_reac+'.dat', 'r')
    #fout = open(name_reac+'_NewBE.dat', 'w')
   # Mnm = ['H','O','C','N','He','Si','Mg','Fe','S','D','Na']

    silent = True
    for line in fin:
        if line[0] == "#":
            continue
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:84].strip()
        v1 = line[84:94].strip()
        v2 = line[94:104].strip()
        v3 = line[104:114].strip() 
        type = line[114:118].strip() 
        

        totreac,Mnm,chgr = calcreacmass([r1,r2])
        
        
        totprod,Mnm,chgp = calcreacmass([p1,p2,p3,p4])
        #if line[0:6].strip() == '5123':
        #    pdb.set_trace()
        if N.sum(totreac-totprod) != 0:
            print 'Atoms not conserved: '
            print 'Reac: ',[r1,r2],totreac, ' Chrg: ',chgr
            print 'Prod: ',[p1,p2,p3,p4],totprod, ' Chrg: ',chgp
        if chgr != chgp:    
            print 'Charge not conserved: '
            print 'Reac: ',[r1,r2],totreac, ' Chrg: ',chgr
            print 'Prod: ',[p1,p2,p3,p4],totprod, ' Chrg: ',chgp
        elif not silent:    
            print 'Reac: ',[r1,r2],totreac, ' Chrg: ',chgr
            print 'Prod: ',[p1,p2,p3,p4],totprod, ' Chrg: ',chgp
    fin.close()
    return

def merge_new(rtype_copy,rtype_bran,branchfile,copyfile,mergename):

    f_branch = open(branchfile+'.dat','r')
    br_dat = f_branch.readlines()
    f_copy = open(copyfile+'.dat','r')
    cp_dat = f_copy.readlines()
    f_branch.close()
    f_copy.close()
    
    morebran = [x+2 for x in range(14)]
    brantot = rtype_bran + morebran
    
    nfile = []
    
    for j in range(len(cp_dat)):
        if cp_dat[j][0] == '#':
            continue
        line = br_dat[j]
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:85].strip()
        v1 = line[85:95].strip()
        v2 = line[95:105].strip()
        v3 = line[105:115].strip() 
        type = int(line[115:119].strip() )
        #pdb.set_trace()
        if type in rtype_copy:
            nfile.append(cp_dat[j])
        elif type in brantot:
            nfile.append(br_dat[j])
        elif type == -23 or type == -24:
            if N.float(v1) == 3.14E-10:
                nfile.append(cp_dat[j])        
            else:
            	#cp_dat[j][94:105].strip()
            	newline = cp_dat[j]
            	nfile.append(newline[0:96]+"%9.2E" % (N.float(v1)/3.14E-10 * N.float(v2))+newline[104:])
        elif type == 23:  # grain- reactions need to be branched
            if N.float(v1) == 3.14E-10:
                nfile.append(cp_dat[j])        
            else:
                #pdb.set_trace()
            	#cp_dat[j][94:105].strip()
            	newline = cp_dat[j]
            	nfile.append(newline[0:105]+"%9.2E" % (N.float(v1)/3.14E-10 * 0.5)+newline[114:]) 
        else:
        	print 'Not identified'
        	print r1,r2,p1,p2,type
    fo = open(mergename+'.dat','w')
    for line in nfile:
        fo.write(line)
    fo.close()
    #pdb.set_trace()
    return

######################################################################


def checkisomers(branchfile):
    # Read in entire new rate file with : delimited.
    spec    = []
    BE      = []
    
    ctt = 0

    fin = open(name_reac+'.dat', 'r')
    #fout = open(name_reac+'_NewBE.dat', 'w')
   # Mnm = ['H','O','C','N','He','Si','Mg','Fe','S','D','Na']

    silent = True
    for line in fin:
        if line[0] == "#":
            continue
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:85].strip()
        v1 = line[85:95].strip()
        v2 = line[95:105].strip()
        v3 = line[105:115].strip() 
        type = line[115:119].strip() 
        

        totreac,Mnm,chgr = calcreacmass([r1,r2])
        
        
        totprod,Mnm,chgp = calcreacmass([p1,p2,p3,p4])    
    return

######################################################################
# Extract all molecules from rreacs
def rspecies(rreacs):
    fin = open(rreacs+'.dat','r')

    species = []

    for line in fin:
        if line[0] == '#': continue

        mols = []
        words = line[:86].split()

        for spec in words:
            if spec in species: continue
            try:
                float(spec)
                continue
            except ValueError:
                species.append(spec)

    fin.close()
    return species
######################################################################
# Identify duplicate species 
def comp_speclist(r1,r2):
    s1 = rspecies(r1)
    s2 = rspecies(r2)
    species_unique = []
    
    bothlist = [s1+s2]
    for mol in bothlist[0]:
        if mol in species_unique: 
            continue
        else:
            species_unique.append(mol)
    species = species_unique
    notmols = ['CRP', 'PHOTON', 'XRAY', 'TEMP', 'C-RAY', 'LYAPHOTON','GRAIN','GRAIN0','GRAIN-','E']


    atomall = []
    chgall = []
    
    for mol in species:
        atoms,Mnm,chgr = calcreacmass([mol])
        atomall.append(atoms)
        chgall.append(chgr)
    for j in range(len(species)):
        for k in range(len(species)):
            if j != k:
                if '(gr' not in species[j] and '(gr' not in species[k]:
                    if species[k] not in notmols and species[j] not in notmols:
                        if not N.any(atomall[j] - atomall[k] != 0): 
                            if chgall[j] == chgall[k]:
                                print 'duplicate: ',species[j],species[k]
        #print mol,atoms
        
    # create a unique subset of molecules that are atomically distinct.  identify isomers in the network.  replace all isomers with a single version except for HNC
    
    return

######################################################################
# Compare my file to see if any rates have been updated.
# Loop through major list, see if reactions match, if rates are same, move on, if rates
# are different ask user.
def comp_rates(originalfile,merged):

    # Read in entire new rate file with : delimited.
    spec    = []
    BE      = []
    
    ctt = 0

    fin = open(merged+'.dat', 'r')
    fin_o = open(originalfile+'.dat', 'r')
    myfile = fin_o.readlines()
    fin_o.close()

    #fout = open(name_reac+'_NewBE.dat', 'w')
   # Mnm = ['H','O','C','N','He','Si','Mg','Fe','S','D','Na']
   
    r1_o = []
    r2_o = []
    p1_o = []
    p2_o = []
    p3_o = []
    p4_o = []
    v1_o = []
    v2_o = []
    v3_o = []   
    type_o = []   

    for j in range(len(myfile)):
        line_orig = myfile[j]
        if line_orig[0] == "#":
            continue
        r1_o.append(line_orig[6:19].strip())
        r2_o.append(line_orig[19:32].strip())
        p1_o.append(line_orig[32:45].strip() )
        p2_o.append(line_orig[45:58].strip() )
        p3_o.append(line_orig[58:71].strip() )
        p4_o.append(line_orig[71:85].strip())
        v1_o.append(line_orig[85:95].strip())
        v2_o.append(line_orig[95:105].strip())
        v3_o.append(line_orig[105:115].strip() )
        type_o.append(line_orig[115:119].strip() )
    
    newfile = []
    for line in fin:
        if line[0] == "#":
            continue
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:85].strip()
        v1 = line[85:95].strip()
        v2 = line[95:105].strip()
        v3 = line[105:115].strip() 
        type = line[115:119].strip() 
        m1 = [1 if x == r1 else 0 for x in r1_o]
        m2 = [1 if x == r2 else 0 for x in r2_o]
        m3 = [1 if x == p1 else 0 for x in p1_o]
        m4 = [1 if x == p2 else 0 for x in p2_o]
        m5 = [1 if x == p3 else 0 for x in p3_o]
        m6 = [1 if x == p4 else 0 for x in p4_o]
        
        indfound = N.where((N.array(m1)==1)&(N.array(m2)==1)&(N.array(m3)==1)&(N.array(m4)==1)&(N.array(m5)==1)&(N.array(m6)==1))[0]
        
        numfind = len(indfound)
        
        if numfind == 1:
            if type == type_o[indfound[0]]:
                if N.float(v1) != N.float(v1_o[indfound[0]]) or N.float(v2) != N.float(v2_o[indfound[0]]) or N.float(v3) != N.float(v3_o[indfound[0]]):
                    print r1,r2,p1, r1_o[indfound[0]], r2_o[indfound[0]], p1_o[indfound[0]],':',v1,v2,v3,'<-', v1_o[indfound[0]], v2_o[indfound[0]], v3_o[indfound[0]]
                    if r1 == 'Si+':
                        pdb.set_trace()
            else:
                newfile.append(line)
        elif numfind == 0:
            #print "Nothing found in original",line
            newfile.append(line)
        else:
            
            if type == '13' or type == '15' or type == '30'  or type == '-1':
                newfile.append(line)
            else:
                print 'too many finds'
                pdb.set_trace()
        
       #pdb.set_trace()
        
        

    fin.close()
    return
    
 
 ######################################################################
def add_power_react(originalfile,merged):

    # Read in entire new rate file with : delimited.
    spec    = []
    BE      = []
    
    ctt = 0

    fin = open(originalfile+'.dat', 'r')
    fin_o = open(merged+'.dat', 'r')
    myfile = fin_o.readlines()
    fin_o.close()

    #fout = open(name_reac+'_NewBE.dat', 'w')
   # Mnm = ['H','O','C','N','He','Si','Mg','Fe','S','D','Na']
   
    r1_o = []
    r2_o = []
    p1_o = []
    p2_o = []
    p3_o = []
    p4_o = []
    v1_o = []
    v2_o = []
    v3_o = []   
    type_o = []   

    for j in range(len(myfile)):
        line_orig = myfile[j]
        if line_orig[0] == "#":
            continue
        r1_o.append(line_orig[6:19].strip())
        r2_o.append(line_orig[19:32].strip())
        p1_o.append(line_orig[32:45].strip() )
        p2_o.append(line_orig[45:58].strip() )
        p3_o.append(line_orig[58:71].strip() )
        p4_o.append(line_orig[71:85].strip())
        v1_o.append(line_orig[85:95].strip())
        v2_o.append(line_orig[95:105].strip())
        v3_o.append(line_orig[105:115].strip() )
        type_o.append(line_orig[115:].strip() )
    
    newfile = []
    for line in fin:
        if line[0] == "#":
            continue
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:85].strip()
        v1 = line[85:95].strip()
        v2 = line[95:105].strip()
        v3 = line[105:115].strip() 
        type = line[115:].strip() 
        m1 = [1 if x == r1 else 0 for x in r1_o]
        m2 = [1 if x == r2 else 0 for x in r2_o]
        m3 = [1 if x == p1 else 0 for x in p1_o]
        m4 = [1 if x == p2 else 0 for x in p2_o]
        m5 = [1 if x == p3 else 0 for x in p3_o]
        m6 = [1 if x == p4 else 0 for x in p4_o]
        
        indfound = N.where((N.array(m1)==1)&(N.array(m2)==1)&(N.array(m3)==1)&(N.array(m4)==1)&(N.array(m5)==1)&(N.array(m6)==1))[0]
        
        numfind = len(indfound)
        
        if numfind == 0:
            print line[:6],r1,r2,p1,p2,p3,p4
            newfile.append(line)
    fin.close()

    fo = open(merged+'_bonus.dat','w')
    for line in myfile:
        fo.write(line)    
    for line in newfile:
        fo.write(line)
    fo.close()

    return
######################################################################
# Calculate molecules masses for reactions that need that info.        
def resetmasses(rtype,merged):
    fin = open(merged+'_bonus.dat','r')
    dat = fin.readlines()
    fin.close()
    molmass = [1.008,15.999,12.011,12.007,4.0026,28.086,24.305,55.845,32.065,2.0,22.990,20.974,35.453,18.998,0.0,0.0]
    
    newfile = []


    for line in dat:
        if line[0] == "#":
            newfile.append(line)
            continue
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:85].strip()
        v1 = line[85:95].strip()
        v2 = line[95:105].strip()
        v3 = line[105:115].strip() 
        type = line[116:].strip() 
        if N.int(type) not in rtype:
            newfile.append(line)
            continue
        totreac,Mnm,chgr = calcreacmass([r1])
        mas = N.sum(totreac*molmass)
        update = line[:95]+"%10.2E" % (mas) + line[105:]
        print update[0:30]+type
        newfile.append(update)

    fo = open(merged+'_masscal.dat','w')
    for line in newfile:
        fo.write(line)
    fo.close()
    
    return


######################################################################
# Calculate molecules masses for reactions that need that info.        
def checkdupes(file):
    fin = open(file,'r')
    dat = fin.readlines()
    fin.close()    
#    newfile = []

    #c = 0
    
    
    for c in range(len(dat)):
        line = dat[c]
        if line[0] == "#":
#            newfile.append(line)
            continue
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:85].strip()
        v1 = line[85:95].strip()
        v2 = line[95:105].strip()
        v3 = line[105:115].strip() 
        type = line[116:].strip() 
            
        for b in range(len(dat)):
            line2 = dat[b]
            r12 = line2[6:19].strip()
            r22 = line2[19:32].strip()
            p12 = line2[32:45].strip() 
            p22 = line2[45:58].strip() 
            p32 = line2[58:71].strip() 
            p42 = line2[71:85].strip()
            v12 = line2[85:95].strip()
            v22 = line2[95:105].strip()
            v32 = line2[105:115].strip() 
            type2 = line2[116:].strip()    
            
            if c != b and r12 == r1 and r22 == r2 and p12 == p1 and p22 == p2 and p32 == p3 and p42 == p4:
                print r1,r2,p1,p2,p3,p4,b,c
                     
        

    return

######################################################################
#  Calculate 6grain.be file assuming same BE
def calcbe(merged):

    s1 = rspecies(merged+'_masscal')
    species_unique = []

    fin = open('6grainbe.inp','r')
    dat = fin.readlines()
    fin.close()
    
    names = [(x.strip()).split()[0] for x in dat]
    vals = [(x.strip()).split()[1] for x in dat]  
    
    storeprint = []
    storemol = []
    storeorigBE = []
    for line in dat:
        if line[0] != '#':
            totreac,Mnm,chgr = calcreacmass([line.split()[0]])
            storeprint.append(totreac)  
            storemol.append(line.split()[0])  
            storeorigBE.append(line.split()[1])   
    
    storeBE = []

    for line in s1:
        if '(gr' in line:
            fnd = len(N.where(N.array(names) == line)[0])
            if fnd == 0:
                new = line.replace('D','H')
                if new in names:
                    BEhere = vals[names.index(new)]
                    storeBE.append("%s%13s" % (line.ljust(13),BEhere.ljust(13)))
                else:
                    totreac,Mnm,chgr = calcreacmass([new])
                    print new
                    for i in range(len(storemol)):
                        if not N.any(totreac - storeprint[i] != 0): 
                            storeBE.append("%s%13s" % (line.ljust(13),BEhere.ljust(13)))
                            print 'found',storemol[i]
            if fnd != 0:
                BEhere = vals[names.index(line)]      
                storeBE.append("%s%13s" % (line.ljust(13),BEhere.ljust(13)))   
                # replace all instances with H, then find same atom pattern
    
    
    fo = open('6grainbe.inp.fulldeut','w')
    for line in storeBE:
        fo.write(line+'\n')
    fo.close()

    return
######################################################################
def correctphoto(merged):
    fin = open(merged+'_masscal.dat','r')
    dat = fin.readlines()
    fin.close()    

    names = []
    for file in os.listdir("xsect"):
        if file.endswith(".photoxs"):
            tm = file.index('.p')
            names.append(file[:tm])
            
    storeprint = []
    for mol in names:
        totreac,Mnm,chgr = calcreacmass([mol])
        storeprint.append(totreac)   
        
        
    #names = [(x.strip()).split()[0] for x in dat]
    #vals = [(x.strip()).split()[1] for x in dat] 
    
    for line in dat:
        if line[0] == "#":
            continue
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:85].strip()
        v1 = line[85:95].strip()
        v2 = line[95:105].strip()
        v3 = line[105:115].strip() 
       # type = line[116:].strip() 
        type  = line.strip().split()[-1]

        if N.int(type) == 30:
            if not os.path.isfile('xsect/%s.photoxs'%r1):
                new = r1.replace('D','H')
                if new in names:
                    print 'found',new
                    shutil.copy2('xsect/%s.photoxs'%new, 'xsect/%s.photoxs'%r1)
                else:
                    totreac,Mnm,chgr = calcreacmass([new])
                    print 'No Direct Match, checkign masses',new
                    for i in range(len(names)):
                        if not N.any(totreac - storeprint[i] != 0): 
                            shutil.copy2('xsect/%s.photoxs'%(names[i]), 'xsect/%s.photoxs'%r1)
    
    return
######################################################################
def reformatchem(merged):
    fin = open(merged+'_masscal.dat','r')
    dat = fin.readlines()
    fin.close() 

    fo = open(merged+'_reformat.dat','w')
    #for line in storeBE:
    
    for line in dat:
        if line[0] == '#':
            fo.write(line)
            continue
        r0 = line[1:6]
        r1 = line[6:19].strip()
        r2 = line[19:32].strip()
        p1 = line[32:45].strip() 
        p2 = line[45:58].strip() 
        p3 = line[58:71].strip() 
        p4 = line[71:85].strip()
        v1 = line[85:95].strip()
        v2 = line[95:105].strip()
        v3 = line[105:115].strip() 
        type = line[115:].strip()         #pdb.set_trace()
        #2X,E8.2,1X,E9.2,1X,E9.2,1X,I3
        new = "%5s %13s%13s%13s%13s%13s%13s%10s%10s%10s %3s\n" % (r0.rjust(5),r1.ljust(13),r2.ljust(13),p1.ljust(13),p2.ljust(13),p3.ljust(13),p4.ljust(13),v1,v2,v3,type)
        fo.write(new)
    fo.close()
      
    return

nlist = ['rreacs_herb0308_LD_GG_thinoC_bg']
#nlist = ['rreacs_herb0308_LgDeut_GG_BE']
nlist = ['rreacs_herb0308_planet']
nlist = ['rreacs_mergenew_reformat_CH2DP_grn']

original = 'rreacs_herb0308_LgDeut_GG_BE_thi_v2'
branch = 'rreacs_normal_isotopized_branch_1118v1'
copy = 'rreacs_normal_isotopized_copy_1118v1'
mergeout = 'rreacs_mergenew'

rtype_copy = [20,21,22,30,42,43,44,45,47,49,-1,-2,-25,0]
rtype_bran = [41,46,48,50,51,99] # and 2-15 inclusive.
rtype_mass = [20,21,23,47]

#20,21,23,47 -- beta column is mass
#-23,-24 is yield
# copy input files for species with 30 (only needs to happen once)
# 

# Combine the two new networks with Copy vs. Branch based on IDs above
# Check phase changes do not change atomic structure
# Compare new and old network for isomers and check whether reaction is branched

for j in range(len(nlist)):
    print 'Analyzing: ',nlist[j]+'.dat'
    #merge_new(rtype_copy,rtype_bran,branch,copy,mergeout)
    #comp_speclist(original,mergeout) # check for isomers
    #comp_rates(original,mergeout) # check old rates vs. new rates
    #add_power_react(original,mergeout)
    #adjust_mass(mergeout)
    #resetmasses(rtype_mass,mergeout)
    #calcbe(mergeout)  # create a 6BE file
    #correctphoto(mergeout) # create missing xsect files for reactions with 30-type
    #reformatchem(mergeout)
    checkdupes(nlist[j]+'.dat')
    
    main(nlist[j]) 
