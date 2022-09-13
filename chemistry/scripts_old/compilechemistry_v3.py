# Compile the abundances into a LIME and plotting friendly single file format
# from a disk_chemistry "runs" directory.
# Dust parameter (here = 'e1') tells the code to look under the correct folder:
# disk_chemistry/runs/MODELNAME/e1/... and does nothing else.
#
# Speed optimized to the best of my ability, but depends on how large your file
# is and the complexity of the molecule you are extracting (how far down it is
# in the file).
# if the output files had a header it would look like this:
# R(m) Z(m) n_H2(number of H2 per m^3) Tgas(K) Tdust(K) abundance-relative-to-H2


import os
import fnmatch
import matplotlib
import matplotlib.pyplot as plt
import optparse, time
import numpy as N
import pdb as pdb
import sys
from scipy import *
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
from matplotlib import rc
from matplotlib import rcParams
from pylab import *
from subprocess import call


def radreturn(bsd,dust,mod,strmol,steps):

    loc1 = bsd+mod+dust
    cmd = 'rm -f '+loc1+'/listI.txt; cd '+loc1+ '; ls -d r*_1.out | awk \'{g=$0; sub(/r/,"",g); split(g, a, "_"); printf("%s\\n",a[1]);}\' > listI.txt'
    call(cmd,shell=True)
    print(loc1)
    f = open(loc1+'/listI.txt')
    radstr = f.readlines()
    f.close()

    nR = len(radstr)
    radval = N.zeros(nR)
    for i in range(nR):
        radstr[i] = radstr[i].rstrip()
        radval[i] = N.float(radstr[i])
        #arrtmp = genfromtxt(baseenvdir+'/1environ.inp.e1.'+radstr[i],skip_header=3)

    sorti = N.argsort(radval)
    radval = radval[sorti]
    radnam = []
    for i in range(nR):
        radnam.append(radstr[sorti[i]])


    # Radnam = names for files.  Radval == radius values in AU taken from file names.
    nZ = 0
    # Find maximum z value from file:
    filetemp = 'r'+radnam[0] +'_'+dust+'_*.out'
    for file in os.listdir(loc1):
        if fnmatch.fnmatch(file, filetemp):
            vals = file[-6:len(file)-4]
            if vals[0] == '_':
                vals = vals[1]
            vals = N.int(vals)
            if vals > nZ:
                nZ = vals
    file = 'r'+radnam[0] +'_'+dust+'_'+str(nZ)+'.out'

    listfound = []
    with open(loc1 +'/'+file) as myFile:
        for num, line in enumerate(myFile, 1):
            if strmol in line:
                listfound.append(num)

    if len(listfound) == 3:
        print('Found molecule in file.')
    else:
        print( 'Error finding molecule: '+strmol+'. Returning.')
        return
    myFile.close()
    molindex = listfound[1]
    nTimes = steps
    Times = N.zeros(nTimes)
    # ------------------------------------------
    # Find Column:
    # First read everything before this point:
    fin = open(loc1 +'/'+file, 'r')
    for i in range(molindex-1):
        fin.readline()
    headerline = fin.readline()
    ist = headerline.split().index(strmol.strip())
    print('Extracting: '+headerline.split()[ist]+'!')
    fin.readline()
    for i in range(nTimes):
        tmp = fin.readline()
        tmp = N.float(tmp[0:12])
        Times[i] = tmp
    fin.close()
    # And now we've found the location.

    # ------------------------------------------
    # Find Physical Information:
    # First read everything before this point:
    listfound = []
    with open(loc1 +'/'+file) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'INITIAL VALUES' in line:
                infind = num
    myFile.close()

    heightind = infind+4
    rhoind = infind+8

    # ------------------------------------------

    cti = 0
    StoreWrite = N.zeros((nR*nZ,7,nTimes))
    StoreAbundance = N.zeros((nR,nZ,nTimes))
    StorePhysical = N.zeros((nR,nZ,3))
    # For each value of z, extract abundance.
    print( "Data dimensions (nr, nz): ("+str(nR)+","+str(nZ)+")")
    start_time = time.time()
    for j in range(nR):
        for i in range(nZ):
            file = 'r'+radnam[j] +'_'+dust+'_'+str(i+1)+'.out'

            fin = open(loc1 +'/'+file, 'r')
            # Get physical info from environment file:
            for g in range(heightind-2):
                fin.readline()

            rauline = fin.readline()
            heightline = fin.readline()
            for g in range(rhoind-heightind-3):
                fin.readline()
            Tgline = fin.readline()
            Tdline = fin.readline()
            rholine = fin.readline()

            rau = N.float(rauline[20:29])
            zau = N.float(heightline[20:29])
            rho = N.float(rholine[20:29])
            Tg = N.float(Tgline[20:26])
            Td = N.float(Tdline[20:26])

            # Density set to zero at top-most zone.
            if i == 0:
                rho = 0.0

            StorePhysical[j,i,0] = rau
            StorePhysical[j,i,1] = zau
            StorePhysical[j,i,2] = rho

            StoreWrite[cti,0,:] = rau
            StoreWrite[cti,1,:] = zau
            StoreWrite[cti,2,:] = rho
            StoreWrite[cti,3,:] = Tg
            StoreWrite[cti,4,:] = Td

            # Rewind the file to the beginning:
            fin.seek(0)

            # Skip intermediate information and go directly to molecular data line:
            for g in range(molindex):
                fin.readline()
            fin.readline() # Garbage --- line
            tv = N.zeros(nTimes)
            for g in range(nTimes):
                tmp = fin.readline()
                tmp = tmp.split()[ist]
                tmp = N.float(tmp.rstrip())
                if isnan(tmp):
                    tmp = 0.0
                tv[g] = tmp
            StoreAbundance[j,i,:] = tv
            StoreWrite[cti,5,:] = tv
            StoreWrite[cti,6,:] = tv
            cti = cti + 1
            fin.close()


    print( "Time elapsed: ", time.time() - start_time, "s")
    return nR,nZ,radnam,radval,loc1,StoreAbundance,StorePhysical,Times,StoreWrite


def writemol(bsd,dust,mod,mol,steps):

    nmol = len(mol)
    # Constant Block:
    mp = 1.67e-24
    cmAU = 1.496e13
    mu = 2.36
    pi = 3.141596
    mau = 1.496e11

    # Converts data to lime format, meters for R, Z, NH2/m^3 for collisional partners, abundance from X_H to X_H2.



    dataformat = '%15.7E %15.7E %15.7E %15.7E %15.7E %15.7E'

    for t in range(nmol):
        molname = ' '+mol[t]+' ' # Pad sides for searching.
        dirtest = bsd+mod+dust+'/lime_'+molname[1:-1]
        testpath = os.path.exists(dirtest)
        if not testpath:
            os.mkdir(dirtest)

        nR,nZ,radnam,radval,loc1,SAb,SPhy,Time,SaveArr = radreturn(bsd,dust,mod,molname,steps)
        print( SaveArr.shape)
        print( Time.shape)

        for j in range(N.size(Time)):
            fout=open(dirtest+'/'+molname[1:-1]+'_time'+str(j)+'.dat','w')
            for i in range(nR*nZ):
                denhere = SaveArr[i,2,j]/(mu*mp)*1e6
                if denhere < 1e10: # simon's density cut off was 10^5 cm^-3, but here 10^4 conservatively.
                    denhere = 0.0
                formatted = dataformat % (SaveArr[i,0,j]*mau,SaveArr[i,1,j]*mau,denhere,SaveArr[i,3,j],SaveArr[i,4,j],SaveArr[i,5,j]*2.0)
                fout.write(formatted+'\n')
            fout.close()

def main():

    # (1) Base directory for disk_chemistry/runs/:
    bsd = '/data/beegfs/astro-storage/groups/henning/schwarz/wedgeford/chemistry/runs/'
    # Subfolder labeled by dust assumed.  ProtoDisk (or GenDisk)
    # already incorporate settling, but always labeled as 'e1'.
    dust = 'e1'

    # (2) Define a list (or single molecule ['CO']) to extract:
    mols =['CO','CO(gr)','CH3OH(gr)','CO2(gr)','H2O(gr)']

#     file = os.listdir(bsd)
#     f = [i for i,x in enumerate(file) if x ==' 6Myr_new_0.03Mslg_0.8m_d200.0_Tgas_ssx_1g0_x31_noH2O']
#
    #file = file[f[0]+1:]
    mod = ['newdustchem_0.03Mslg_0.5m_ssx/']

    # file1 = '/home/kschwarz/MasterChemistry/rspecies_herb0308_CDR_FullGr.dat'
    # infile = open(file1,"r")
    # mols=[]
    # for line in infile:
    # 	if line.count('C') > 0 and line.count('#') == 0 and line.count('gr') == 1 and line.count('Cl')  != line.count('C'):
    # 		mols.append(line.strip('\n'))#need to strip line ending
    # infile.close()

    #file1 = '/nfs/ebergin1/kamberrs/disk_chemistry/MasterChemistry/rspecies_herb0308_CDR_FullGr.dat'
    #infile = open(file1,"r")
    #for line in infile:
    # 	if line.count('O') > 0 and line.count('PHOTON') == 0 and line.count('gr') == 0:
    # 		mols.append(line.strip('\n'))#need to strip line ending
    #infile.close()


    # Next define the name of the run created with runchemistry.py.
    # This will be the name of the folder under disk_chemistry/runs/
    # which you wish to extract.  Then run "writemol" which will then
    # look at all the successful chemical output in your runs
    # directory, and then create a subfolder "lime_NAME" containing
    # N = 50 files for each time step in the chemistry.  I typically
    # use the file in the 45th time step, e.g. "CO_time45.dat."
    #
    # If looking at more than one chemical run just include them below
    # in sequence.


   # steps = np.ones(len(mod),dtype=int)*180
    #for file in os.listdir(bsd):
    # 	if file.count('drl300') == 1:
    # 		writemol(bsd,dust,file+'/',mols,180)

    for n in range(len(mod)):
	       writemol(bsd,dust,mod[n],mols,180)

    return

main()
