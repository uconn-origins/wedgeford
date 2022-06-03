#!/usr/bin/env python

import pdb as pdb
import os, time
import matplotlib
import matplotlib.pyplot as plt
import optparse,sys
import numpy as N
from scipy import *
from pylab import *
from matplotlib import rc
import fnmatch

############################################################
# Code looks for all 20 wavelengths and replaces erroneous midplane values with value of X-ray
# field above midplane for found files.  Returns successful file list, energy levels and type
# of dust used to generate radiation fields.

def fixmidplane(basedir,bsname):
    numwav = 9
    tag = 'Output/dir_'+bsname+'/'
    endst=''

    lamid = N.zeros(numwav)
    filelist = []
    fileexample = ''
    oldline = ''

    for i in range(numwav):
        fname = "%s_e*_E%02i*.txt" % (bsname,i+1)
        print( fname)
        fout = "%sfx_%s_E%02i.txt" % (basedir+tag,bsname,i+1)
        foundfile = ''
        dups = 0
        for file in os.listdir(basedir+tag):
            if fnmatch.fnmatch(file,fname):
                if dups == 1:
                    time2=os.path.getmtime(basedir+tag+file)
                    if time2>time1:
                        foundfile = file
                        time1=time2
                    lamid[i]=lamid[i]+1.0
                else:
                    foundfile = file
                    lamid[i]=lamid[i]+1.0
                    time1 = os.path.getmtime(basedir+tag+foundfile)
                    dups = 1
                fileexample = file
        if foundfile != '':
            f=0
            filelist.append(fout)
            f = open(basedir+tag+foundfile, 'r')
            fo = open(fout,'w')
            ct = 0
            for line in f:
                if ct <= 5:
                    cline = line # f.readline()
                    fo.write(cline)
                else:
                    cline = line
                    splint = cline.split()
                    zv = N.float(splint[1])
                    xv = N.float(splint[0])
                    pv = N.float(splint[2])
                    if zv/xv <= 1e-6:   # If point is at midplane!
                        if oldline == '':
                            print( "Height values are in increasing, not decreasing order. Please reformat model and try again.\n")
                            sys.exit()
                        splold = oldline.split()
                        pvold = N.float(splold[2])
                        lineout = "%13.4f%13.4f%13.5e\n" % (xv,zv,pvold)
                        if zv > N.float(splold[1]):
                            print( "Height values are in increasing, not decreasing order. Please reformat model and try again.\n")
                            sys.exit()
                        fo.write(lineout)
                    else:
                        fo.write(cline)
                    oldline = line
                ct += 1

            f.close()
            fo.close()

    Eval = N.array(N.where(lamid>0))+1
    Eval = Eval[0]

    # Identify dust setting used in initial run from file name:
    if '_e1_' in fileexample:
        dust = 'e1'
    elif '_e0p1_' in fileexample:
        dust = 'e0p1'
    elif '_e0p01_' in fileexample:
        dust = 'e0p01'
    elif '_e0p001_' in fileexample:
        dust = 'e0p001'
    else:
        print( "Could not identify dust type, exiting.\n")
        sys.exit()

    return filelist, Eval, dust

############################################################

def readUV(filelist):

    c = 2.99792458e10
    dims = genfromtxt(filelist[0],skip_header=6).shape
    photar = N.zeros((dims[0],dims[1],len(filelist)))
    for i in range(len(filelist)):
        photar[:,:,i] = genfromtxt(filelist[i],skip_header=6)
        photar[:,2,i] = photar[:,2,i]*c
        photar[np.where(photar < 0.0)] = 0.0
    return photar

############################################################
def movingaverage(interval, window_size):
    window= N.ones(int(window_size))/float(window_size)
    return N.convolve(interval, window, 'same')
############################################################

def applyspectra(basedir,fluxarr,Eval,L_UV=1.0):

    angs_2_erg  = 1.6021765e-9
    dTW = 55.0*3.08568e18
    h = 6.6260755e-27
    c = 2.99792458e10 # cm/s

    dims = fluxarr.shape
    outputE = N.arange(0, 2010.0-930, 10)+930.0
    lamlist = genfromtxt(basedir+'/wvlnfile_uv.txt')
    #lameval = [N.append(x) for
    lameval = []
    for j in range(len(Eval)):
        iL = N.where(lamlist[:,0] == Eval[j])
        iL = iL[0]
        iL = iL[0]
        lameval.append(lamlist[iL,1])

    # Input file in erg/cm^2/s/A at TW Hya distance.
    specdir = basedir+'/pfiles/'
    template = genfromtxt(specdir+'uv_twh_smooth.dat')
    #template = genfromtxt(specdir+'uv_hd163_smooth.dat')



    template = template[template[:,0] < 2050.0,:]

    en = (h*c/(template[:,0]*1e-8))
    fluxeseval = template[:,1]/en  # phot/cm^2/s/A  - was in erg/cm^2/s/A observed.
    lumsum = 4.0*N.pi*dTW**2*fluxeseval  # phot/s/A <- energy output at each wavelength per A.

    iBelLy = N.argmin(abs(template[:,0]-1192.0))
    iAbvLy = N.argmin(abs(template[:,0]-1232.0))

    lumcont = lumsum*1.0
    lumcont[iBelLy:iAbvLy+1] = (lumsum[iBelLy]+lumsum[iAbvLy])/2.0

    lymline = lumsum[iBelLy:iAbvLy+1]#-(lumsum[iBelLy]+lumsum[iAbvLy])/2.0

    lym_x = template[iBelLy:iAbvLy+1,0]
    Lum_Lya = N.trapz(lymline,x=lym_x) # phot/s in lyman alpha
    Lum_cont = N.trapz(lumcont,x=template[:,0])

    midrng_y = lumcont[template[:,0] < 1600.0]
    midrng_x = template[template[:,0] < 1600.0,0]
    pa = N.polyfit(N.log10(midrng_x),N.log10(midrng_y),1)
    cont_fit = pa[0]*N.log10(template[:,0])+pa[1]
    cont_fit = 10**cont_fit

    movingy = movingaverage(lumcont,  200)

    #fig = plt.figure(1)
    #fig.clf()
    #plt.semilogy(template[:,0],lumcont)
    #plt.semilogy(template[iBelLy:iAbvLy+1,0],lumsum[iBelLy:iAbvLy+1],'r')
    #plt.semilogy(template[:,0],cont_fit)
    #plt.semilogy(template[:,0],movingy)
    #plt.xlim([1150,1300]) #(1200-1230)
    #plt.ylim([1e37,3e39])
    #plt.show()

    #pdb.set_trace()
    expandphot = N.zeros((dims[0],dims[1],len(outputE)))

    for i in range(len(outputE)):
        lam = outputE[i]
        if lam <= lamlist[N.min(Eval)-1,1]:
            expandphot[:,:,i] = fluxarr[:,:,0]
        elif lam >= lamlist[N.max(Eval)-1,1]:
            expandphot[:,:,i] = fluxarr[:,:,-1]
        else:
            #pdb.set_trace()
            Earr = N.concatenate(([lam],lameval))
            Earr=N.sort(Earr)
            iE=N.where(Earr == lam)
            iE = iE[0]
            Ebel = Earr[iE[0]-1]
            Eabv = Earr[iE[0]+1]
            Edist = Eabv-Ebel
            iEbel = N.where(lameval == Ebel)
            iEabv = N.where(lameval == Eabv)
            iEbel = iEbel[0]
            iEabv = iEabv[0]
            A1 = fluxarr[:,:,iEbel]
            A1 = A1[:,:,0]
            A2 = fluxarr[:,:,iEabv]
            A2 = A2[:,:,0]
            expandphot[:,0:1,i] = fluxarr[:,0:1,0]
            expandphot[:,:,i] = A1*N.float(abs(Eabv-lam)/Edist) + A2*N.float(abs(Ebel-lam)/Edist)  # Linear interpolation between two nearest computed wavelengths.
            #pdb.set_trace()
        if lam == 1210.0:
            scx = 10**(N.log10(Lum_Lya)-44.0)/10.0 # phot/s/A (integrated over wavelength previously).
            expandphot[:,0:1,i] = fluxarr[:,0:1,0]
            expandphot[:,2,i] = expandphot[:,2,i]*scx
        else:
            scx = 10**(N.log10(N.interp(lam,template[:,0],cont_fit))-44.0) # phot/s/A @ specified wavelength - already has A dependence from input spectrum.
            expandphot[:,0:1,i] = fluxarr[:,0:1,0]
            expandphot[:,2,i] = expandphot[:,2,i]*scx
        #print lam, scx

    #pdb.set_trace()
    return expandphot,outputE

############################################################
# Eval = Energies in keV RT was evaluted, i.e, 1,3,7,10,15,20
# calfluxarr = Calibrated X-ray spectrum at the above wavelengths.

def format4chem(calfluxarr,outE,basedir,bsname,dust):
    dims = calfluxarr.shape

    tag = 'Output/dir_'+bsname+'/'
    #fo = open(basedir+tag+'Spectra_'+bsname+'.dat','w')

    datarr = N.zeros((dims[0],2+dims[2]))

 # ---------------------------------
 # Extract the radii values.  Precision for Rad values out of
 # Tom's code not good enough - go to inputs
    dir_look4inp = basedir+'Output/dir_'+bsname+'/inputs/'

    fname = "inp_*.txt"
    foundfile = ''
    dups = 0
    for file in os.listdir(dir_look4inp):
        if fnmatch.fnmatch(file,fname):
            foundfile=file
            break
    inputtoToms = genfromtxt(dir_look4inp+foundfile,skip_header=13)

    radvals = inputtoToms[:,1]
    tn = N.size(radvals)
    ct = 1
    radstore = N.zeros(tn)
    radstore[0]=radvals[0]
    for i in range(tn-1):
        checkpresence = N.size(N.where(radstore == radvals[i+1]))
        if checkpresence <= 0:
            radvals[i+1]
            radstore[ct] = radvals[i+1]
            ct = ct + 1


    fz = N.where(radstore == 0)
    fz = fz[0]
    radstore = radstore[0:fz[0]]

    nr = radstore.shape
    nr = nr[0]
    radstore = N.sort(radstore)
    nz = int(tn/nr)
    print( 'Model size, Nr: %i, Nz: %i \n' % (nr,nz))
    nwv = 108

    rawdata = N.zeros((nr*(2+nwv),nz+1))
    checkflag = N.zeros(nwv)
#   Reformat scaled output
    for g in range(nr):
        d = N.where(radvals == radstore[g])
        d = d[0:nz]
        rawdata[g*(2+nwv),:] = N.concatenate(([radstore[g]],N.zeros(nz)))
        dathere = inputtoToms[d,:]
        rawdata[g*(2+nwv)+1,:] = N.concatenate(([0],dathere[0,:,2])) #[0,(transpose(beth(1,d(0)+(nr*nz):d(0)+nz-1+(nr*nz))))]
        for p in range(nwv):
            photdat = calfluxarr[d,2,p]
            rawdata[p+(2+nwv)*g+2,:] = N.concatenate(([outE[p]],photdat[0,:]))

    fo = open(basedir+tag+'AllOut_'+bsname+'.dat','w')
    for g in range(nr):
        fo.write(''.join("%13.5f" % (x) for x in rawdata[g*(2+nwv),:])+"\n")
        fo.write(''.join("%13.5f" % (x) for x in rawdata[g*(2+nwv)+1,:])+"\n")
        for p in range(nwv):
            lineout=rawdata[p+(2+nwv)*g+2,1:-1]
            fo.write('%13.5f' % (rawdata[p+(2+nwv)*g+2,0]) +''.join("%13.5e" % (x) for x in lineout)+"\n")
    fo.close()

    # Now write the pretty version:

    fo = open(basedir+tag+'uv_photons_'+bsname+'_'+dust+'.dat','w')
    fo.write('Results from Bethell code with run: '+bsname+'\n')
    for g in range(nr):
        fo.write('Radius(AU) %9.4f' % radstore[g] + '\n')
        fo.write('z(AU)  '+''.join("%10.4f" % (x) for x in rawdata[g*(2+nwv)+1,1:])+'\n')
        fo.write('Wavelength A               Photons/cm2/s/A\n')
        for p in range(nwv):
            lineout=rawdata[p+(2+nwv)*g+2,1:]
            fo.write('%7.2f' % (rawdata[p+(2+nwv)*g+2,0]) +''.join("%10.2e" % (x) for x in lineout)+"\n")
    fo.close()

    return
############################################################
def compionrates(calfluxarr,Eval):
    Nsec = 30.0  # Number of secondary ionizations per keV photon.
    x = Eval*1000.0/15.4
    sigH2 = (45.57*(1 - 2.003*x**-0.5 - 4.806/x + 50.577*x**-1.5 - 171.044*x**-2 + 231.608*x**-2.5 - 81.885*x**-3)/Eval**3.5)*1e-24 # Barns -> cm^2

    ZetaH2= N.trapz(calfluxarr[:,2,:]*sigH2*Nsec*Eval,x=Eval)

    xe = Eval*1e3/24.58
    ahe = N.zeros(6)
    ahe[0] = -4.7416
    ahe[1] =  14.8200
    ahe[2] = -30.8678
    ahe[3] =  37.3584
    ahe[4] = -23.4585
    ahe[5] =  5.9133

    seqterm = 0.0
    for i in range(6):
        seqterm = seqterm + ahe[i]/xe**(i/2.0)
    sigHe = 733.0/(Eval**(7.0/2.0))*(1.0+seqterm)*1e-24

    Nsec = 24.0 # Secondary Helium ionizations < H2.
    ZetaHe= N.trapz(calfluxarr[:,2,:]*sigHe*Nsec*Eval,x=Eval)

    return ZetaHe, ZetaH2
############################################################

def main():
# INPUTS:
# User inputs directory name pointing to xrayrt2dcyl and desired X-ray luminosity/template spectrum.
# Possible input is standard template or flaring template.
    #basedir = '/Shire1/cleeves/IDLWorkspace71/xrayrt2dcyl/'
# The code expects the Z-values to be in decreasing order (top down).

    basedir = os.getcwd()
    if basedir[-1] != '/':
        basedir = basedir + '/'

    parser = optparse.OptionParser()
    parser.add_option("-n", "--runname", dest="bsname",type='string', help="Name of UV run after 'dir_' prefix [REQUIRED]. ",default=False)
    parser.add_option("-L", "--Lscale",  dest="L",type="float", help="UV Luminosity relative to TW Hya [REQUIRED].", default=1.0)
    (options, args) = parser.parse_args()

    if options.bsname == False:
        print( 'Need name of input run.')
        parser.print_help()
        sys.exit()

    if options.L == False:
        print( 'Need UV luminosity scale factor.')
        parser.print_help()
        sys.exit()

    bsname = options.bsname
    list,Eval,dust = fixmidplane(basedir,bsname)
    #pdb.set_trace()
    if N.size(Eval) < 2:
        print( 'Not enough wavelengths for interpolation (need at least 2). Exiting.\n')
        sys.exit()
    fluxarr = readUV(list)  # Read in fixed X-ray files.
    calfluxarr,outE = applyspectra(basedir,fluxarr,Eval,L_UV=options.L)  # Scale UV spectrum using model template.
    format4chem(calfluxarr,outE,basedir,bsname,dust) # Calculate ionization rates, output spectra to a file Spectra_..txt.
    print( "Model formatted successfully.\n")

    return

main()



#
# Code checks to see if all desired wavelengths are present, and if not, asks if the user would like to continue and interpolate over missing.
