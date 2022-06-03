#!/usr/bin/env python

import pdb as pdb
import os, time
import matplotlib
import matplotlib.pyplot as plt
import optparse,sys
import numpy as N
from scipy import *
from matplotlib import rc
import fnmatch

############################################################
# Code looks for all 20 wavelengths and replaces erroneous midplane values with value of X-ray
# field above midplane for found files.  Returns successful file list, energy levels and type
# of dust used to generate radiation fields.

def fixmidplane(basedir,bsname):
    numwav = 20
    tag = 'Output/dir_'+bsname+'/'
    endst=''

    lamid = N.zeros(numwav)
    filelist = []
    fileexample = ''
    oldline = ''

    for i in range(numwav):
        fname = "%s_e*_E%02i_inp_*.txt" % (bsname,i+1)
        # print fname
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
        #pdb.set_trace()
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
                    if zv/xv <= 3e-5:   # If point is at midplane!
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
        #sys.exit()
    dust = 'e0p1'
    return filelist, Eval, dust

############################################################

def readXrays(filelist):

    c = 2.99792458e10
    dims = genfromtxt(filelist[0],skip_header=6).shape

    photar = N.zeros((dims[0],dims[1],len(filelist)))
    for i in range(len(filelist)):
        photar[:,:,i] = genfromtxt(filelist[i],skip_header=6)
        photar[:,2,i] = photar[:,2,i]*c

    return photar

############################################################

def applyspectra(basedir,fluxarr,Eval,L_XR=10**29.5,flare=False,altspec=''):

    kev_2_erg  = 1.6021765e-9

    # Read in the X-ray template spectra:
    # Input file is in flux-units, photons/cm^2/s/KeV normalized to L=1 erg/s

    dims = fluxarr.shape
    specdir = basedir+'/pfiles/'

    if altspec == '':

        if flare:
            template = genfromtxt(specdir+'flared_template.txt')
        else:
            template = genfromtxt(specdir+'quiescent_templ_fine.txt')
        readcol = 1

        # Evaluate scale factor from flat spectrum, L=10^44 phot/s.
        lxsum = N.trapz(template[:,readcol]*template[:,0]*kev_2_erg,x=template[:,0])  # Integrate to calculate total luminosity in template spectrum ~ 1.
        scales = template[:,readcol]/lxsum # Normalize so that the area under the spectrum = 1 and that the integrated luminosity = 1 erg/s.
        scales = scales*L_XR*1e-44  # Scale factors at the specified wavelengths under Etemp (keV).
        print( scales*L_XR)


    else:

        template = genfromtxt(specdir+altspec)
        readcol = 1
        scales = template[:,readcol]*1e-44
        print( scales)

    Etemp = template[:,0]

    ct=0
    for E in Eval:
        scx = N.interp(E,Etemp,scales)
        fluxarr[:,2,ct] = fluxarr[:,2,ct]*scx
        ct+=1

    return fluxarr

############################################################
# Eval = Energies in keV RT was evaluted, i.e, 1,3,7,10,15,20
# calfluxarr = Calibrated X-ray spectrum at the above wavelengths.

def formatoutput(calfluxarr,Eval,basedir,bsname):
    dims = calfluxarr.shape

    tag = 'Output/dir_'+bsname+'/'
    fo = open(basedir+tag+'Spectra_'+bsname+'.dat','w')

    datarr = N.zeros((dims[0],2+dims[2]))

    for i in range(dims[0]):
	    spectra = calfluxarr[i,2,:]
	    specstr = ''.join(" %13.6e " % (x) for x in spectra)

	    lineout = "%13.4f%13.4f" % (calfluxarr[i,0,0],calfluxarr[i,1,0])+specstr+"\n"
	    #pdb.set_trace()
	    datarr[i,:] = N.concatenate((calfluxarr[i,0:2,0],calfluxarr[i,2,:]))
	    fo.write(lineout)
    fo.close()

    ZetaHe,ZetaH2 = compionrates(calfluxarr,Eval)

    fo = open(basedir+tag+'ZetaHe_'+bsname+'.dat','w')
    for i in range(dims[0]):
	    lineout = "%13.4f%13.4f%13.4e" % (calfluxarr[i,0,0],calfluxarr[i,1,0],ZetaHe[i])+"\n"
	    fo.write(lineout)
    fo.close()

    fo = open(basedir+tag+'ZetaH2_'+bsname+'.dat','w')
    for i in range(dims[0]):
	    lineout = "%13.4f%13.4f%13.4e" % (calfluxarr[i,0,0],calfluxarr[i,1,0],ZetaH2[i])+"\n"
	    fo.write(lineout)
    fo.close()

    return ZetaHe,ZetaH2
############################################################
# Eval = Energies in keV RT was evaluted, i.e, 1,3,7,10,15,20
# calfluxarr = Calibrated X-ray spectrum at the above wavelengths.

def format4chem(calfluxarr,Eval,basedir,bsname,dust,Lx):
    dims = calfluxarr.shape

    tag = 'Output/dir_'+bsname+'/'
    fo = open(basedir+tag+'Spectra_'+bsname+'.dat','w')

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

    nwv = 20

    rawdata = N.zeros((nr*(2+nwv),nz+1))
    checkflag = N.zeros(nwv)
#   Reformat scaled output
    for g in range(nr):
        d = N.where(radvals == radstore[g])
        d = d[0:nz]
        rawdata[g*(2+nwv),:] = N.concatenate(([radstore[g]],N.zeros(nz)))
        dathere = inputtoToms[d,:]
        if dathere[0,0,2]-dathere[0,1,2] < 0.0:
            dathere = dathere[0,::-1,:]
            rev = True
        else:
        	dathere = dathere[0,:,:]
        	rev = False
        #pdb.set_trace()
        rawdata[g*(2+nwv)+1,:] = N.concatenate(([0],dathere[:,2])) #[0,(transpose(beth(1,d(0)+(nr*nz):d(0)+nz-1+(nr*nz))))]
        for p in range(nwv):
            # Check if wavelength was computed.  Otherwise, interpolate:
            #print p+1
            if p+1 in Eval:
                vh = N.where(Eval == p+1)
                vh=vh[0]
                vh=vh[0]
                photdat = calfluxarr[d,2,vh]
                if rev:
                    rawdata[p+(2+nwv)*g+2,:] = N.concatenate(([p+1],photdat[0,::-1]))
                else:
                    rawdata[p+(2+nwv)*g+2,:] = N.concatenate(([p+1],photdat[0,:]))
            else:
                checkflag[p] = 1
                if p == 0:
                    photdat = calfluxarr[d,2,0]
                elif p == 19:
                    photdat = calfluxarr[d,2,-1]
                else:
                    Earr = N.concatenate(([p+1],Eval))
                    Earr=N.sort(Earr)
                    iE=N.where(Earr == p+1)
                    #pdb.set_trace()
                    Ebel = Earr[iE[0]-1]
                    Eabv = Earr[iE[0]+1]
                    iEa=N.where(Eval == Eabv)
                    iEb=N.where(Eval == Ebel)
                    #print iEa[0],iEb[0],calfluxarr.shape
                    Edist = Eabv-Ebel
                    #pdb.set_trace()
                    photdat = calfluxarr[d,2,iEb[0]]*N.float(abs(Eabv-(p+1)))/Edist + calfluxarr[d,2,iEa[0]]*N.float(abs((p+1)-Ebel))/Edist
                if rev:
                    photdat = photdat[0,::-1]
                else:
                    photdat = photdat[0,:]
                rawdata[p+(2+nwv)*g+2,:] = N.concatenate(([p+1],photdat))
        #pdb.set_trace()
    for p in range(nwv):
        if checkflag[p] == 1:
            print( "Missing wavelength!  Interpolating E=%5i keV." % (p+1)  )


    fo = open(basedir+tag+'AllOut_'+bsname+'.dat','w')
    for g in range(nr):
        fo.write(''.join("%13.5f" % (x) for x in rawdata[g*(2+nwv),:])+"\n")
        fo.write(''.join("%13.5f" % (x) for x in rawdata[g*(2+nwv)+1,:])+"\n")
        for p in range(nwv):
            lineout=rawdata[p+(2+nwv)*g+2,1:-1]
            fo.write('%13.5f' % (rawdata[p+(2+nwv)*g+2,0]) +''.join("%13.5e" % (x) for x in lineout)+"\n")
    fo.close()


    # Now write the pretty version:


    fo = open(basedir+tag+'xray_phot_'+bsname+'_Lx'+('%4.1f'%(N.log10(Lx))).strip()+'.dat','w')
    fo.write('Results from Bethell code with run: '+bsname+'\n')
    for g in range(nr):
        fo.write('Radius(AU) %9.4f' % radstore[g] + '\n')
        fo.write('z(AU)  '+''.join("%10.4f" % (x) for x in rawdata[g*(2+nwv)+1,1:])+'\n')
        fo.write('Energy(keV)               Photons/cm2/s/keV\n')
        for p in range(nwv):
            lineout=rawdata[p+(2+nwv)*g+2,1:]
            fo.write('%7.3f' % (rawdata[p+(2+nwv)*g+2,0]) +''.join("%10.2e" % (x) for x in lineout)+"\n")
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
    parser.add_option("-n", "--runname", dest="bsname",type='string', help="Name of X-ray run after 'dir_' prefix [REQUIRED]. ",default=False)
    parser.add_option("-L", "--Lxray",  dest="LX",type="float", help="Integrated stellar X-ray luminosity [REQUIRED].", default=1e30)
    parser.add_option("-s", "--spectrum", dest="flarespec", help="Spectral shape type. Quiescent TT: False [DEFAULT] or Flaring TT: True.", default=False)
    parser.add_option("-a", "--altsp", dest="altsp", help="Alternative Spectral Model", default='')

    (options, args) = parser.parse_args()

    if options.bsname == False:
        print( 'Need name of input run.')
        parser.print_help()
        sys.exit()

    if options.LX == False:
        print( 'Need stellar X-ray luminosity.')
        parser.print_help()
        sys.exit()

    bsname = options.bsname
    list,Eval,dust = fixmidplane(basedir,bsname)
    if N.size(Eval) < 2:
        print( 'Not enough wavelengths for interpolation (need at least 2). Exiting.\n')
        sys.exit()
    fluxarr = readXrays(list)  # Read in fixed X-ray files.
    calfluxarr = applyspectra(basedir,fluxarr,Eval,L_XR=options.LX,flare=options.flarespec,altspec=options.altsp)  # Scale x-ray spectrum using model templates.  Normalized to have Lx.
    ZetaHe,ZetaH2 = formatoutput(calfluxarr,Eval,basedir,bsname) # Calculate ionization rates, output spectra to a file Spectra_..txt.
    format4chem(calfluxarr,Eval,basedir,bsname,dust,options.LX) # Calculate ionization rates, output spectra to a file Spectra_..txt.
    print( "Model formatted successfully.\n")

    return

main()



#
# Code checks to see if all desired wavelengths are present, and if not, asks if the user would like to continue and interpolate over missing.
