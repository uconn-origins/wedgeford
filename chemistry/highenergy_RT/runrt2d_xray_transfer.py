#!/usr/bin/env python

# Xray Transfer Wrapper
# Last Modified 10/22/2012
# 1) Edit parameters.txt file - save with unique name
    # 1a) Create run file with that name
# 2) Run xrayrt2dcyl on the parameters file in parallel using Gnu Parallel.
# 3) Save output to a file.

import os
import stat
import sys
import time
# import shutil
import optparse
import numpy as N
from scipy import genfromtxt
import pdb as pdb


################################################################################
def feq(a,b, eps=1e-7):
    if abs(a-b)<eps:
        return 1
    else:
        return 0

################################################################################
def create_inps(basedir, rawinp, Eval, inpname,facm, spacing='linear'):
#    Create inputs with position dependent dust.

    mh = 1.6733e-24 			# m_hydrogen.
    mu = 2.36 					# Mean Molecular Weight

    epsopts = [1e-3,1e-2,1e-1,1.0]  # Epsilon values

    #diskmod = genfromtxt(basedir+'MyDisk/'+rawinp,skip_header=2)
    #f = open(basedir+'MyDisk/'+rawinp, 'r')
    #tmp = f.readline()


    #pdb.set_trace()


    header = """\
Column 1: index
Column 2: r [AU]
Column 3: z [AU]
Column 4: nH [cm^-3]
Column 5: epsilon [dust:gas ratio, rel. to ISM value]
Column 6: sca, dust scattering x-section [cm^-2 per H]
Column 7: abs, dust absorption x-section [cm^-2 per H]
Column 8: g, scattering phase g=<Cos>
Column 9: f_H, atomic H fraction
Column 10: tg, gas temperature
Column 11: f_H2O, water mol. fraction
Column 12: f_OH, OH mol. fraction
 Col1           Col2           Col3           Col4           Col5           Col6           Col7           Col8           Col9          Col10          Col11          Col12 \n"""

    diskmod = genfromtxt(basedir+'MyDisk/'+rawinp,skip_header=2)
    numlines = N.size(diskmod[:,0])
    epsilon = 0.0
    sca = 0.0
    abs = 0.0
    g = 0.0

    old_wd = os.getcwd() # <--- store the initial wd
    os.chdir(basedir)
    try:
        with open('MyDisk/'+inpname+'.txt','w') as outfile:
            outfile.write(header)
            for q in range(0,numlines):
                sca = 0.0
                g = 0.0
                fH = 0.0
                #facm = 1.0
                densityhere = diskmod[q,2]/(mu*mh)*facm #*2.5  # 2.5x massive disk  #*0.68

                if densityhere < 1e-30:
                    densityhere = 0.0

                eps_calc = 100.0/(diskmod[q,2]/diskmod[q,3])
                #if densityhere > 1e-20 and diskmod[q,3] > 1e-30:

                if q > 0 and diskmod[q,1] > diskmod[q-1,1]:
                	densityhere = 0.0  # set top of zone to zero

                #    pdb.set_trace()
                if densityhere < 1e-30:
                    eps_calc = 1.0
                if diskmod[q,3] < 1e-30:
                    eps_calc = 0.0

                local_epsilon = epsopts[N.argmin(N.abs(eps_calc-N.array(epsopts)))]  # Match your value of epsilon to the closest tabulated value of epsilon.

                if N.isnan(local_epsilon):
                    local_epsilon = 1.0

#                if spacing not in ['linear', 'log']:
#                    print 'Invalid spacing value, exiting'
#                    sys.exit(1)
#                else:
#                    absval = read_xray_xsect(basedir,local_epsilon, Eval, spacing)

                absval = read_xray_xsect_prec(basedir, eps_calc, Eval)

                #print local_epsilon,absval
	            #### Do not change epsilon in the input file.  The x-ray opacites already
	            #    incorporate the reduced opacity as a result of dust settling in the surface
	            #    layers, such that the opacities are representative of a gas+dust mix
	            #    (mostly gas).

                line = "%5i%15.7f%15.7f%15.7e%15.7f%15.7e%15.7e%15.7f%15.7e%15.7f%15.7e%15.7e\n" % \
                   (q+1,diskmod[q,0],diskmod[q,1],densityhere,1.0,sca,absval,g,fH,diskmod[q,4],0.0,0.0)
                outfile.write(line)
    finally:
    	os.chdir(old_wd)  # <--- back to the initial wd




################################################################################
def edit_parameters(basedir, dust, Eval, runprefix, inpfile, Rmax=1000., Nphotons=100000.,
    thetamin=55.0, thetamax=90.0, spacing='linear'):

    fin = open(os.path.join(basedir, 'parameters_xray.txt'), 'r')
    outname = 'parameters_%s_E%02d_%s' % (dust,Eval,runprefix)
    fout = open(os.path.join(basedir, outname), 'w')

    runname = '%s_%s_E%02i' % (runprefix, dust, Eval)

    for line in fin:
        if 'Name of run' in line:               # Unique
            fout.write('%s%s\n' % (line[:20], runname))
        elif 'Disk file (no suf)' in line:      # Input file with disk structure
            fout.write('%s%s\n' % (line[:20], inpfile))
        elif 'Wavelength (A)' in line:      # X-ray energy
            fout.write('%s%.2f\n' % (line[:20], Eval))
        elif '1 code unit in AU' in line:       # max size of disk
            fout.write('%s%.1f\n' % (line[:20], Rmax))
        elif 'Num photon packets' in line:      # number of packets
            fout.write('%s%.1f\n' % (line[:20], Nphotons))
        elif 'Source min theta' in line:        # minimum angle (>=50-60?)
            fout.write('%s%.1f\n' % (line[:20], thetamin))
        elif 'Source max theta' in line:        # maximum angle (89?)
            fout.write('%s%.1f\n' % (line[:20], thetamax))
        elif 'Overwrite sca cm2/H' in line:     #
            fout.write('%s%.2e\n' % (line[:20], 6.65E-25))
        elif 'Overwrite epsilon' in line:       # write epsilon
            fout.write('%s%.4f\n' % (line[:20], 1.0))
        elif 'Overwrite abs cm2/H' in line:     # from Tom's table
            fout.write('%s%.4f\n' % (line[:20], 0.0))
        else:
            fout.write(line)

    fout.close()
    fin.close()


################################################################################
def read_xray_xsect(basedir, dust, Eval, spacing):
    fdata = open(os.path.join(basedir, 'bethell_xsects_' + spacing + '.dat'), 'r')
    found_epsilon = False
    xsect = 0.0

    for line in fdata:

        if 'energy' in line:
            continue

        if 'epsilon' in line:
            readeps = float(line.split()[1])
            if feq(readeps, dust):
                found_epsilon = True
                continue

        if not found_epsilon: continue
        words = [float(x) for x in line.split()]
        if Eval == words[0]:
            xsect = words[1]

            return xsect


    fdata.close()

################################################################################
# get more specific absorption opacities by interpolating the bethell tables!
#

def read_xray_xsect_prec(basedir, dust, Eval):
    fdata = open(os.path.join(basedir, 'bethell_xsects_linear.dat'), 'r')

    xsect = 0.0

    alldat = N.zeros((20,4))

    epses = []
    es = []
    epct = -1
    en = 0

    for line in fdata:

        if 'energy' in line:
            continue

        if 'epsilon' in line:
            readeps = float(line.split()[1])
            epses.append(readeps)
            epct += 1
            found_epsilon = True
            continue

        if not found_epsilon: continue
        if en >= 20:
            en = 0
            continue

        words = [float(x) for x in line.split()]
        #print line
        ehere = words[0]
        if epct == 0:
            es.append(ehere)
        xsect = words[1]
        alldat[en,epct] = xsect
        en+=1

    fdata.close()

    dustonly = (alldat[:,3]-alldat[:,0])/(epses[3]-epses[0])
    gasonly = alldat[:,3] - dustonly
    xsect_dust = gasonly + dustonly*dust
    absv = N.interp(Eval,es,xsect_dust)
    return absv


############################################################
def run_idl(rawinp, basedir, dust, Eval, runprefix, numphots, maxR, spacing, thmn, thmx,fmul):
    # Create series of startup files to be piped to IDL via parallel.

    inpname = 'inp_%s_E%02d_%s' % (dust,Eval,runprefix)

    #inpname = 'inp_%s_E%02d_%s' % (dust,Eval,runprefix)

    if spacing:
        create_inps(basedir, rawinp, Eval, inpname, fmul,spacing='log')
        edit_parameters(basedir, dust, Eval, runprefix, inpname, Rmax=maxR, Nphotons=numphots, thetamin=thmn, thetamax=thmx, spacing='log')
    else:
        create_inps(basedir, rawinp, Eval, inpname, fmul,spacing='linear')
        edit_parameters(basedir, dust, Eval, runprefix, inpname, Rmax=maxR, Nphotons=numphots, thetamin=thmn, thetamax=thmx, spacing='linear')

    startup_file = 'startup_%s_E%02d_%s' % (dust,Eval,runprefix)
    outname = 'parameters_%s_E%02d_%s' % (dust,Eval,runprefix)

    fstart = open(os.path.join(basedir, startup_file), 'w')
    fstart.write('uvrt2dcyl, INPUT_FILE=\"%s\", /NO_DEVICE \n' % (outname))
    fstart.write('exit \n')
    fstart.close()

############################################################
def make_folders(basedir,runprefix):
    # Folder where to save output

    folddir1 = os.path.join(basedir, 'Output','dir_'+runprefix)
    if not os.path.exists(folddir1):
        commandrun = 'mkdir %s' % (folddir1)
        os.system(commandrun)

    folddir = os.path.join(basedir, 'Output','dir_'+runprefix,'inputs')
    if not os.path.exists(folddir):
        commandrun = 'mkdir %s' % (folddir)
        os.system(commandrun)

    folddir = os.path.join(basedir, 'Output','dir_'+runprefix,'params')
    if not os.path.exists(folddir):
        commandrun = 'mkdir %s' % (folddir)
        os.system(commandrun)
    return folddir1

############################################################

def main():
    start_time = time.time()

    linearE = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
    logE =    [1.0, 1.2, 1.3, 1.6, 1.8, 2.1, 2.5, 2.9, 3.3, 3.9,  4.5,  5.2,  6.0,  7.0,  8.1,  9.5,  11.0, 12.8, 14.8, 17.2, 20.0]

    parser = optparse.OptionParser()
    parser.add_option("-o", "--output", dest="runoutnm",type='string', help="NAME to save run as [REQUIRED]",default=False)
    parser.add_option("-i", "--rawinp", help="File prefix with disk structure [REQUIRED].", default=None)
    parser.add_option("-d", "--dust", dest="dust", help="Dust settling parameter eps = 100/fg. Options: e1,e0p1,e0p01", metavar="DUST", default='e1')
    parser.add_option("-n", "--nproc", dest="nproc", type="int", help="Number of processes to run at once", default=1)
    parser.add_option("-p", "--photons", dest="numphot", type="float", help="Number of photon packets", default=1000.0)
    parser.add_option("-E", "--energy", help="Specify energy values at which code is run (keV). Default: All 1-20 keV", default=None)
    parser.add_option("-r", "--maxrad", dest="maxrad",type="float", help="Maximum radius, AU",default=500.0)
    # Minimum and maximum angular extent of model, 0deg corresponds to orthogonal to the disk.
    parser.add_option("-t", "--angextent", dest="thetaval",type='string',help="Min&Max Angular Extent of RT, degrees e.g. [36,90]", default='30.0,89.9')
    parser.add_option("-f", "--facmul", dest="facmul",type='float',help="gasmassmultiplier", default=1.0)

    (options, args) = parser.parse_args()
    logE = False
    if options.runoutnm == False:
        print( 'Need NAME to save run as')
        parser.print_help()
        sys.exit()
    runprefix=options.runoutnm
    if options.energy:
        Evals = [float(x) for x in options.energy.split(',')]
    else:
        if logE:
            Evals = logE
        else:
            Evals = linearE

    thetanums = [float(x) for x in options.thetaval.split(',')]

    thmn = thetanums[0]
    thmx = thetanums[1]

    # Location of rt2dcyl folder:
    basedir = os.path.abspath('.')
    if basedir[-1] !=  '/':
        basedir = basedir + '/'
    #basedir = '/Shire1/cleeves/IDLWorkspace71/uvrt2dcyl/'

    PROCESSES = options.nproc
    print( 'Running with %d processes\n' % PROCESSES)

    goflag = 1

    for E in Evals:
        dirtest = basedir+'Output/dir_'+options.runoutnm
        testpath = os.path.exists(dirtest)
        # If path exists (this run is old) *and* you are evaluating wavelength E=1 ask nicely first.
    	# if testpath is True and E == 1:
	#		goflag = input("A previous run exists with this name, continue and overwrite? (Yes = 1, No = 0): ")

        # if not os.path.isfile(basedir+'run_idl'):
        #     print("Making preliminary File")
        #     runidl = open(basedir+'run_idl', 'w')
        #     print(>>runidl, "#!/bin/tcsh")
        #     print(>>runidl, "source /etc/profile.d/idl_setup.csh")
        #     print(>>runidl, "idl $1 > $1.out")
        #     runidl.close()
        #     os.chmod('run_idl', stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH | stat.S_IROTH)

        if goflag == 1:
            folddir1 = os.path.join(basedir, 'Output','dir_'+runprefix)
            run_idl(options.rawinp, basedir, options.dust, E, runprefix,options.numphot, options.maxrad, logE,thmn,thmx,options.facmul)

    if goflag == 1:

        # ------------------------------------------------------------------------------------------
        # Execute the commands!
        folddir = make_folders(basedir,runprefix)
        print( "Beginning Xray Transfer")
        startup_file = 'startup_%s_E??_%s' % (options.dust,options.runoutnm)
        commandrun = 'ls %s | parallel -j%s ./run_idl' % (startup_file,PROCESSES)
        os.system(commandrun)
        print( "Ending Xray Transfer")

    	# ------------------------------------------------------------------------------------------
    	#  Let's move everything and clean up:
        fromhere = os.path.join(basedir, 'Output')
        commandrun = 'mv %s/%s_%s_*.txt %s' % (fromhere,runprefix,options.dust,folddir)
        os.system(commandrun)
        commandrun = 'mv %sparameters_%s_*_%s %s' % (basedir,options.dust, \
            runprefix,folddir+'/params')
        os.system(commandrun)
        commandrun = 'mv %sstartup_%s_*_%s* %s' % (basedir,options.dust, \
            runprefix,folddir)
        os.system(commandrun)
        outdir = os.path.join(basedir, 'Output','dir_'+runprefix,'inputs')
    	#inpfile = 'inp_%s_%s' % (options.dust,runprefix)
        inpfile = 'inp_%s_E%02d_%s' % (options.dust,E,runprefix)
        commandrun = 'cp %sMyDisk/%s.txt %s' % (basedir,inpfile,outdir)
        os.system(commandrun)
        commandrun = 'rm -rf %s/Output/%s_%s*.sav' % (basedir,runprefix,options.dust)
        os.system(commandrun)
    else:
        print( "Exiting.")

    print( "Time elapsed: ", time.time() - start_time, "s")

############################################################
if __name__ == '__main__':
    main()
