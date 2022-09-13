#!/usr/bin/env python

# RT Transfer Wrapper
# Last Modified 2/26/14, LIC
# 1) Edit parameters.txt file - save with unique name
    # 1a) Create run file with that name
# 2) Run rt2dcyl on the parameters file in parallel using Gnu Parallel.
# 3) Save output to a file.
# To run in terminal: 
# './runrt2d_uv_transfer.py -o'+s+curr_name_uv+s+'-i'+s+output_filename1+s+'-n 15 -p 300000 -r 200 -E 1,2,3,4,5,7,9 -k pfiles/'+dir_name+'/'

import os
import stat
import sys
import time
import optparse
import pdb as pdb
from scipy import genfromtxt
import numpy as N

################################################################################
def feq(a,b, eps=1e-7):
    if abs(a-b)<eps:
        return 1
    else:
        return 0

################################################################################
# Create Parameters Files:
#
def edit_parameters(basedir, dust, Eval, runprefix, inpfilename, linearlam, Rmax=1000., Nphotons=100000.,
    thetamin=36.0, thetamax=89.5):

    fin = open(os.path.join(basedir, 'parameters_uv.txt'), 'r')
    outname = 'parameters_%s_E%02d_%s' % (dust,Eval,runprefix)
    fout = open(os.path.join(basedir, outname), 'w')

    runname = '%s_%s_E%.1f' % (runprefix, dust, Eval)

    # Check if doing Lyman Alpha:
    if N.abs(linearlam[int(Eval-1)]-1215.668) < 10.0:
    	lymflag = 1.0
    else:
    	lymflag = 0.0

    for line in fin:
        if 'Name of run' in line:               # Unique
            fout.write('%s%s\n' % (line[:20], runname))
        elif 'Disk file (no suf)' in line:      # Input file with disk structure
            fout.write('%s%s\n' % (line[:20], inpfilename))
        elif 'Wavelength (A)' in line:      # Wavelength
            fout.write('%s%.2f\n' % (line[:20], linearlam[int(Eval-1)]))
        elif '1 code unit in AU' in line:       # max size of disk
            fout.write('%s%.1f\n' % (line[:20], Rmax))
        elif 'Num photon packets' in line:      # number of packets
            fout.write('%s%.1f\n' % (line[:20], Nphotons))
        elif 'Source min theta' in line:        # minimum angle (>=50-60?)
            fout.write('%s%.1f\n' % (line[:20], thetamin))
        elif 'Source max theta' in line:        # maximum angle (89?)
            fout.write('%s%.1f\n' % (line[:20], thetamax))
        elif 'Overwrite abs cm2/H' in line:     # from Tom's table
            fout.write('%s%.1f\n' % (line[:20], 0.0))
        elif 'Overwrite sca cm2/H' in line:
            fout.write('%s%.1f\n' % (line[:20], 0.0))
        elif 'Overwrite g=<cos>' in line:
            fout.write('%s%.1f\n' % (line[:20], 0.0))
        elif 'Overwrite epsilon' in line:
            fout.write('%s%.1f\n' % (line[:20], 0.0))
        elif 'Doing Ly-a?' in line:
            fout.write('%s%i\n' % (line[:20], lymflag))
        else:
            fout.write(line)

    fout.close()
    fin.close()

################################################################################
# Read Albedo files for each type of dust (# read from model file header).
def read_albedo(basedir, albedodir, numberdust):

    line = "albedo01.dat"
    filedest = os.path.join(basedir, albedodir, line)
    tmp = genfromtxt(filedest,skip_header=2)
    print( filedest)
    print( "fileshape:",tmp.shape)
    nlines = N.size(tmp[:,0])

    albedoarr = N.zeros((nlines,6,numberdust))
    for j in range(0,numberdust):
    	line = "albedo%02i.dat" % (j+1)
    	filedest = os.path.join(basedir, albedodir, line)
    	print( "Reading: ",filedest)
    	albedoarr[:,:,j] = genfromtxt(filedest,skip_header=2)

    return albedoarr

################################################################################
def create_inps(basedir,rawinp,inpname,dust,albedoarr,fg,lamval,dusthole,gashole):
#    Create inputs with position dependent dust.

    mh = 1.6733e-24 			# m_hydrogen.
    mu = 2.3 					# Mean Molecular Weight

    diskmod = genfromtxt(basedir+'MyDisk/'+rawinp,skip_header=2)
    f = open(basedir+'MyDisk/'+rawinp, 'r')
    tmp = f.readline()
    itmp = tmp.index(' types of dust')
    numberdust = int(tmp[itmp-5:itmp])

    print( 'Number of dust species used:',numberdust)
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
    opacity_lam = N.zeros((numberdust)) + 1e4
    numlines = N.size(diskmod[:,0])
    epsilon = float(dust.replace('e','').replace('p','.'))
    sca = 0.0
    abs = 0.0
    g = 0.0

    lamstore = albedoarr[:,0,0]*1e4
    ind = N.argmin(N.abs(lamstore-lamval))

    #scaopa = /fg*(mu*mh)
    #absopa = albedoarr[ind,2,:]/fg*(mu*mh)
    #g_cos  = albedoarr[ind,5,:]

    old_wd = os.getcwd() # <--- store the initial wd
    os.chdir(basedir)
    try:
        with open('MyDisk/'+inpname+'.txt','w') as outfile:
            outfile.write(header)
            for q in range(0,numlines):
                #print diskmod[q,4:4+numberdust]
                #pdb.set_trace()
                #sca = sum(diskmod[q,4:4+numberdust]*scaopa)
                #abs = sum(diskmod[q,4:4+numberdust]*absopa)
                #g = sum(diskmod[q,4:4+numberdust]*g_cos)
                fg = diskmod[q,2]/diskmod[q,3]
                if diskmod[q,2] < 1e-25:
                    fg = 100.0
               # pdb.set_trace()
                sca = sum(diskmod[q,5:5+numberdust]*albedoarr[ind,3,:]/fg*(mu*mh)) # converts to per H2
                abs = sum(diskmod[q,5:5+numberdust]*albedoarr[ind,2,:]/fg*(mu*mh)) # converts to per H2
                g = sum(diskmod[q,5:5+numberdust]*albedoarr[ind,5,:])

                #sca = sum([0,1,0,0]*scaopa)
                #abs = sum([0,1,0,0]*absopa)
                #g = sum([0,1,0,0]*g_cos)

                #print diskmod[q,0],diskmod[q,1],diskmod[q,4+numberdust],diskmod[q,5+numberdust]
                #pdb.set_trace()

                if diskmod[q,0] < dusthole:
                	sca = 0.0
                	abs = 0.0
                	g = 0.0

                #if diskmod[q,0] < gashole:
                #	fH = 0.0
                dustdensityhere = diskmod[q,3]
                #fg = 100.0
                #gasden_fromdust = dustdensityhere*fg/(mu*mh)         # Gas density from dust profile
                #densityhere = gasden_fromdust

                densityhere = diskmod[q,2]/(mu*mh)                   # Actual gas density

                if densityhere < 1e5:
                	densityhere = 0.0
                if q > 0 and diskmod[q,1] > diskmod[q-1,1]:
                	densityhere = 0
                	#fH = 0.0

                # model sca wanted in dust scat cross section per H.
                # albedo file is dust sca cross / gram dust (kappa_d)
                # for rho_g and rho_d  ideally want rho_d * kappa_d * dz
                # model is n_g * kappa_H * dz
                # what is kappa_H?
                # rho_d * kappa_d = n_g * kappa_H
                # kappa_H = rho_d/n_g * kappa_d = rho_d/rho_g * mh*mu * kappa_d = mu_mg/ fg * kappa_d

                #eps_calc = 100.0/(diskmod[q,2]/diskmod[q,3])
                #if N.isnan(eps_calc):
                #   eps_calc = 1.0

                fH = diskmod[q,5+numberdust]
                line = "%5i%15.7f%15.7f%15.7e%15.7f%15.7e%15.7e%15.7f%15.7e%15.7f%15.7e%15.7e\n" % \
                   (q+1,diskmod[q,0],diskmod[q,1],densityhere,1.0,sca,abs,g,fH,diskmod[q,4],0.0,0.0)
                outfile.write(line)
    finally:
    	os.chdir(old_wd)  # <--- back to the initial wd


 ################################################################################

def run_idl(rawinp, basedir, dust, Eval, runprefix, numphots, maxR, opacdir, fg, dusthole, gashole, linearlam, thmn, thmx):
    # Create series of startup files to be piped to IDL via parallel.

    # How many dust species?
    f = open(basedir+'MyDisk/'+rawinp, 'r')
    tmp = f.readline()
   # numberdust = int(tmp[61:66])
    itmp = tmp.index(' types of dust')
    numberdust = int(tmp[itmp-5:itmp])
    print( numberdust)
    # Create input model file names for each wavelength:
    inpname = 'inp_%s_E%02d_%s' % (dust,Eval,runprefix)
    albedodir = opacdir
    albedoarr = read_albedo(basedir, albedodir, numberdust)

    create_inps(basedir,rawinp,inpname,dust,albedoarr,fg,linearlam[int(Eval-1)],dusthole,gashole)
    edit_parameters(basedir, dust, Eval, runprefix, inpname,linearlam, Rmax=maxR, Nphotons=numphots,thetamin=thmn, thetamax=thmx)

    startup_file = 'startup_%s_E%02d_%s' % (dust,Eval,runprefix)
    outname = 'parameters_%s_E%02d_%s' % (dust,Eval,runprefix)

    # Create IDL scripts that run the models:
    fstart = open(os.path.join(basedir, startup_file), 'w')
    fstart.write('uvrt2dcyl, INPUT_FILE=\"%s\", /NO_DEVICE \n' % (outname))
    fstart.write('exit \n')
    fstart.close()

############################################################
def make_folders(basedir,runprefix):
    # Folder where to save output
    def folder_exists(folder):
        if not os.path.exists(folder):
            commandrun = 'mkdir %s' % (folder)
            os.system(commandrun)
    
    folddir = os.path.join(basedir,'Output','dir_'+runprefix)
    folder_exists(folddir)
    folddir1 = os.path.join(folddir,'inputs')
    folddir2 = os.path.join(folddir, 'params')
    
    folder_exists(folddir1)
    folder_exists(folddir2)
    
    #folddir1 = os.path.join(basedir, 'Output','dir_'+runprefix)
    #if not os.path.exists(folddir1):
       # commandrun = 'mkdir %s' % (folddir1)
        #os.system(commandrun)

    #folddir = os.path.join(basedir, 'Output','dir_'+runprefix,'inputs')
   # if not os.path.exists(folddir):
        #commandrun = 'mkdir %s' % (folddir)
        #os.system(commandrun)

    #folddir = os.path.join(basedir, 'Output','dir_'+runprefix,'params')
    #if not os.path.exists(folddir):
        #commandrun = 'mkdir %s' % (folddir)
        #os.system(commandrun)
    
    return folddir2

############################################################

def main():

    # Begin timer:
    start_time = time.time()

    # Wavelengths of transfer are hardcorded here and in the RT compilation script later. If this is
    # changed, the corresponding values in wvlnfile also must be changed.
    #linearE = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    #linearlam = [930.0,1000.0,1100.0,1216.0,1300.0,1400.0,1600.0,1800.0,2000.0]
    
    #loads the wavelengths in the wvlnfile txt
    linearE = list(np.loadtxt('wvlnfile.txt')[:,1])
    linearlam = list(np.loadtxt('wvlnfile_uv.txt')[:,1])

    # General configuration parameters.
    parser = optparse.OptionParser()
    parser.add_option("-o", "--output", dest="runoutnm",type='string', help="NAME to save run as [REQUIRED]",default=False)
    parser.add_option("-i", "--rawinp", dest="rawinp", help="File off-which to make inps", default='none')
    parser.add_option("-d", "--dust", dest="dust", help="dust", metavar="DUST", default='e1')
    parser.add_option("-n", "--nproc", dest="nproc", type="int", help="Number of processes to run at once", default=1)
    parser.add_option("-p", "--photons", dest="numphot", type="float", help="Number of photon packets", default=100.0)
    parser.add_option("-E", "--energy", help="Energy values to run at", default=None)
    parser.add_option("-r", "--maxrad", dest="maxrad",type="float", help="Maximum radius, AU",default=1000.0)
    parser.add_option('--dusty', dest='dustvar', help='Use position dependent dust, if false then constant.', default=True)
    parser.add_option("-k", "--opac", dest="opacdir", help="Albedo-format directory with opacities", default='none')
    parser.add_option("-g", "--gIId", dest="gas2dust", type="float", help="Gas to Dus ratio", default=100.0)
    parser.add_option("-f", "--dsthl", dest="dusthole", type="float", help="DustOpacityHole? (AU)", default=0.0)
    parser.add_option("-a", "--gashl", dest="gashole", type="float", help="GasHole? (AU)", default=0.0)
    # Minimum and maximum angular extent of model, 0deg corresponds to orthagonal to the disk.
    parser.add_option("-t", "--angextent", dest="thetaval",type='string',help="Min&Max Angular Extent of RT, degrees e.g. [36,90]", default='35.0,89.9')
    (options, args) = parser.parse_args()

    if options.runoutnm == False:
        print( 'Need NAME to save run as.  Here are all of your options.')
        parser.print_help()
        sys.exit()
    runprefix=options.runoutnm

#  If only running a subset of energy values (wavelengths).
    if options.energy:
        Evals = [float(x) for x in options.energy.split(',')]
    else:
        Evals = linearE

    thetanums = [float(x) for x in options.thetaval.split(',')]

    thmn = thetanums[0]
    thmx = thetanums[1]

    # Location of rt2dcyl folder:
    basedir = os.path.abspath('.')
    if basedir[-1] !=  '/':
        basedir = basedir + '/'
        
    ## setup for parallel job call
    PROCESSES = options.nproc
    print( 'Running with %d processes\n' % PROCESSES)

    goflag = 1

    for E in Evals:
        dirtest = basedir+'Output/dir_'+options.runoutnm

        testpath = os.path.exists(dirtest)

        # If path exists (this run is old) *and* you are evaluating wavelength E=1 ask nicely first.
        if testpath is True and E == 1:
            goflag = input("A previous run exists with this name, continue and overwrite? (Yes = 1, No = 0): ")

        if not os.path.isfile(basedir+'run_idl'):
            print( "Making preliminary File")
            runidl = open(basedir+'run_idl', 'w')
            print("#!/bin/tcsh",file=runidl)
            print("source /etc/profile.d/idl_setup.csh",file=runidl)
            print("idl $1 > $1.out",file=runidl)
            runidl.close()
            os.chmod('run_idl', stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH | stat.S_IROTH)

        if goflag == 1:
            # Remove any old run-data and then startup code.
            folddir1 = os.path.join(basedir, 'Output','dir_'+runprefix)
            # Create inputs:
            run_idl(options.rawinp, basedir, options.dust, E, runprefix,options.numphot, options.maxrad, \
                options.opacdir,options.gas2dust,options.dusthole,options.gashole,linearlam,thmn,thmx)

    if goflag == 1:

        # ------------------------------------------------------------------------------------------
        # Execute the commands!
        folddir = make_folders(basedir,runprefix)
        print( "Beginning UV Transfer")
        startup_file = 'startup_%s_E??_%s' % (options.dust,options.runoutnm)
        commandrun = 'ls %s | parallel -j%s ./run_idl' % (startup_file,PROCESSES)
        os.system(commandrun)
        print( "Ending UV Transfer")

        # ------------------------------------------------------------------------------------------
        # Let's move everything and clean up:
        fromhere = os.path.join(basedir, 'Output')
        commandrun = 'mv %s/%s_%s_*.txt %s' % (fromhere,runprefix,options.dust,folddir)
        os.system(commandrun)

        # House cleaning!
        commandrun = 'mv %s/parameters_%s_*_%s %s' % (basedir,options.dust, \
            runprefix,folddir+'/params')
        os.system(commandrun)
        commandrun = 'mv %s/startup_%s_*_%s* %s' % (basedir,options.dust, \
            runprefix,folddir)
        os.system(commandrun)
        outdir = os.path.join(basedir, 'Output','dir_'+runprefix,'inputs')
        commandrun = 'cp %sMyDisk/inp_%s*_%s.txt %s' % (basedir,options.dust,runprefix,outdir)
        os.system(commandrun)
        commandrun = 'rm -rf %s/Output/%s_%s*.sav' % (basedir,runprefix,options.dust)
        os.system(commandrun)
        # ------------------------------------------------------------------------------------------

    else:
    		print( "Exiting.")

    print( "Time elapsed: ", time.time() - start_time, "seconds" )


############################################################
if __name__ == '__main__':
    main()
