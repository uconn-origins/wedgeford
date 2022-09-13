#!/usr/bin/env python

## 1) delete the calculated photodissociaton files from the last run
## 2) extract the radius/dust values from the environ file
## 3) run ./disk_chemistry
## 4) create radii.<dust>.<radius> file
## 5) move <radius>_<dust>_<zone>.out files to runs/

import os, shutil, time, sys, subprocess
import optparse
import numpy as N
import logging
import multiprocessing as MP
import pdb as pdb
import stat
import fnmatch

############################################################
def main():
    start_time = time.time()
    basedir = os.getcwd()
    runsdir = "runs/"

    parser = optparse.OptionParser()
    parser.add_option("-o", "--output", dest="output", help="NAME to save run as [REQUIRED]", metavar="NAME")
    parser.add_option("-d", "--dust", dest="dust", help="pick a single (or list) of DUST values to run at", metavar="DUST", default='e1')
    parser.add_option("-i", "--iofil", dest="iofil", help="Name of *end* of 0iofile.", default='.pl')
    parser.add_option("-n", "--nproc", dest="nproc", type="int", help="Number of processes to run at once", default=1)
    parser.add_option("-e", "--environroot", dest="environroot", help="ENVIRON dir prefix", default='environ/')
    parser.add_option("-r", "--rads", dest="rads", help="Radii to run at", default=None)
    parser.add_option("--min", type=float, dest="min", help="Minimum radius to run at", default=0.0)
    parser.add_option("--max", type=float, dest="max", help="Maximum radius to run at", default=100.0)
    parser.add_option("--minphi", type=float, dest="minphi", help="Minimum phi (az) to run at", default=0.0)
    parser.add_option("--maxphi", type=float, dest="maxphi", help="Maximum phi (az) to run at", default=180.0)
    parser.add_option("--crtag", dest="crtag", help="CR ionization rate", default='w98')
    parser.add_option("-t", type=float, dest="timeout", help="Time-out limit", default=10000)    
    (options, args) = parser.parse_args()

    

    if options.environroot[-1] == '/':
        environroot = basedir + '/'+ options.environroot
        environroot = options.environroot
    else:
        environroot = basedir + '/'+ options.environroot + '/'
        environroot = options.environroot + '/'

    try:
        outdir = runsdir + options.output + '/'
    except TypeError:
        print 'Need NAME to save run as'
        #parser.print_help()
        #sys.exit()

    listnm = []
    listvl = []
    minphi = options.minphi
    maxphi = options.maxphi

    for file in os.listdir(environroot):
        if fnmatch.fnmatch(file,"planet*.out"):
            foundfile = file 
            i_b = file.index('loT_')
            i_e = file.index('.out')
            listnm.append(file[i_b+4:i_e])
            listvl.append(N.float(file[i_b+4:i_e]))

    i = N.argsort(listvl)
    listnm = [listnm[t] for t in i]
    listvl = [listvl[t] for t in i]
    
    angles = [listnm[i] for i in N.where((N.array(listvl)>=minphi)*(N.array(listvl)<=maxphi)==True)[0]]
    angles = angles[::-1]
    
    #print 'hello'    
    crtag = options.crtag    
    
    for b in angles:
        environdir = environroot+'slice__'+b+'_'+crtag+'/'
        
        i_rb = environroot[8:].index('net_')
        i_re = environroot[8:].index('AU') 
        i_ab = environroot[8:].index('dot')
       # i_ae = environroot[8:].index('_w') 
        i_ae = environroot[8:].index('_loT') 
        
        iofile = '0io.'+environroot[8:][i_rb+4:i_re]+environroot[8:][i_ab+3:i_ae]+options.iofil

        os.system('cp -f '+environroot+iofile+' '+environdir)

        outdir = runsdir+options.output + '/' +'slice__'+b+'_'+crtag+'/'
        
        dustvals = options.dust.split(',')
        if 'e1' not in dustvals and 'e0p1' not in dustvals and 'e0p01' not in dustvals:
            print 'Need appropriate DUST value to run at'
            parser.print_help()
            sys.exit() 

        testpath = os.path.exists(outdir + options.dust + '/')
        if not testpath:    
            command = 'mkdir -p ' + outdir + options.dust
            command = 'mkdir -p ' + outdir + options.dust + '/BadRads'
            os.system(command)

        os.system('cp 3abunds.inp '+outdir)
        os.system('cp 5flags.inp '+outdir)
        os.system('cp 2times.inp '+outdir)
        os.system('cp 6grainbe.inp '+outdir)        
        rundc = open(basedir+'/run_dc', 'w')
        print >>rundc, "./disk_chemistry " + environdir + iofile + " " + environdir + "1environ.inp."+options.dust+".$1 > code.out.r$1.out"
        print >>rundc, "mv " + basedir + '/r$1_' + options.dust + '_*.out ' + outdir + options.dust + '/'
        print >>rundc, "mv " + basedir + '/code.out.r$1.out ' + outdir + options.dust + '/'
        print >>rundc, "./streetsweeper.py -o " + outdir + options.dust + '/' + ' -e ' + environdir + "1environ.inp."+options.dust+".$1"
        print >>rundc, "echo R=$1 completed."
        rundc.close()
        os.chmod('run_dc', stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH | stat.S_IROTH)
    
## Delete the calculated photo-xsect files so that the correct one is created at runtime
        for file in os.listdir(basedir + '/xsect'):
            if '.calc' in file: 
                os.remove(basedir + '/xsect/' + file)

# Only look at environ files, sort by radius and reverse the sort so that the smallest (slowest) radii are last
        inpfiles = []
        rads = []
        for dust in dustvals:
            envfiles = [x for x in os.listdir(environdir) if '1environ' in x and dust in x]
            envfiles = N.array([x for x in envfiles if options.min <= float(x.split('.')[3]+'.'+x.split('.')[4]) <= options.max])
            envnums = N.array([float(x.split('.')[-2] + '.' + x.split('.')[-1]) for x in envfiles])
            envfiles = envfiles[envnums.argsort()][::-1]
            
            # Find closest radii to entered values
            if options.rads is not None:
                for val in options.rads.split(','):
                    rads.append(str(envnums[abs(envnums-float(val)).argmin()]))
                envfiles = [x for x in envfiles if (x.split('.')[3]+'.'+x.split('.')[4]) in rads]

            inpfiles += envfiles
            
        radstr = [x.split('e1.')[-1] for x in envfiles]    
        runrad = open(basedir+'/run_rad', 'w')
        for rst in radstr:
            print >>runrad, rst
        runrad.close()
        
## setup for parallel job call
        PROCESSES = options.nproc
        print 'Running with %d processes\n' % PROCESSES    
        commandrun = 'cat run_rad | parallel -j%s --timeout %f ./run_dc' % (PROCESSES,options.timeout)
        os.system(commandrun)
#    commandrun = 'cat run_rad | parallel -j%s ./run_sweep' % (PROCESSES)
#    os.system(commandrun)

        print "Time elapsed: ", time.time() - start_time, "s"

    return
############################################################
if __name__ == '__main__':
    main()
