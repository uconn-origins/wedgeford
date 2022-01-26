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

############################################################
def main():
    start_time = time.time()
    basedir = os.getcwd()
    runsdir = "runs/"

    parser = optparse.OptionParser()
    parser.add_option("-o", "--output", dest="output", help="NAME to save run as [REQUIRED]", metavar="NAME")
    parser.add_option("-d", "--dust", dest="dust", help="pick a single (or list) of DUST values to run at", metavar="DUST", default='e1')
    parser.add_option("-i", "--iofil", dest="iofil", help="Name of 0io.* file.", default='0io.e1')
    parser.add_option("-n", "--nproc", dest="nproc", type="int", help="Number of processes to run at once", default=1)
    parser.add_option("-e", "--environ", dest="environ", help="ENVIRON dir prefix", default='environ/')
    parser.add_option("-r", "--rads", dest="rads", help="Radii to run at", default=None)
    parser.add_option("--min", type=float, dest="min", help="Minimum radius to run at", default=0.0)
    parser.add_option("--max", type=float, dest="max", help="maximum radius to run at", default=1000.0)
    (options, args) = parser.parse_args()

    if options.environ[-1] == '/':
        environdir = basedir + '/'+ options.environ
        environdir = options.environ
    else:
        environdir = basedir + '/'+ options.environ + '/'
        environdir = options.environ + '/'

    try:
        outdir = runsdir + options.output + '/'
    except TypeError:
        print 'Need NAME to save run as'
        parser.print_help()
        sys.exit()

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
    print >>rundc, "./disk_chemistry " + environdir + options.iofil + " " + environdir + "1environ.inp."+options.dust+".$1 > code.out.r$1.out"
    print >>rundc, "mv " + basedir + '/r$1_' + options.dust + '_*.out ' + outdir + options.dust + '/'
    print >>rundc, "mv " + basedir + '/code.out.r$1.out ' + outdir + options.dust + '/'
    print >>rundc, "./streetsweeper.py -o " + outdir + options.dust + '/' + ' -e ' + environdir + "1environ.inp."+options.dust+".$1"
    print >>rundc, "echo R=$1 completed."
    rundc.close()
    os.chmod('run_dc', stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH | stat.S_IROTH)
    
#    rundc = open(basedir+'/run_sweep', 'w')
#    print >>rundc, "./streetsweeper.py -o " + outdir + options.dust + '/' + ' -e ' + environdir + "1environ.inp."+options.dust+".$1"
#    rundc.close()
#    os.chmod('run_sweep', stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IXOTH | stat.S_IROTH)

    #pdb.set_trace()
    # later make this smarter -- check if zone exists, if yes move, if no, log that the zone could not be moved.

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

        inpfiles = envfiles #changed from inpfile += envfiles
        
    radstr = [x.split('e1.')[-1] for x in envfiles]    
    runrad = open(basedir+'/run_rad', 'w')
    for rst in radstr:
        print >>runrad, rst
    runrad.close()
    
## setup for parallel job call
    PROCESSES = options.nproc
    print 'Running with %d processes\n' % PROCESSES    
    commandrun = 'cat run_rad | parallel -j%s ./run_dc' % (PROCESSES)
    os.system(commandrun)
    os.system('cp *.rout '+outdir+'/e1/')  
#    commandrun = 'cat run_rad | parallel -j%s ./run_sweep' % (PROCESSES)
#    os.system(commandrun)

    print "Time elapsed: ", time.time() - start_time, "s"

############################################################
if __name__ == '__main__':
    main()
