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

############################################################
def before_zone(envfile, outdir):
## for each environ file, extract the dust and radius values
    words = envfile.split(".")
    dust = words[2]
    if len(words) == 5:
        radius = words[3] + '.' + words[4]
    else:
        radius = words[3]

## Create directory to save run in if it doesn't already exist
    try:
        os.mkdir(outdir)
    except OSError:
        pass

    try:
        os.mkdir(outdir + dust + '/', 0755)
    except OSError:
        pass

    return radius, dust

############################################################
def after_zone(radius, dust, basedir, outdir, nzones, retcode):
## move the output files to the runs directory if finished, otherwise delete them
    zone = 1
    while zone <= nzones:
        if not retcode:     # radii completed
            try:
                shutil.move(basedir + '/r' + radius + '_' + dust + '_' + str(zone) + '.out', outdir + dust + '/')
            except shutil.Error:    # File already exists
                shutil.copy(basedir + '/r' + radius + '_' + dust + '_' + str(zone) + '.out', outdir + dust + '/')
                os.remove(basedir + '/r' + radius + '_' + dust + '_' + str(zone) + '.out')
            except:
                logging.error('R = ' + radius + ', Z = ' + str(zone) + ' was not moved')
        else:       # radii didn't copmlete
            try:
                os.remove(basedir + '/r' + radius + '_' + dust + '_' + str(zone) + '.out')
            except:
                logging.error('R = ' + radius + ' last zone completed was %d' % (zone-1))
                zone = nzones + 1
        zone += 1

# If selfshield.test files are created, move them to the appropriate directory
    try:
        shutil.move(basedir + '/selfshield_' + dust + '_' + radius + '.test', outdir + dust + '/')
    except IOError:
        pass

#    return 'R = %s, %s Done' % (radius, dust)
    print 'R = %s, %s Done' % (radius, dust)

############################################################
def run_zone(inpfile, basedir, environdir, runsdir, outdir, nzones, iofilenm, rnrate):

    # Set timeout for each zone to 3 hours
    start = time.time()
    TIMEOUT = 3.0 * 60 * 60


    #shutil.copyfile(basedir+'/5flags_new.inp',basedir+'/5flags.inp')
    # return
    radius, dust = before_zone(inpfile, outdir)

    iofile =  environdir + iofilenm
    envfile = environdir + '1environ.inp.' + dust + '.' + radius
    cmd = basedir + '/disk_chemistry ' + iofile + ' ' + envfile

    ## run ./disk_chemistry and save the output to a code.out file
    Log = open(outdir + dust + '/code.out.r' + radius + '.' + dust, 'w')
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

#    while (p.poll() == None) and (time.time() < start + TIMEOUT):
#        pass

    # None if timeout, 0 if finished, 1 if killed
    retcode = p.poll()
#    if retcode == None:
#        print 'kill'
#        os.kill(p.pid, 9)
#        retcode = 1

    stdoutdata, stderrdata = p.communicate()
    Log.write(stdoutdata)
    Log.close()
    
#    return after_zone(radius, dust, basedir, outdir, nzones, retcode)
    after_zone(radius, dust, basedir, outdir, nzones, retcode)

############################################################
def main():
    start_time = time.time()

# General configuration parameters.  Probably don't need to be changed once set
    # nzones = 60
    t = os.getcwd()

    basedir = os.getcwd()

    runsdir = basedir + "/runs/"
    runsdir = "runs/"

    parser = optparse.OptionParser()
    parser.add_option("-o", "--output", dest="output", help="NAME to save run as [REQUIRED]", metavar="NAME")
    parser.add_option("-d", "--dust", dest="dust", help="pick a single (or list) of DUST values to run at", metavar="DUST", default='e1')
    parser.add_option("-i", "--iofil", dest="iofil", help="Name of 0io.* file.", default='0io.e1')
    parser.add_option("-n", "--nproc", dest="nproc", type="int", help="Number of processes to run at once", default=1)
    parser.add_option("-e", "--environ", dest="environ", help="ENVIRON dir prefix", default='environ/')
    parser.add_option("-t", "--timeout", type=int, dest="timeout", help="TIMEOUT (s) for entire model", default=86400)
    parser.add_option("-r", "--rads", dest="rads", help="Radii to run at", default=None)
    parser.add_option("--min", type=float, dest="min", help="Minimum radius to run at", default=0.0)
    parser.add_option("--max", type=float, dest="max", help="maximum radius to run at", default=400.0)
    parser.add_option("-z","--nzones", type=int, dest="nzones",help="Number of zones per radius.",default=45)
    parser.add_option("-u","--radn", type=float, dest="rnrate", help="radionuclide ionization rate", default=1e-18)
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


    # print 'opening flags'
    # Modify the 5flags.inp file to contain RN ionrate
    # print basedir+'5flags.inp'
    #shutil.copyfile(,)
    os.system("rm -f "+basedir+'/5flags_prev.inp')
    os.system("rm -f "+basedir+'/5flags_new.inp')
    os.system("cp "+basedir+'/5flags.inp '+basedir+'/5flags_prev.inp')
    ffo = open(basedir+'/5flags_new.inp','w')
    ff = open(basedir+'/5flags.inp', 'r')
    for line in ff:    
        # line = ff.readline()
        tmp = line.split(" ")
        # tmp = tmp[0]
        # print tmp[0]+'r'
        if tmp[0] != 'RadNuc':
            ffo.write(line)
        else:
            temp = 'RadNuc        = %8.2E\n' % (options.rnrate)
            ffo.write(temp)
    # print 'closing flags'
    ff.close()
    ffo.close()
    os.system("cp "+basedir+'/5flags_new.inp '+basedir+'/5flags.inp')

    nzones = options.nzones
    dustvals = options.dust.split(',')
    if 'e1' not in dustvals and 'e0p1' not in dustvals and 'e0p01' not in dustvals:
        print 'Need appropriate DUST value to run at'
        parser.print_help()
        sys.exit()

# Logging
    logging.basicConfig(level=logging.ERROR, filename=options.output + '.errorlog', 
                        filemode='w', format='%(asctime)s %(levelname)s: %(message)s', 
                        datefmt='%Y-%m-%d %H:%M:%S')

## Delete the calculated photo-xsect files so that the correct one is created at runtime
    for file in os.listdir(basedir + '/xsect'):
        if '.calc' in file: 
            os.remove(basedir + '/xsect/' + file)

# only look at environ files, sort by radius and reverse the sort so that the smallest (slowest) radii are last
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

## setup for parallel job call
    PROCESSES = options.nproc
    print 'Running with %d processes\n' % PROCESSES
    pool = MP.Pool(processes=PROCESSES)

## run the models
    results = []
    try:
        for inpfile in inpfiles:
            result = pool.apply_async(run_zone, args=(inpfile, basedir, environdir, runsdir, outdir, nzones, options.iofil,options.rnrate))
            results.append(result)
        pool.close()
        pool.join()
    except KeyboardInterrupt:
        print 'CTRL-C pressed, exiting...'
        pool.terminate()
    except MP.TimeoutError:
        print 'Model did not finish and timed out after %d s' % options.timeout

    print "Time elapsed: ", time.time() - start_time, "s"

############################################################
if __name__ == '__main__':
    main()
