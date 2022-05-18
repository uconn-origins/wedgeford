#!/usr/bin/env python

import os, shutil, time, sys
import optparse
import numpy as N
import pdb as pdb

############################################################
def main():

    parser = optparse.OptionParser()
    parser.add_option("-o", "--output", dest="output", help="NAME to save run as [REQUIRED]", metavar="NAME")
    parser.add_option("-e", "--environ", dest="environ", help="environment file", default='environ/')
    (options, args) = parser.parse_args()

    environ = options.environ

    if options.output[-1] != '/':
        outdir = options.output + '/'
    else:
        outdir = options.output

    envfile = environ.split('/')[-1]     # Extract name of environment file.
    dust = outdir.split('/')[-2]         # Extract dust value.
    rad = envfile[envfile.find(dust)+3:] # Extract radius value.

    f = open(environ,'r')
    junk = f.readline()
    junk = f.readline()
    nz = f.readline().strip()
    f.close()
    testfinish = os.path.exists(outdir+"r"+rad+"_"+dust+"_"+nz+".out")
    teststart = os.path.exists(outdir+"r"+rad+"_"+dust+"_1.out")

    if teststart:
        if not testfinish:
            os.system("mv "+outdir+"r"+rad+"_"+dust+"*.out "+outdir+"BadRads/")
#            os.system("echo R="+rad+" finished unsuccessfully.")
            print( "R="+rad+" finished unsuccessfully.")
        else:
#            os.system("echo R="+rad+" finished successfully.")
            print( "R="+rad+" finished successfully.")

    return

############################################################
if __name__ == '__main__':
    main()
