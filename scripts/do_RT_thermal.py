from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy import interpolate
import sys
sys.path.append('../')
from models import *


p0 = 'default'
outdir = '/out/m1_test/'
pdef = new_model(p0)
m0 = initialize_model(pdef,outdir=outdir)
save_setup_fig(m0)

prep_thermal_transfer(m0,nphot=100000)
print('ready to do thermal montecarlo!')

do_thermal_transfer(m0,nt=8)
print('thermal monte carlo done!')

save_dustRT(m0)