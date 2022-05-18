from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy import interpolate
import sys
sys.path.append('../')
from models import *

outdir = './out/m1_test/'
m0 = load_model(outdir)

prep_he_transfer(m0)
print('ready to do high energy radiation montecarlo!')

do_he_transfer(m0,nt=8)
print('high energy monte carlo done!')

save_heRT(m0)