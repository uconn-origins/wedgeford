from . make_model import model_
from . outputs import *
from . units import *
from . write_radmc_files import *
from . make_plots import *
from . dust import *
from . observe import *
from . import_fargo import *
from .templates import param_library
import os
import numpy as np
import pickle

from pylab import *

from pathlib import Path

p = Path().resolve()
parent_dir = str(p).split('wedgeford')[0] + 'wedgeford/'
models_dir = parent_dir+'/models/'

np.seterr(all='ignore')


def new_model(key='default'):
    par_dict = param_library.library[key]
    return par_dict.copy() #important if you are making several models at once

def initialize_model(params,outdir='/m1_test/'):
    model1 = model_(params,outdir=outdir)
    # define dictionary
    os.chdir(model1.outdir)
    pickle.dump(params,open('pars.pkl','wb'))
    file_list = ['amr_grid.inp', 'wavelength_micron.inp','stars.inp']
    func_list = [write_grid, write_wavelength,write_star]
    for file, func in zip(file_list,func_list):
        if os.path.exists(file) != True:
            func(model1)
    output1 = out(model1)
    print('Assigned model directory:' + model1.outdir)
    return output1

def load_model(outdir='/m1_test/'):
    os.chdir(parent_dir+outdir)
    params = pickle.load( open( 'pars.pkl', "rb" ) )
    model1 = model_(params,outdir=outdir)
    print('Loading from model directory:' + model1.outdir)
    file_list = ['amr_grid.inp', 'wavelength_micron.inp','stars.inp']
    func_list = [write_grid, write_wavelength,write_star]
    for file, func in zip(file_list,func_list):
        if os.path.exists(file) != True:
            func(model1)
    output1 = out(model1)
    return output1