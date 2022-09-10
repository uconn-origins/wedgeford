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
    """ generates a copy of a parameter dictionary listed in the templates directory in param_library.py
    Parameters:
    -----------
    key: name of the default model dictionary being used in the template directory
    
    Returns: a copy of the dictionary that can be modified before inititalizing a new model
    """
    par_dict = param_library.library[key]
    return par_dict.copy() #important if you are making several models at once

def initialize_model(params,outdir='/m1_test/'):
    """ initializes a new model in a directory based on input parameters 
    writes a grid, wavelength, and stars input file if there is not one already
    pkls the list of parameters the model was initialized with
    
    Parameters:
    ----------
    params: dict, dictionary of parameters to pass into the model_ class
    
    outdir: str, directory where model will be initialized
    
    Returns: output class object see outputs.py
    
    """
    model1 = model_(params,outdir=outdir)
    os.chdir(model1.outdir)
    pickle.dump(params,open('pars.pkl','wb'))
    file_list = ['amr_grid.inp', 'wavelength_micron.inp']
    func_list = [write_grid, write_wavelength]
    for file, func in zip(file_list,func_list):
        if os.path.exists(file) != True:
            func(model1)
    output1 = out(model1)
    print('Assigned model directory:' + model1.outdir)
    return output1

def load_model(outdir='/m1_test/',read_previous=False):
    """ loads in a model from a directory location 
    by passing the corresponding pars.pkl file into parameters of a model class object
    
    Parameters:
    ----------
    outdir: str, directory where model will be initialized
    
    read_previous: boolean, default False, if True will read in data from existing .inp and .dat files
    
    Returns: output class object, see outputs.py
    
    """
    os.chdir(parent_dir+outdir)
    if os.path.exists('pars.pkl'):
        params = pickle.load( open( 'pars.pkl', "rb" ))
        model1 = model_(params,outdir=outdir)
        print('Loading from model directory:' + model1.outdir)
        
        file_list = ['amr_grid.inp', 'wavelength_micron.inp']
        func_list = [write_grid, write_wavelength]
        for file, func in zip(file_list,func_list):
            if os.path.exists(file) != True:
                func(model1)
        output1 = out(model1)
        
        if read_previous == True:
            prev_files = ['dust_density.inp','dust_temperature.dat','gas_temperature.inp']
            functions = [output1.read_rho, output1.read_Tdust, output1._fieldReader]
            for file, func in zip(prev_files, functions):
                func(file)
                
        return output1
    else:
        print('Could not locate parameters (pars.pkl) in',parent_dir+outdir)
        return 0
    


    
    