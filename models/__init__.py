from matplotlib import rc
from matplotlib import cm
from models.make_model import *
from models.radiation import *
from models.units import *
from models.outputs import *
from models.make_plots import *
import numpy as np


from pathlib import Path

p = Path().resolve()
parent_dir = str(p)
models_dir = parent_dir+'/models/'

np.seterr(divide='ignore')

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=False)
rc('mathtext', fontset = 'stix')
rc('axes', linewidth = 1.25)

SMALL_SIZE = 9
MEDIUM_SIZE = 10
BIGGER_SIZE = 11
BIGGEST_SIZE = 12

rc('font', size=SMALL_SIZE)          # controls default text sizes
rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=SMALL_SIZE, direction='in')    # fontsize of the tick labels
rc('ytick', labelsize=SMALL_SIZE, direction='in')    # fontsize of the tick labels
rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
rc('figure', titlesize=BIGGEST_SIZE) 

colors = cm.get_cmap('plasma',10)
plkw={'lw':'2','color':'black'} #plot a thick black line
plkw2 = {'lw':1,'color':'gray','ls':'dashed'}  #plot a thin gray dashed line
plkw3 = {'lw':1.5} # plot a thicker line
plsty = {'base':plkw,'ann':plkw2,'line':plkw3} #basic line, basic annotation, basic line, no color specified.