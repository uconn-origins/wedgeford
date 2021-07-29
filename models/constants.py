# units and conversions
Msun = 1.9891e33 # g
Rsun = 69.6e9 #cm
Lsun  = 3.8525e33 #erg/s
AU = 1.49598e13 #cm
yr = 3.14e7 #seconds
cmtokm = 10**(-5) #converts cm to km
R_mu = 36149835 # (R/mu) for 2.4 g/mol in PPD
Gconv = 6.6743e-8 # cgs
S0conv = (Msun/(AU**2)) #S0cgs = S0code*S0conv
sigsb = 5.67e-5 #Stefan Boltzmann constant
c = 3e10 # speed of light

from matplotlib import rc
from matplotlib import cm

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=False)
rc('mathtext', fontset = 'stix')

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14
BIGGEST_SIZE = 16

rc('font', size=SMALL_SIZE)          # controls default text sizes
rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
rc('axes', labelsize=BIGGEST_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
rc('figure', titlesize=BIGGEST_SIZE) 


colors = cm.get_cmap('plasma',10)
plkw={'lw':'2','color':'black'} #plot a thick black line
plkw2 = {'lw':1,'color':'gray','ls':'dashed'}  #plot a thin gray dashed line
plkw3 = {'lw':1.5} # plot a thicker line
plsty = {'base':plkw,'ann':plkw2,'line':plkw3} #basic line, basic annotation, basic line, no color specified.