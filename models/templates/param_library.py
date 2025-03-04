from pylab import *
### library of pre-set parameters, can put your own in here if you plan to use it a lot.
library = {}

library['default'] = {'Ms': 1, 'Rs': 2.0, 'Ts': 5000, 'accrate':1e-7,'f':0.01, 'Mdisk': 0.1, 'Mfrac': [0.01,0.02],'R0':[5,5], 'Rout':[100,50], 'H0':[1,0.3], 'p':[-1,-1], 'Rdisk':[125,125],'Min': 1e-6, 'Rc':125, 'rho_amb':1e-25, 'rho_0': 3e-22,'theta_min': 25,'exf':0.25,'Rmax':1.5e4, 'd2g': 0.01, 'shock':False,'N':[180,90,48], 'min':[0.1,pi/16.,0], 'max':[400,pi/2.,2*pi], 'spacing':['log','lin','lin'],'rho_si':3.1518, 'amax_ism': 1.0, 'amin': [0.005,0.005], 'amax': [1,1e3], 'apow': [3.5,3.5],'cr_model': 'ssx','zetacr': 1.3e-17, 'LX': 1e30, 'G0':1, 'viscous_heating':False}

library['l1527_a'] = {'Ms': 0.2, 'Rs': 2.0, 'Ts': 4000.0, 'accrate': 1e-07, 'f': 0.01, 'Mdisk': 0.010499999999999999, 'Mfrac': [0.005, 0.005], 'R0': [1, 1], 'H0': [1, 0.3], 'p': [-1, -1], 'Rdisk': [105, 105], 'Min': 3e-06, 'Rc': 105, 'rho_amb': 1e-25, 'rho_0': 3e-22, 'theta_min': 15, 'exf': 0.25, 'Rmax': 15000.0, 'd2g': 0.01, 'shock': True, 'N': [120, 90, 48], 'min': [0.1, 0.19634954084936207, 0], 
'max': [400, 1.5707963267948966, 6.283185307179586], 'spacing': ['log', 'lin', 'lin'], 'rho_si': 3.1518, 'amin_chem': 0.06, 'amax_ism': 1.0,'amin': [0.005, 0.005], 'amax': [1, 1000.0], 'apow': [3.5,3.5], 'cr_model': 'ssx', 'zetacr': 1.3e-17, 'LX': 1e+30, 'G0': 1, 'viscous_heating': True}

library['classII_disk'] = {'Ms': 0.5, 'Rs': 2.0, 'Ts': 4000.0, 'accrate': 1e-08, 'f': 0.01, 'Mdisk': 0.1, 'Mfrac': [0.005, 0.005], 'R0': [1, 1], 'H0': [1, 0.2], 'p': [-1, -1], 'Rdisk': [60, 60], 'Min': 0, 'Rc': 100, 'rho_amb': 1e-25, 'shock': False, 'N': [160, 90, 48], 'min': [0.1, 0.19634954084936207, 0], 'max': [400, 1.5707963267948966, 6.283185307179586], 'spacing': ['log', 'lin', 'lin'], 'rho_si': 3.1518, 'amin_chem': 0.06, 'amax_ism': 1.0,'amin': [0.005, 0.005], 'amax': [1, 1000.0], 'apow': [3.5,3.5], 'cr_model': 'ssx', 'zetacr': 1.3e-17, 'LX': 1e+30, 'G0': 1, 'viscous_heating': True}


library['gma'] = {'Ms': 1.1, 'Rs': 1.9, 'Ts': 4350, 'accrate':8e-9,'f':0.01, 'Mdisk': 0.2, 'Mfrac': [0.0005,0.003],'R0':[1,1], 'H0':[1,0.2], 'p':[-1,-1], 'Rdisk':[175,175],'Min': 1e-6, 'Rc':90, 'rho_amb':1e-36, 'rho_0': 3e-22,'theta_min': 25,'exf':0.25,'Rmax':1.5e4, 'd2g': 0.01, 'shock':False,'N':[300,180,90], 'min':[1,0.7854,0], 'max':[650,pi/2.,2*pi], 'spacing':['log','lin','lin'],'rho_si':1.675, 'amin_chem':0.06, 'amax_ism': 1.0, 'amin': [0.005,0.005], 'amax': [1,1e3], 'apow': [3.5,3.5],'cr_model': 'ssx','zetacr': 1.3e-17, 'LX': 1e30, 'G0':0, 'viscous_heating':False, 'stream_frac':0.75, 'nstreams':6}


library['fiducial'] = {'Ms': .2, 'Rs': 3.0, 'Ts': 3485, 'accrate':8.5e-9,'f':0.01, 'Mdisk': 0.2, 'Mfrac': [0.0005,0.003],'R0':[1,1], 'H0':[1,0.2], 'p':[-1,-1], 'Rdisk':[175,175],'Min': 1e-6, 'Rc':90, 'rho_amb':1e-36, 'rho_0': 3e-22,'theta_min': 25,'exf':0.25,'Rmax':1.5e4, 'd2g': 0.01, 'shock':False,'N':[300,180,90], 'min':[1,0.7854,0], 'max':[650,pi/2.,2*pi], 'spacing':['log','lin','lin'],'rho_si':1.675, 'amin_chem':0.06, 'amax_ism': 1.0, 'amin': [0.005,0.005], 'amax': [1,1e3], 'apow': [3.5,3.5],'cr_model': 'ssx','zetacr': 1.3e-17, 'LX': 1e30, 'G0':0, 'viscous_heating':False, 'stream_frac':0.75, 'nstreams':6}

library['classI'] = {'Ms': 0.8, 'Rs': 2.0, 'Ts': 4000, 'accrate':1.0e-7,'f':0.01, 'Mdisk': 0.2, 'Mfrac': [0.0005,0.003],'R0':[0.1,0.1], 'H0':[1,0.2], 'p':[-1,-1], 'Rdisk':[60,60],'Min': 1e-6, 'Rc':90, 'rho_amb':1e-26, 'rho_0': 3e-22,'theta_min': 25,'exf':0.25,'Rmax':1.5e4, 'd2g': 0.01, 'shock':False,'N':[300,180,90], 'min':[0.05,pi/16.,0], 'max':[650,pi/2.,2*pi], 'spacing':['log','lin','lin'],'rho_si':1.675, 'amin_chem':0.06, 'amax_ism': 1.0, 'amin': [0.005,0.005], 'amax': [1,1e3], 'apow': [3.5,3.5],'cr_model': 'ssx','zetacr': 1.3e-17, 'LX': 1e30, 'G0':0, 'viscous_heating':False, 'stream_frac':1, 'nstreams':1}

library['TW_Hya'] = {'Ms': 0.8, 'Rs': 1.11, 'Ts': 4000.0, 'accrate': 2.5e-9, 'f': 0.5, 
                    'Mdisk': 0.023, 'Mfrac': [0.0005, 0.005], 'R0': [0.1, 0.1], 'H0': [1, 0.2], 'p': [-1, -1], 'Rdisk': [35, 35], 
                    'Min': 0, 'Rc': 100, 'rho_amb': 1e-25, 'shock': False, 'N': [160, 90, 48],
                    'min': [0.1, pi/4, 0], 'max': [200, 1.5707963267948966, 6.283185307179586],
                    'spacing': ['log', 'lin', 'lin'], 'rho_si': 3.1518, 'amin_chem': 0.06, 'amax_ism': 1.0,
                    'amin': [0.005, 0.005], 'amax': [1, 1000.0], 'apow': [3.5,3.5], 
                    'cr_model': 'ssx', 'zetacr': 1.3e-17, 'LX': 1e+30, 'G0': 1, 'viscous_heating': False}

library['HL_Tau'] = {'Ms': 1.7, 'Rs': 7.0, 'Ts': 4000.0, 'accrate': 1.67e-7, 'f': 0.01, 
                    'Mdisk': 0.2, 'Mfrac': [0.002, 0.008], 'R0': [0.1, 0.1], 'H0': [1, 0.2], 'p': [-1, -1], 'Rdisk': [60, 60], 
                    'Min': 1e-6, 'Rc': 65, 'rho_amb': 1e-25, 'theta_min': 20, 'shock': False, 'N': [160, 90, 48],
                    'min': [0.1, 0.25, 0], 'max': [200, 1.5707963267948966, 6.283185307179586],
                    'spacing': ['log', 'lin', 'lin'], 'rho_si': 3.1518, 'amin_chem': 0.06, 'amax_ism': 1.0,
                    'amin': [0.005, 0.005], 'amax': [1, 1000.0], 'apow': [3.5,3.5], 
                    'cr_model': 'ssx', 'zetacr': 1.3e-17, 'LX': 1e+30, 'G0': 1, 'viscous_heating': False}
