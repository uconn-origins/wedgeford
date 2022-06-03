#!/usr/bin/env python

import os
import matplotlib
import matplotlib.pyplot as plt
import optparse, time
import numpy as N
import pdb as pdb
from scipy import *
from matplotlib.mlab import griddata
from matplotlib import rc
from matplotlib import rcParams
from pylab import *
from subprocess import call
import fnmatch
import shutil

##################################################################################

def zetaeval(fineN,model):
    labels = ['m02','w98','ssm','ssx','ttm','ttx']
#Rates from Cleeves et al. 2014    
    allrts = N.array([[6.8e-16,3.7e-18,0.423,210.0],
        [2.0e-17,9.4e-19,0.021,260.0],
        [1.1e-18,3.0e-19,-0.00,260.0],
        [1.6e-19,7.5e-20,-0.01,250.0],
        [7.0e-21,4e-21,-0.01,290.0],
        [1.1e-22,3e-23,-0.02,490.0]])


# Updated cosmic ray rates as of April 2013:
#     allrts = N.array([[6.8e-16,3.7e-18,0.423,210.0],
#         [2.0e-17,9.4e-19,0.021,260.0],
#         [1.3e-18,3e-18,0.00,190.0],
#         [2e-19,8e-19,-0.01,230.0],
#         [1.0e-20,2e-19,-0.03,270.0],
#         [3.0e-22,2e-19,-0.03,270.0]])

    

    modin = labels.index(model)   
    zp = allrts[modin,0]
    ze = allrts[modin,1]
    al = allrts[modin,2]
    co = allrts[modin,3]
    mumult = 2.36
    piv = 1e20
    F = zp*ze*mumult/(ze*mumult*(fineN/piv)**al+zp*(exp(fineN*2.0*1.6733e-24/co)-1))
    return F

##################################################################################
    
def mainCR(bd):
     mp = 1.67e-24 
     cmAU = 1.496e13
     mu = 2.36 # Mean Mol. Wt.
     rn = 7.6e-19 # Umebayashi et al. 2008 SL + LL rates.    
     
     if bd[-1] == '/': bd = bd[:-1]  # trim off trailing slash if necessary.

     labels = ['m02','w98','ssm','ssx','ttm','ttx']
     baseenvdir = bd + '/'
     
     for e in range(6):
          tag = labels[e]

          destenvdir = bd + '_'+tag+'/'
                   
          pthfg = os.path.exists(destenvdir)
          if pthfg == False: 
                os.makedirs(destenvdir)   
          dust = 'e1' 
     
          dataformat = '%12.5E%14.5E%14.5E%14.5E%14.5E%14.5E%14.5E%14.5E'
          dataformat = '%12s%14.5E%14s%14s%14s%14s%14s%14s'
          dataformat = '%12s%14.5E%14s%14s%14s%14s%14.5E%14.5E'
          dataformat2 = '%12s%14.5E%14s%14s%14s%14s%14.5E%14.5E%14s'
          
          cmd = 'rm -f listI.txt; cd '+baseenvdir+'; ls -d 1env* | awk \'{newFile=$0; sub(/1environ.inp.'+dust+'./,"",newFile); printf("%s\\n",newFile);}\' > listI.txt'
          call(cmd,shell=True)
          
          f = open(baseenvdir+'listI.txt')
          radstr = f.readlines()
          f.close()
          
          fin = open(baseenvdir +'/1environ.inp.e1.'+radstr[0].rstrip(), 'r')
          nZ = np.size(fin.readlines())-3
          fin.close()
          
          nR = len(radstr)
          radval = N.zeros(nR)
          for i in range(nR):
               radstr[i] = radstr[i].rstrip()
               radval[i] = N.float(radstr[i])       
               arrtmp = genfromtxt(baseenvdir+'/1environ.inp.e1.'+radstr[i],skip_header=3)
          
               ncolv_down = N.zeros(nZ)
               ncolv_up = N.zeros(nZ)
               Ncol = 1e10 # starting column density
               z_last = 0.0
          
               for k in range(nZ):
                    rhoorig = arrtmp[k,1]
                    # fx = N.interp(radval[i],fixes[0:-2,0],fixes[0:-2,1])  
                    fx = 1.0
                    rho = rhoorig*fx#*0.5
                    z_here = N.float(arrtmp[k,5])
                    dx = (z_here-z_last) # in cmAU already.
                    Ncol = rho/(mu*mp)*dx + Ncol
                    ncolv_down[k] = Ncol           
                    z_last = z_here
               
               z_last = 0.0
          
               Ncol = ncolv_down[nZ-1]*2.0 # Go through both sides of the disk
               negNcol = 0.0
               for k in range(nZ):
                    ncolv_up[nZ-1-k] = 2.0*ncolv_down[nZ-1] - ncolv_down[nZ-1-k]

               zeta = 0.5*zetaeval(ncolv_down,tag) + 0.5*zetaeval(ncolv_up,tag)   # 'm02'
               #pdb.set_trace()       
               fin = open(baseenvdir+'/1environ.inp.e1.'+radstr[i], 'r')
               fout = open(destenvdir+'/1environ.inp.e1.'+radstr[i], 'w')
               ct = 0
               v = 0

               for line in fin:
               
                    if ct > 2:  # Skip header, modify values below.  Here modifying rho_gas.
                         temp = line.rstrip()
                         info = temp.split('  ')
                         rhoorig = N.float(info[1])
                         Nrz = N.float(info[6])
                         if Nrz < 0:
                              Nrz = 0.0
                         #fx = N.interp(radval[i],fixes[0:-2,0],fixes[0:-2,1])  # Interpolate scaling factor.
                         fx = 1.0
                         rhonew = rhoorig*fx#*0.5
                         if len(info)==9:
                         	formatted = dataformat2 % (info[0],rhonew,info[2],info[3],info[4],info[5],Nrz,zeta[v],info[8])
                         else:
                         	formatted = dataformat % (info[0],rhonew,info[2],info[3],info[4],info[5],Nrz,zeta[v])
                         fout.write(formatted+'\n')
                         v = v+1
                    else:
                         fout.write(line)
                    ct = ct + 1     
               fin.close()
               fout.close()

     return

##################################################################################

def mainTG(baseenvdir,destenvdir,Tgasfile):
     mp = 1.67e-24 
     cmAU = 1.496e13
     mu = 2.36 # Mean Mol. Wt.
     
     constcrrt = [1e-16,1e-17,1e-19,1e-20]
     for e in range(1):
          t = os.getcwd()
          Tgasmod = genfromtxt(baseenvdir+Tgasfile,skip_header=0)
          
          pthfg = os.path.exists(destenvdir)
          if pthfg == False: 
                os.makedirs(destenvdir)   
          dust = 'e1' 
  
          dataformat = '%12s%14.5E%14.5E%14s%14s%14s%14s%14.5E'
          dataformat2 = '%12s%14.5E%14.5E%14s%14s%14s%14s%14.5E%14s'
          
          cmd = 'rm -f listI.txt; cd '+baseenvdir+'; ls -d 1env* | awk \'{newFile=$0; sub(/1environ.inp.'+dust+'./,"",newFile); printf("%s\\n",newFile);}\' > listI.txt'
          call(cmd,shell=True)
          
          f = open(baseenvdir+'listI.txt')
          radstr = f.readlines()
          f.close()
          
          fin = open(baseenvdir +'/1environ.inp.e1.'+radstr[0].rstrip(), 'r')
          nZ = np.size(fin.readlines())-3
          fin.close()
          
          nR = len(radstr)
          radval = N.zeros(nR)
          for i in range(nR):
               #pdb.set_trace()
               radstr[i] = radstr[i].rstrip()
               radval[i] = N.float(radstr[i])  

               indv = N.where(N.abs(Tgasmod[:,0]-radval[i])/radval[i] < 1e-5)

               Tgas_rho = Tgasmod[indv,2]               
               TgasZ = Tgasmod[indv,1]
               TgasT = Tgasmod[indv,-1]
               TgasT =TgasT[0]
               TgasZ = TgasZ[0]
    
               fin = open(baseenvdir+'/1environ.inp.e1.'+radstr[i], 'r')
               fout = open(destenvdir+'/1environ.inp.e1.'+radstr[i], 'w')
               ct = 0
               v = 0

               for line in fin:
               
                    if ct > 2:  # Skip header, modify values below.  Here modifying rho_gas.
                         temp = line.rstrip()
                         info = temp.split('  ')
                         
                         rhoorig = N.float(info[1])
                         zval = N.float(info[4])
                         indz = N.argmin(N.abs(TgasZ-zval))
                         crval = N.float(info[7])
                         #crval = constcrrt[e]
                         #crval = 0.0
                         
                         NewTGas = TgasT[indz]
                         if NewTGas > 4000.0:
                              rhoorig= 0.0
                              

                         if NewTGas == 0:  
                              rhoorig = 0.0
                         if rhoorig < 1e-24:
                              rhoorig = 0.0
                         #rhoorig = rhoorig*0.55
                         if len(info)==9:                   
                         	formatted = dataformat2 % (info[0],rhoorig*0.75,TgasT[indz],info[3],info[4],info[5],info[6],crval,info[8])
                         else:
                         	formatted = dataformat % (info[0],rhoorig*0.75,TgasT[indz],info[3],info[4],info[5],info[6],crval)
                         fout.write(formatted+'\n')
                         v = v+1
                         #pdb.set_trace()
                    else:
                         fout.write(line)
                    ct = ct + 1     
               fin.close()
               fout.close()
     return

##################################################################################

def Generate_Environment(maindir,runnm):

    basedir = maindir+'/environ/'+runnm+'/'
    destdir = maindir+'/environ/'+runnm+'_Tgas'+'/'
    tgfile  = runnm+'_model_Tgas_SB.txt'
    
    # Incorporate gas temperatures into the env files:        
    mainTG(basedir,destdir,tgfile)
    
    # Incorporate CR ionization rates.
    mainCR(destdir)

    return
    
##################################################################################


def Tgas_Comp(maindir,runnm,uvfile,shortrunnm):    
    c = 2.99792458e10
    mh = 1.67e-24
    mu = 2.36

    # Read in model information:
    mod = genfromtxt(maindir+'/environ/'+runnm+'/gas_'+shortrunnm+'.out',skip_header=2)
    #mod = genfromtxt(maindir+'/environ/'+runnm+'/gas_gas_parameters_'+runnm+'.out',skip_header=2)
    
    # Output file names:
    outname1      =  maindir+'/environ/'+runnm+'/'+runnm+"_Tgas_SB.txt"
    outname2      =  maindir+'/environ/'+runnm+'/'+runnm+"_model_Tgas_SB.txt"
    finp = maindir + '/' + uvfile

    h = 6.6260755e-27
    c = 2.99792458e10
    
    nel = N.size(mod[:,0])
    lvec = N.linspace(1,108,108)*10+920
    evec = h*c/(lvec*1e-8)
    allflux = N.zeros((nel,3))

    fi = open(finp,'r')
    
    nlin = len(fi.readlines())
    fi.seek(0) 
    junk = fi.readline()
    numblocks = (nlin-1)/111
    
    ct = 0
    
    Habing = 1.83590e+08  # Photons/cm^2/s between 930-2000
    G0 = 2.7e-3 # erg/cm^2/s^-1

    for i in range(numblocks):
        tmp = fi.readline()
        tmp = tmp.split()
        radval = N.float(tmp[1])
        tmp = fi.readline()
        tmp = tmp.split()
        zstr = tmp[1:]
        zvals = N.array([N.float(x) for x in zstr]) # Convert to numbers
        junk = fi.readline()  
        
        photstore = N.zeros((108,len(zvals)))
        lam = N.linspace(1,108,108)
        
        # Read in the blocks at each radius with UV photon fluxes.
        for w in range(108):
            tmp = fi.readline()
            tmp = tmp.split()
            phots = N.array([N.float(x) for x in tmp])
            photstore[w,:] = phots[1:] # at a given wavelength.  Flux times ionization cross section.
 
        # photstore[28,:]=(photstore[27,:]+photstore[29,:])/2.0  # interpolate lya
                
        for k in range(len(zvals)):    
            allflux[ct,0] = radval
            allflux[ct,1] = zvals[k]
            allflux[ct,2] = N.trapz(photstore[:,k]*evec,x=lvec,axis=0)/G0 # in Habings
            ct = ct + 1   
    fi.close()
    
    NH2 = mod[:,2]/(mu*mh)
    NH = 2.0*NH2
    x = N.log10(NH/1e5)
    x[x>9.2920]=9.29199  # Max density where this works is ~1e14 cm^-3.

    delT = 420.2*(allflux[:,2]/1e3)**(1.05-0.113*x)*1.0/(10**(0.486*x-0.014*x**2))
    
    #delT = 420.2*(allflux[:,2])**(1.05-0.113*x)*1.0/(10**(0.486*x-0.014*x**2))

    # Masking lyman alpha:
    #delT = 420.2*(allflux[3720:3735,2])**(1.05-0.113*x[3720:3735])*1.0/(10**(0.486*x[3720:3735]-0.014*x[3720:3735]**2))
    
    
    delT[allflux[:,2]<=0.0]=0.0
    delT[NH<1e-2] = 0.0
    Tgas = mod[:,4]+delT
    i = N.where(Tgas >4200.0)
    Tgas[i] = 4200.0
    
    i = N.where(mod[:,2]<1e-30)
    mod[i,2] = 0.0
    mod[i,3] = 0.0
     
    gaparray = allflux*1.0
    gaparray[:,2] = Tgas
    
    savef = 1
    if savef == 1: 
        N.savetxt(outname1, gaparray)
        tg = N.reshape(Tgas,(N.size(Tgas),1))
        newarray = N.concatenate((mod,tg),axis=1)
        fmtstr = '%11.6f\t%11.6f\t%13.5e\t%13.5e\t%11.5f\t%11.5f\t%11.5f\t%11.5f\t%11.5f\t%11.5f\t%11.5f'
        N.savetxt(outname2, newarray, fmt=fmtstr, delimiter='\t')

    return
    
##################################################################################
def Create_0io(dir,runnm,uv,xr,isrf,slr,st,rre,rsp,ioname):

    finp = dir+'/environ/'+runnm+'/'+ioname
    fi = open(finp,'w')
    fi.write('# input & output files:\n')
    fi.write('{:<45}'.format(rsp)+'# file with species\n')
    fi.write('{:<45}'.format(rre)+'# file with reactions\n')
    fi.write('{:<45}'.format(uv)+'# file with uvfield\n')
    fi.write('{:<45}'.format(xr)+'# file with xrfield\n')
    fi.write('{:<45}'.format(isrf)+'# file with ISRF\n')
    fi.write('{:<45}'.format(slr)+'# Radionuclide Ion Rate\n')
    fi.write('{:<45}'.format(st)+'# 2D abundance file ending\n')
    fi.close()
    
    for file in os.listdir(dir+'/environ/'):
        if fnmatch.fnmatch(file,runnm+'*'):
            foundfile = file
            if file != runnm:
                shutil.copy2(finp,dir+'/environ/'+foundfile+'/'+ioname)

    return    

##################################################################################
# Incorporate additional gas temperature and CR information into environment files.
# This code will expect that the UV file is in the main disk chemistry directory and 
# that the environment file contains a copy of the relevant gas_Model.out file. 
#
# The script can also produce 0io files and copy into all relevant directories.

# Location of disk chemical model:
diskchem_dir = '/Nirgal1/kamberrs/disk_chemistry/MasterChemistry'

#cp -r /Nirgal1/kamberrs/torus/environ/Ms2.2Md0.0* .
#10micron_0.03Mslg_0.99m_d60.0/
# Name(s) of torus environment output:
mass = ['0.1Ms']
extraflag =['']
#lg = ['lg_0.0m','lg_0.9m','lg_0.99m','lg_0.9m','lg_0.9m','lg_0.9m','lg_0.99m','lg_0.99m','lg_0.99m']
#Routmm = [200.,200.,200.,158.,123.,60.,158.,123.,60.]
lg = ['lg_0.8m','lg_0.1m','lg_0.2m','lg_0.3m','lg_0.4m','lg_0.5m','lg_0.6m','lg_0.7m','lg_0.8m','lg_0.9m','lg_0.99m']
Routmm = [200.,200.,200.,200.,200.,200.,200.,200.,200.,200.,200.]

# 0.003Mslg_0.0m_d200.0
drift = False
#lg = ['lg_0.9m','lg_0.675m','lg_0.45m','lg_0.225m','lg_0.0m']
#lg = ['lg_0.45m','lg_0.225m','lg_0.0m']
#srt_time_vals = ['5.00E+05','1.50E+06','6.00E+06','1.32E+06','2.15E+06','2.97E+06','6.00E+06','5.001D+05','1.32E+06','1.50E+06','2.146D+06','2.97E+06','6.00E+06'
srt_time_vals = ['None','5.0E+05','1.3E+06','1.5E+06','2.1E+06','3.0E+06']#,'None'] #round
srt_time_vals = ['None','None','None','None','None','None','None','None','None','None','None']
env_orig = []#np.empty(len(mass)*len(lg),dtype=string)
short_env_orig = []
uv_model = []#np.empty_like(env_orig)
slr_model = []
srt_time = []


#env_orig = ['0.1Mslg_0.9m','0.1Mslg_0.675m','0.1Mslg_0.45m','0.1Mslg_0.225m','0.1Mslg_0.0m']
xr_vals = ['29.0','30.0','31.0']
isrf_vals = ['1.0']#,'300.0','3000.0']
rreac = 'rreacs_herb0308_CDR_FullGr.dat'
rspec = 'rspecies_herb0308_CDR_FullGr.dat'

# UV photon files:
for m in range(len(mass)):
	for l in range(len(lg)):
		env_orig.append(extraflag[m]+mass[m]+lg[l]+'_d'+str(Routmm[l]))
		short_env_orig.append(mass[m]+lg[l]+'_d'+str(Routmm[l]))
		uv_model.append('uv_photons_'+extraflag[m]+mass[m]+lg[l]+'_d'+str(Routmm[l])+'_uv_e1.dat')
		slr_model.append('None')
		srt_time.append(srt_time_vals[l])

# Additional files (if generating 0io files):
#xr_model = ['xray_phot_0.1Mslg_0.9m_xray_Lx29.0.dat','xray_phot_0.1Mslg_0.675m_xray_Lx29.0.dat','xray_phots_0.1Mslg_0.45m_xray_Lx29.0.dat','xray_phot_0.1Mslg_0.225m_xray_Lx29.0.dat','xray_phot_0.1Mslg_0.0m_xray_Lx29.0.dat']
#isrf_model = ['isrf_1_half_0.1Mslg_0.9m.dat','isrf_1_half_0.1Mslg_0.675m.dat','isrf_1_half_0.1Mslg_0.45m.dat','isrf_1_half_0.1Mslg_0.225m.dat','isrf_1_half_0.1Mslg_0.0m.dat']
#ioname = '0io.1G0.x29.inp'

for x in range(len(xr_vals)):
	for e in range(len(env_orig)):
		xr_model = 'xray_phot_'+short_env_orig[e]+'_xray_Lx'+xr_vals[x]+'.dat'
		for i in range(len(isrf_vals)):
			isrf_model = 'isrf_'+isrf_vals[i]+'_half_'+env_orig[e]+'.dat'
			isrf_model = 'None'
			ioname = '0io.'+str(int(float(isrf_vals[i])))+'G0.x'+str(int(float(xr_vals[x])))+'.inp'
			Tgas_Comp(diskchem_dir,env_orig[e],uv_model[e],short_env_orig[e])
			Generate_Environment(diskchem_dir,env_orig[e])
			Create_0io(diskchem_dir,env_orig[e],uv_model[e],xr_model,isrf_model,slr_model[e],srt_time[e],rreac,rspec,ioname)
# for j in range(len(env_orig)):
#     Tgas_Comp(diskchem_dir,env_orig[j],uv_model[j])
#     Generate_Environment(diskchem_dir,env_orig[j])
#     Create_0io(diskchem_dir,env_orig[j],uv_model[j],xr_model[j],isrf_model[j],slr_model[j],srt_time[j],rreac,rspec,ioname)
#     print 'Compiled environment files for model %s successfully.'%(env_orig[j])

