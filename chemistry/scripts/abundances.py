import os
import matplotlib
import matplotlib.pyplot as plt
import optparse, time
import numpy as np
from scipy import *
from matplotlib.mlab import griddata
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
import pdb as pdb
import sys

cmAU=6.68458712e-14 #convert from cm to AU
mAU=6.68458712e-12 #convert from meters to AU

def round_to_n(x,n): #stole this from stack overflow
#original .out file only goes to 3 decimal places
    return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)

def readfile(filename,radius):
    infile=open(filename,"r")
    z=[]
    abun=[]
    for line in infile:
        linesplit=line.split()
        if eval(linesplit[0]) == radius:
            abun.append(eval(linesplit[5])*eval(linesplit[2]))
            z.append(eval(linesplit[1]))
            

    infile.close()

    return abun,z

def round_to_n(x,n): #stole this from stack overflow

    return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)

def main():
    timesteps = 50
    model = ['test21']
    path = '/Users/kamberrs/diskchemmay25/runs/'
    spec = ['CO','YO','CZ'] 

    radius = 7.5323600E+12 #50 AU
    #radius = 1.5199360E+11 #1 AU
    #radius = 5.2285200E+13 #349 AU
    #radius = 3.0668000E+13 #204 AU
    time = 45

    abun=np.zeros((len(spec)*len(model),45),dtype=np.float64)
    z=np.zeros((len(spec)*len(model),45),dtype=np.float64)


    for j in range(len(model)):
        for i in range(len(spec)):
            filename=path+model[j]+'/e1/lime_'+spec[i]+'/'+spec[i]+'_time'+str(time)+'.dat'
            abun[i+(j*len(spec)),:],z[i+(j*len(spec)),:]=readfile(filename,radius)
    

    z=z*mAU
  
 
#    sys.exit("Mathmatical!")

    ratio=np.zeros(15,dtype=np.float64)
    ratio=60*abun[1,:]/abun[0,:]

    err=0.0005
    x=[0,20]
    shade={'lw':0.0,'edgecolor':None}#kwargs for fill_between

    #print 60*abun[1,40]/abun[0,40], z[1,40]

    #x = [13.49,12.32,11.16,10.015,8.88,7.75,6.62]#UV
    #y = [2.27E+05,9.35E+04,4.15E+04,1.55E+04,3.53E+03,5.92E+02,6.10E+01]#UV
    #y = [3.67E+06,4.72E+05,9.51E+03,2.93E+02,1.09E+02,3.45E+00,1.72E+00]#xray
    x = [13.491,12.321,11.162,10.015,8.878,7.750,6.629]#xray
    y = [8.57200e-16,1.52200e-15,2.57000e-15,4.12800e-15,6.30900e-15,9.17200e-15,1.26900e-14]#gas density
    y = [1.28600e-18,2.28200e-18,3.85400e-18,6.19200e-18,9.46300e-18,1.37600e-17,1.90300e-17]#dust density
    y = [53.00000,45.10000,40.70000,38.2,35.5,32.9,31.1]#dust temp

    thick = 2
    plt.plot(x,y,lw=thick)
    #plt.semilogy(z[0,:],abun[0,:],lw=thick,label=spec[0],color='purple',marker='o')
    #plt.semilogy(z[1,:],60*abun[1,:],lw=thick,label=spec[1],color='orange',marker='o')
    #plt.semilogy(z[2,:],500*abun[2,:],lw=thick,label=spec[2],color='cyan',marker='o')
    #plt.semilogx(time,500*abun[3,:],lw=thick,label=spec[3],color='#3385D6')
    #plt.legend((spec[0],spec[1],spec[2]),loc=0,prop=dict(size=15))
    #plt.fill_between(x,1-err,1+err,facecolor='0.7',**shade)
    #plt.xlim(0,20)
    #plt.ylim(1.01,0.99)
    plt.xlabel('Z (AU)',size=20)
    #plt.ylabel('n(X) (m$^{-3}$)',size=20)
    plt.ylabel('T dust (K)')
    plt.show()
          
main()
