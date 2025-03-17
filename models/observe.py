from pylab import *
#import radmc3dPy as rpy
import numpy as np
import os
from models.units import *
from models.outputs import *

class observe:
    def __init__(self,output,incl=0,wav=870,dpc=140,npix=300,sizeau=300,beam_au=[30,20],fwhm=[None,None]):
        pars = {'incl':incl,'wav':wav,'npix':npix,'sizeau':sizeau}
        self.dpc = dpc
        if fwhm[0] != None:
            self.fwhm = fwhm
        else:
            self.fwhm = [206265*(AU/(dpc*pc)) * i for i in beam_au]
        self.I = {}
        self.Ic = {}
        #self.i_do = self.make_continuum_image(output,envelope=False,**pars)
        #self.i_bo = self.make_continuum_image(output,envelope=True,**pars)
        #self.I = {'disk': self.i_do.image, 'both':self.i_bo.image, 'env': self.i_bo.image - self.i_do.image}
        #self.Ic = {'disk': self.i_do.imConv(fwhm=self.fwhm, pa=0., dpc=self.dpc).image, 'both':self.i_bo.imConv(fwhm=self.fwhm, pa=0., dpc=self.dpc).image}
        #self.Ic['env'] = self.Ic['both'] - self.Ic['disk']
        self.make_continuum_image(output,str(pars['wav']),**pars)
        
    def make_continuum_image(self,output,key,**params):
        model=output.m
        write_main(model,scat=1,lines=0)
        if os.getcwd() != model.outdir:
            os.chdir(model.outdir)
        rpy.image.makeImage(**params)
        im=rpy.image.readImage(binary=False)
        self.nfreq = im.nfreq
        self.nx = im.nx
        self.ny = im.ny
        self.y = im.y/AU
        self.x = im.x/AU
        self.I[key] = im
    
    def convolve(self,key,beam_au='',pa=0.):
        if len(beam_au) < 1:
            fwhm = self.fwhm
        else:
            fwhm = [206265*(AU/(self.dpc*pc)) * i for i in beam_au]
        self.Ic[key] = self.I[key].imConv(fwhm=fwhm, pa=pa, dpc=self.dpc).image
        
    def z_slice(self,I,z=0):
        zs = {}
        for key in I.keys():
            imarr = I[key]
            if imarr.ndim > 2:
                imarr = imarr[:,:,int(self.nfreq/2)-1]
            index = int(self.ny/2) -1 + int(0.5*self.ny*z/np.amax(self.y))
            zs[key] = imarr[index,:]
        return zs
    
    def r_slice(self,I,r=0):
        rs = {}
        for key in I.keys():
            imarr = I[key]
            if imarr.ndim > 2:
                imarr = imarr[:,:,int(self.nfreq/2)-1]
            index = int(self.nx/2) -1 + int(0.5*self.nx*r/np.amax(self.x))
            rs[key] = imarr[:,index]
        return rs
    
