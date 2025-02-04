from pylab import *
#import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata

from matplotlib import rc
from matplotlib import cm

from sys import platform

from . units import *
from . outputs import *

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=False)
#rc('mathtext', fontset = 'stixsans')
rc('axes', linewidth = 1.25)

SMALL_SIZE = 8
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
rc('pdf',fonttype=42)

colors = cm.get_cmap('plasma',10)
plkw={'lw':'2','color':'black'} #plot a thick black line
plkw2 = {'lw':1,'color':'gray','ls':'dashed'}  #plot a thin gray dashed line
plkw3 = {'lw':1.5} # plot a thicker line
plsty = {'base':plkw,'ann':plkw2,'line':plkw3} #basic line, basic annotation, basic line, no color specified.

##### plotting convenience functions #####
def plot_slice(output,rho=None,log=True,average = True, **plot_params):
    model = output.m
    R_CYL,Z_CYL = model.make_rz()
    if rho is not None:
        if rho.ndim > R_CYL[:,:,0].ndim:
            if average == True:
                rho2d = np.average(rho,axis=-1) #azimuthal average 
            else:
                rho2d = rho[:,:,0]
        else:
            rho2d = rho
        if np.shape(rho2d) != np.shape(R_CYL[:,:,0]):
            rho2d = rho2d.T
        else:
            rho2d = rho2d
        ax = gca()
        if log == True:
            im=ax.contourf(R_CYL[:,:,0],Z_CYL[:,:,0], np.log10(rho2d),**plot_params)
        else:
            im=ax.contourf(R_CYL[:,:,0],Z_CYL[:,:,0], rho2d,**plot_params)
        return im

def plot_contour(output,rho=None,thresh=None,log=False,average=True,**plot_params):
    model = output.m
    R_CYL,Z_CYL = model.make_rz()
    if rho is not None:
        if rho.ndim > R_CYL[:,:,0].ndim:
            if average == True:
                rho2d = np.average(rho,axis=-1)
            else:
                rho2d = rho[:,:,0]
        else:
            rho2d = rho
        if np.shape(rho) != np.shape(R_CYL[:,:,0]):
            rho2d = rho2d.T
        else:
            rho2d = rho2d
        ax = gca()
        if log == True:
            cs = ax.contour(R_CYL[:,:,0],Z_CYL[:,:,0], np.log10(rho2d),**plot_params,levels=thresh)
        else:
            cs = ax.contour(R_CYL[:,:,0],Z_CYL[:,:,0], rho2d,**plot_params,levels=thresh)
        ax.clabel(cs, cs.levels)
        return cs
     
        
def plot_components(output,fluid=0,rlim=400):
    model = output.m
    f,ax= subplots(1,3,constrained_layout=True,dpi=100)
    fluids = ['gas','small dust', 'large dust']
    f.suptitle(fluids[fluid],fontsize=11)
    f.set_size_inches(6,2)
    components = {'Disk': model.rho_disk(fluid=fluid),'Envelope': model.rho_embedded(fluid=fluid) - model.rho_disk(fluid=fluid), 'Disk + Envelope': model.rho_embedded(fluid=fluid)}
    if fluid == 2:
        components['Envelope'] = None
    for a,c in zip(ax, components.keys()):
        sca(a)
        im=plot_slice(output,rho=components[c],log=True,levels=np.arange(model.rhomin,model.rhomax,0.5))
        a.set_title(c,fontsize=11)
        a.set_xlim(model.grid['min'][0],model.grid['max'][0])
        a.yaxis.set_major_locator(MultipleLocator(int(rlim/4)))
        a.xaxis.set_major_locator(MultipleLocator(int(rlim/4)))
        if a != ax[0]:
            a.yaxis.set_major_formatter(NullFormatter())
        if c == 'Disk':
            for j in np.arange(1,4):
                a.plot(model.r, model.H(model.r)*j,color='white',lw=1,label=str(j)+'H') #plot disk scale heights
        elif c == "Envelope":
            for j in np.linspace(np.radians(model.env['theta_min']),pi/2.,8):
                streamline = model.stream(th0=j,shock=False)
                r = streamline['path'][0]
                th = streamline['path'][1]
                l = a.plot(r*np.sin(th),r*np.cos(th),color='white',ls='dashed',lw=1) #plot streamlines
        elif c == "Disk + Envelope" and model.env['shock'] == True:
            xpts = []
            ypts = []
            for j in np.linspace(np.radians(model.env['theta_min']),pi/2.,30):
                streamline = model.stream(th0=j,shock=True)
                r = streamline['path'][0]
                th = streamline['path'][1]
                Ts = streamline['Tg'] - model.T(r)
                s = np.argmax(Ts)
                xpts.append(r[s]*np.sin(th[s]))
                ypts.append(r[s]*np.cos(th[s]))
            a.plot(xpts,ypts,color='white',lw=1)
        a.set_ylim(0,rlim)
        a.set_xlim(0,rlim)
    colorbar(im,ax=ax,label=r'$\rho \ \mathrm{[g \ cm^{-3}]}$')
    ax[1].set_xlabel('R [au]')
    ax[0].set_ylabel('z [au]')



def plot_midplaneT(output):
    Tdust_1 = output.T['dust']
    Tdust2D_1 = output.calc_T2D(fluid='dust')
    Tdust_mid_r_1 = np.average(Tdust2D_1[:,-10:],axis=-1)
    f, ax = subplots(1)
    ax.loglog(output.r, Tdust_mid_r_1)


def plot_velocities(output,rlim=400):
    model = output.m
    f,ax= subplots(1,3,constrained_layout=True,dpi=100)
    f.set_size_inches(6,2)
    vectors = [r'$v_r$',r'$v_{\theta}$',r'$v_{\rm in}$']
    keys = ['r','theta','phi']
    f.set_size_inches(6,2)
    vfield = model.v_embedded()
    vkep = model.v_disk()['phi']
    components = dict(zip(vectors, [np.abs(vfield['r']),np.abs(vfield['theta']), np.sqrt(vfield['r']**2 + vfield['theta']**2)]))
    for a,c in zip(ax,components.keys()):
        sca(a)
        a.yaxis.set_major_locator(MultipleLocator(int(rlim/4)))
        a.xaxis.set_major_locator(MultipleLocator(int(rlim/4)))
        if a != ax[0]:
            a.yaxis.set_major_formatter(NullFormatter())
        dv_field = (components[c]-vkep)/vkep
        levels = np.linspace(-1.4,1.4,11)
        im=plot_slice(output,rho=dv_field,log=False,levels=levels, cmap = 'Spectral')
        for j in np.linspace(np.radians(model.env['theta_min']),pi/2.,8):
            streamline = model.stream(th0=j,shock=False)
            r = streamline['path'][0]
            th = streamline['path'][1]
            l = a.plot(r*np.sin(th),r*np.cos(th),color='white',ls='dashed',lw=1) #plot streamlines
        a.set_title(c,fontsize=11)
        a.set_xlim(1,rlim)
        a.set_ylim(1,rlim)
    colorbar(im,ax=ax,label=r'$\delta v \ \mathrm{[v_k]}$')
    ax[1].set_xlabel('R [au]')
    ax[0].set_ylabel('z [au]')
    
    

def plot_dustRT(output,rlim=400): #convenience function to plot the dust radiative transfer real quick
    model = output.m
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    X,Z = model.make_rz()
    X = X[:,:,0]
    Z = Z[:,:,0]
    f,ax = subplots(1,2,constrained_layout=True,dpi=100)
    f.set_size_inches(6,4)
    if output.rho == {}:
        output.read_rho()
    if output.T == {}:
        output.read_Tdust()
    toplot_rho = np.log10(output.calc_rho2D('dust')).T
    toplot_T = output.calc_T2D('dust').T
    
    c = ax[0].contourf(X, Z, toplot_rho,levels=np.arange(model.rhomin,model.rhomax,1),extend='both')
    ax[0].set_xlabel('r [au]')
    ax[0].set_ylabel('z [au]')
    ax[0].set_title('Input Model Dust Density',fontsize=11)
    cb = colorbar(c,ax=ax[0],location='bottom',aspect=10,ticks=np.array(c.levels[::2]).astype(int))
    cb.set_label(r'$\rho$ [$\mathrm{g/cm^{-3}}$]')
    
    smooth_T = toplot_T
    c2 = ax[1].contourf(X, Z, smooth_T,levels=np.linspace(5,125,19),cmap='twilight_shifted',extend='both')
    c3 = ax[1].contour(X, Z, smooth_T,levels=[10,20,30,60],colors='black',linewidths=1,linestyles=['solid','dashed','dashdot','dotted'])
    c4 = ax[0].contour(X, Z, smooth_T,levels=[10,20,30,60],colors='white',linewidths=1,linestyles=['solid','dashed','dashdot','dotted'])
    ax[1].plot(model.r, model.H(model.r),color='white',lw=1) #plot disk scale height
    fmt = {}
    for l in c3.levels:
        fmt[l] = str(int(l))
    ax[1].clabel(c3, c3.levels,fmt=fmt ,inline=True,fontsize=8)
    ax[0].clabel(c4, c4.levels,fmt=fmt ,inline=True,fontsize=8)
    
    streamline = model.stream(th0=np.radians(model.env['theta_min']),shock=False)
    r = streamline['path'][0]
    th = streamline['path'][1]
    ax[1].plot(r*np.sin(th),r*np.cos(th),color='white',ls='dashed',lw=1)
    
    ax[1].set_xlabel('r [au]')
    ax[1].set_title('Output Dust Temperature',fontsize=11)
    cb2 = colorbar(c2,ax=ax[1],location='bottom',ticks=np.array(c2.levels[::3]).astype(int),aspect=10)
    cb2.set_label(r'$T \ \mathrm{[K]}$')
    
    for a in ax:
        a.set_xlim(0,rlim)
        a.set_ylim(0,rlim)
        a.xaxis.set_major_locator(MultipleLocator(int(rlim/5)))
        a.yaxis.set_major_locator(MultipleLocator(int(rlim/5)))
    
    
def plot_flux_components(output,fname=''):
    model = output.m
    if os.getcwd() != model.outdir:
        os.chdir(model.outdir)
    accrate = model.star['accrate']
    wav,freq = output.wav, output.freq
    wav,fnu_bb = calc_blackbody_spectrum(model,wav = wav)
    wav,fnu_uv = calc_accretion_spectrum_uv(model,wav=wav,accrate=accrate)
    
    if model.rad['xray'] == True:
        wav,fnu_xray = calc_accretion_spectrum_xray(model,wav=wav,accrate=accrate)
    else:
        fnu_xray = np.zeros_like(fnu_bb)
        
    flya_peak = model_Lya_line(model,wav=wav,fLya=model.rad['fLya'])
    fnu_uv[np.argmin(np.abs(wav - lam_lya))] += flya_peak
    
    fnu_tot = fnu_bb + fnu_uv + fnu_xray
    
    f,ax= subplots(1,constrained_layout=True,dpi=100)
    f.set_size_inches(6,3)
    
    ax.loglog(wav,fnu_bb*freq,label='Star',lw=2,color='navy')
    ax.loglog(wav, fnu_uv*freq,label='UV Accretion',lw=2,color='skyblue')
    
    if model.rad['xray'] == True:
        ax.loglog(wav, fnu_xray*freq, label='X-ray Accretion',lw=2,color='violet')
        ax.set_xlim(xray_min,1e3)
    else:
        ax.set_xlim(xray_max,1e3)
        
    ax.loglog(wav,  freq*fnu_tot, label='total',color='black',ls='dashed',lw=1)
    ax.set_ylim(np.amax(freq*fnu_bb)*1e-8 , np.amax(freq*fnu_bb)*5)
    ax.set_xlabel(r'$\log \ \lambda \ \mathrm{[\mu m]}$')
    ax.set_ylabel(r'$\log \ \nu F_{\nu} \ \mathrm{[ergs \ cm^{-2} \ s^{-1}]}$')
    ax.set_title(r'$T_{\rm s}$' + ' = {} K,'.format(model.star['Ts']) + r' $\ \dot{M}_{\rm s}$' + '= {} '.format(model.star['accrate']) + r'$\mathrm{M_{\odot} \ yr^{-1}}$',size=12)
    ax.legend(loc=1)
    
    
    
def plot_shock_model(output):
    model=output.m
    colors=cm.get_cmap('twilight_r',20)
    thstart = np.radians(model.env['theta_min'])
    f,ax = subplots(2,4,constrained_layout=True,dpi=100)
    f.set_size_inches(6,5)
    Rc = model.env['Rc']
    Rd = model.disk['Rdisk'][0]
    Rin = model.env['Rc']*np.sin(thstart)**2
    halfstream = (thstart + pi/2)/2.
    Rhalf = model.env['Rc']*np.sin(halfstream)**2
    for j,k in zip(np.linspace(thstart+0.01,pi/2-0.01,20),range(20)):
        streams = model.stream(th0=j,shock=True)
        rland = streams['path'][0][-1]
        if j < halfstream:
            a = ax[0,:]
        else:
            a = ax[1,:]
        r = streams['path'][0]
        th = streams['path'][1]
        z = r*np.cos(th)
        rp = r*np.sin(th)
        rhod = model.rho_disk_profile(r=r*np.sin(th),z=r*np.cos(th))
        #Tcrit = 130*((streams['rho']+rhod)/(mu*mh*1e10))**(0.3)
        cs0 = model.cs(rland)
        a[0].loglog(z, streams['rho'] + rhod,color=colors(k))
        a[1].loglog(z,(streams['Tg']),color=colors(k))
        a[2].loglog(z,np.sqrt(streams['v'][0]**2+streams['v'][1]**2)/cs0,color=colors(k))
        a[3].loglog(z, r,color=colors(k))
        a[0].set_ylim(model.env['rho_amb']*1e6,10**(model.rhomax-1))
        a[1].set_ylim(5,5e3)
        a[2].set_ylim(0.05,50)
        a[3].set_ylim(1,2*Rc)
        for aa in a:
            aa.set_xlim(1,2*Rc)
    for aa in ax[0,:]:
        aa.axvline(model.H(Rc),ls='dotted',lw=1,color='gray')
        aa.axvline(10*model.H(Rc),ls='dotted',lw=1,color='gray')
    for aa in ax[1,:]:
        aa.axvline(model.H(Rc),ls='dotted',lw=1,color='gray')
        aa.axvline(10*model.H(Rc),ls='dotted',lw=1,color='gray')
        
    ax[-1,0].set_xlabel('z [au]')
    ax[-1,1].set_xlabel('z [au]')
    ax[-1,2].set_xlabel('z [au]')
    ax[-1,3].set_xlabel('z [au]')
    
    ax[0,0].set_title(r'$\rho_{\rm in} \ \mathrm{[g \ cm^{-3}]}$')
    ax[0,1].set_title(r'$T_{\rm gas}$ [K]')
    ax[0,2].set_title(r'$v_{\rm in} \ \mathrm{[c_s]}$')
    ax[0,3].set_title(r'$r \ \mathrm{[au]}$')
    ax[0,0].set_ylabel(r'$R_{\rm land}$'+' < {} au'.format(int(Rhalf)))
    ax[1,0].set_ylabel(r'$R_{\rm land}$'+' > {} au'.format(int(Rhalf)))
    f.suptitle('Shock Model: ' + r'$R_{\rm in}$ =' + '{} au,'.format(int(Rin))+ r'$R_{\rm c}$ =' + '{} au,'.format(int(Rc)) +
              r'$R_{\rm d}$ =' + '{} au,'.format(int(Rd)))
    norm = mpl.colors.Normalize(vmin=np.degrees(thstart), vmax=90)
    colorbar(cm.ScalarMappable(norm=norm, cmap='twilight_r'),location='bottom',ax=ax[0,:],label=r'streamline $\theta_0$ [deg]',shrink=0.5,aspect=10)
    
def plot_streams_polar(output,rlim=400):
    model = output.m
    np0 = len(model.phi)
    nstream = model.env['nstreams']
    dphi_frac = model.env['stream_frac']
    npstream = int(dphi_frac*np0)/nstream #how many phi0 per stream
    stream_phi = np.array([(np.arange(0,npstream,1) + i*np0/nstream).astype(int) for i in range(nstream)]).flatten()
    thstart = np.radians(model.env['theta_min'])

    dphi = np.gradient(model.coords[2])[0]
    subset_phi = model.phi[stream_phi]
    index_phi = stream_phi
    
    for p0 in subset_phi:
        r_s = np.array([])
        th_s = np.array([])
        phi_s = np.array([])
        streamline = model.stream(th0=thstart,p0=p0,shock=False)
        r_s = np.append(r_s,streamline['path'][0])
        th_s = np.append(th_s,streamline['path'][1])
        phi_s = np.append(phi_s,streamline['path'][2])
    
        x_s = r_s*np.sin(th_s)*np.cos(phi_s)
        y_s = r_s*np.sin(th_s)*np.sin(phi_s)
    
        plot(x_s,y_s,lw=1,color='black')
    axis('scaled')
    xlim(-rlim,rlim)
    ylim(-rlim,rlim)


def plot_kappa(output,fluid=1,filename='',**kwarg):
    model = output.m
    lam,kabs,kscat,g= read_kappa(model,fluid=fluid,filename=filename)
    loglog(lam,kabs,lw=2,**kwarg,label=filename+ r': $\kappa_{abs}$')
    loglog(lam,kscat,lw=2,ls='--',**kwarg,label=filename + r': $\kappa_{sca}$',alpha=0.5)
    ax=gca()
    ax.set_xlabel(r'$\lambda \ \mathrm{[\mu m]}$')
    ax.set_ylabel(r'$\kappa \ \mathrm{[cm^{2} \ g^{-1}]}$')
    ax.set_xlim(1.2e-4,1.2e4)
    ax.set_ylim(1e-4,)
    
def plot_opacities(output,filenames = None,xray=False):
    f,ax = subplots(1,constrained_layout=True,dpi=100)
    f.set_size_inches(6,3)
    colors = ['teal','tomato','skyblue','goldenrod']
    if filenames is not None:
        pops = {}
        for dust,c in zip(filenames,colors):
            pops[dust] = c
    else:
        pops = {'dust-1':'teal','dust-2':'tomato'}  
    for dust in pops.keys():
        filename = dust
        if xray == True:
            filename = dust + 'x'
            header = False
        plot_kappa(output,filename=filename,color=pops[dust])
    ax.set_title('Dust Population Opacities',fontsize=11)
    ax.legend()
    
def save_setup_fig(output,rlim=400):
    from matplotlib.backends.backend_pdf import PdfPages
    os.chdir(output.m.outdir)
    with PdfPages('disk_setup.pdf') as pp:
        plot_components(output,0,rlim=rlim)
        f = gcf()
        f.savefig(pp,format='pdf')
        plot_components(output,1,rlim=rlim)
        f = gcf()
        f.savefig(pp,format='pdf')
        plot_components(output,2,rlim=rlim)
        f = gcf()
        f.savefig(pp,format='pdf')
        plot_velocities(output,rlim=rlim)
        f = gcf()
        f.savefig(pp,format='pdf')
        plot_flux_components(output)
        f = gcf()
        f.savefig(pp,format='pdf')
        plot_opacities(output)
        f = gcf()
        f.savefig(pp,format='pdf')
        if output.m.env['shock'] == True:
            plot_shock_model(output)
            f = gcf()
            f.savefig(pp,format='pdf')
            
def save_dustRT(output):
    from matplotlib.backends.backend_pdf import PdfPages
    os.chdir(output.m.outdir)
    with PdfPages('dustRT.pdf') as pp:
        plot_dustRT(output,rlim=output.m.grid['max'][0])
        f = gcf()
        f.savefig(pp,format='pdf')
        plot_dustRT(output,rlim=output.m.grid['max'][0]/5)
        f = gcf()
        f.savefig(pp,format='pdf')

    
def plot_radiation_field(output,rlim=400, field='uv',levels=None):
    model = output.m
    if field not in output.J.keys():
        output.calc_Jint(field=field)
    f,ax = subplots(1, constrained_layout=True,dpi=100)
    f.set_size_inches(3,2)
    sca(ax)
    if levels is None:
        levels = np.linspace(1,12,12)
    im = plot_slice(output,rho=np.average(4*pi*output.J[field]['J_phot'],axis=-1),log=True,levels=levels,cmap='BuPu_r')
    cb=colorbar(im,ax=ax)
    
    cb.set_label(r'$\log$ photons [$\mathrm{cm^{-2} \ s^{-1}}$]')
    ax.set_ylabel('z [au]')
    
    ax.set_title(field)
    ax.set_xlabel('r [au]')
    ax.set_xlim(0,rlim)
    ax.set_ylim(0,rlim)
        
        
def plot_all_temp(output,rlim=200):
    pp = {'levels':np.linspace(0.5,3,51),'cmap':'twilight_shifted','extend':'both'}
    f,ax = subplots(2,2,constrained_layout=True)
    f.set_size_inches(6,6)
        
    for a, temp in zip([ax[0,0],ax[0,1],ax[1,0],ax[1,1]], ['uv','shock','dust','gas']):
        sca(a)
        toplot_T = output.calc_T2D(temp)
        im=plot_slice(output,rho=toplot_T.T,**pp)
        title('T [{}]'.format(temp))
        a.axis('scaled')
        a.set_xlim(0,rlim)
        a.set_ylim(0,rlim)
    

    cb = colorbar(im,ax=ax,location='top',shrink=0.75)
    cb.set_label(r'$\log$ T [K]')
    ticks = np.log10([5,10,25,50,100,250,500,1000,2000,4000])
    cb.set_ticks(ticks)
    cb.set_ticklabels(['5','10','25','50','100','250','500','1000','2000','4000'])
    
    ax[0,0].set_ylabel('z [au]')
    ax[1,0].set_ylabel('z [au]')
    ax[1,0].set_xlabel('r [au]')
    ax[1,1].set_xlabel('r [au]')

        
def save_heRT(output):
    from matplotlib.backends.backend_pdf import PdfPages
    os.chdir(output.m.outdir)
    with PdfPages('heRT.pdf') as pp:
        plot_radiation_field(output,rlim=output.m.grid['max'][0],field='uv')
        f = gcf()
        f.savefig(pp,format='pdf')
        
        if output.m.rad['xray'] == True:
            plot_radiation_field(output,rlim=output.m.grid['max'][0],field='xray')
            f = gcf()
            f.savefig(pp,format='pdf')
        plot_all_temp(output,rlim=output.m.grid['max'][0])
        f = gcf()
        f.savefig(pp,format='pdf')
        plot_opacities(output,xray=output.m.rad['xray'])
        f = gcf()
        f.savefig(pp,format='pdf')
        
