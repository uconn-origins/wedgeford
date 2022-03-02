from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy.interpolate import griddata
from models.units import *

parent_dir = os.getcwd()
models_dir = parent_dir+'/models/'

class model:
    def __init__(self,stellar_params, disk_params, envelope_params,grid_params, dust_params, RT_params, outdir = '/m1_test/'):
        self.star = stellar_params
        self.disk = disk_params
        self.env = envelope_params
        self.grid = grid_params
        self.dust = dust_params
        self.rad = RT_params
        self.outdir = parent_dir+outdir
        self.coords = [0,0,0]
        self.rhomin = np.log10(self.env['rho_amb'])
        self.rhomax = -11
        
        try:
            os.mkdir(self.outdir)
        except:
            print('directory exists - will overwrite current model if you write to it!')
        for dim in [0,1,2]:
            N,gmin,gmax,spacing = [self.grid['N'][dim], self.grid['min'][dim], self.grid['max'][dim], self.grid['spacing'][dim]]
            if spacing == 'log':
                self.coords[dim] = np.logspace(np.log10(gmin),np.log10(gmax),N+1)
            elif spacing == 'lin':
                self.coords[dim] = np.linspace(gmin,gmax,N+1)
            else:
                print('Defaulting to linear spacing')
                self.coords[dim] = np.linspace(gmin,gmax,N+1) # these are cell faces
            dim += 1
        self.r = (self.coords[0][1:] + self.coords[0][:-1]) / 2. # these are cell centers
        self.theta = (self.coords[1][1:] + self.coords[1][:-1]) / 2.
        self.phi = (self.coords[2][1:] + self.coords[2][:-1]) / 2.
        
        
        params = rpy.params.radmc3dPar()
        if os.path.exists(self.outdir+ 'problem_params.inp') != True:
            params.loadDefaults() #loads default parameter dictionary
        else:
            params.readPar()
        params.ppar['mstar'] = [self.star['Ms']*Msun] # sets the paramaters to the ones passed to the model
        params.ppar['rstar']= [self.star['Rs']*Rsun]
        params.ppar['tstar'] = [self.star['Ts']*1.0]
        params.ppar['accrate'] = self.star['accrate']*Msun/yr
        params.ppar['starsurffrac'] = self.star['f']
        self.radmcpar = params.ppar
        
    def T(self,R):
        return self.star['Ts']*np.sqrt((self.star['Rs']*Rsun)/(R*AU))
    
    def Rsub(self): #dust sublimation radius at T = 1500 K
        Tsub = 1500.
        T_R = self.T(R=self.r)
        try:
            rsub = np.amin(self.r[T_R>=Tsub])
        except:
            rsub = 0.0
        return rsub
        
    def cs(self,R):
        return 5862*np.sqrt(self.T(1))*(R**(-0.25))
    
    def vk(self,R):
        return sqrt((Gconv*self.star['Ms']*Msun)/(R*AU))
    
    def H(self,R,fluid=0):
        def H_gas(R): #computes the pressure scale height for temp. distribution.
            return self.cs(R)*R/(self.vk(R))
        H0 = self.disk['H0'][fluid-1]
        R0 = self.disk['R0'][fluid-1]
        FI = self.disk['fi'][fluid-1]
        if fluid == 0:
            return H_gas(R)
        else:
            return H_gas(R)*H0
        #elif fluid == 1 and H0 == None: 
            #return H_gas(R)
        #elif fluid > 0 and H0 > H_gas(R0):
            #return H_gas(R)
        #else:
            #return H0*(R/R0)**(1+FI)
        
    def make_grid(self):
        return np.meshgrid(self.r,self.theta,self.phi)
    
    def make_rz(self): #non-uniform spherical grid mapped onto cylindrical coordinates
        R,THETA,PHI = self.make_grid()
        R_CYL = np.clip(R*np.sin(THETA),a_min=1e-5,a_max=None)
        Z_CYL = R*np.cos(THETA)
        return R_CYL, Z_CYL
        
    def sig_profile(self,R,fluid=0): #fluid = 0 for gas, 1 for small dust (follows gas), 2 for large dust (settled)
        i = max(fluid,1)-1
        Rd = self.disk['Rdisk'][i]
        Rout = self.disk['Rout'][i]
        p = self.disk['p'][i]
        R0 = self.disk['R0'][i]
        frac = [1.-self.disk['Mfrac'][0]-self.disk['Mfrac'][1],self.disk['Mfrac'][0],self.disk['Mfrac'][1]]
        Mtot = self.disk['Mdisk']*frac[fluid]*Msun #assigns mass of disk component
        def sig_r(R):
            sig = (R/R0)**(p)*np.exp(-R/Rd)*np.exp(-(R/Rout)**2)
            sig[R<R0] = 0.
            return sig
        def sig_int(R):
            if len(np.shape(R)) < 2:
                sig = sig_r(R)
                return np.sum(sig*(R*AU)*(np.gradient(R))*AU)
            else:
                r = np.ravel(np.unique(R))
                sig = sig_r(r)
                return np.sum(sig*(r*AU)*np.gradient(r)*AU)
        sig_0 = (Mtot)/(2.*pi*sig_int(R=self.r))
        return sig_0*sig_r(R)
    
    def rho_midplane(self,R,fluid=0):
        rho_mid = np.clip(self.sig_profile(R,fluid=fluid)/(sqrt(pi*2.)*self.H(R,fluid=fluid)*AU),a_min=self.env['rho_amb'],a_max=None)
        return rho_mid
    
    def rho_disk_profile(self,r,z,fluid=0):
        rho_mid = self.rho_midplane(r,fluid=fluid)
        h_mid = self.H(r,fluid=fluid)
        rho_at_pt = rho_mid*np.exp(-0.5*(z/h_mid)**2)
        return rho_at_pt
    
    def rho_disk(self,fluid=0):
        R_CYL,Z_CYL = self.make_rz()
        rho_mid = self.rho_midplane(R_CYL,fluid=fluid)
        self.rhomax = np.log10(np.amax(rho_mid))
        h_mid = self.H(R_CYL,fluid=fluid)
        rho_vol = rho_mid*np.exp(-0.5*(Z_CYL/h_mid)**2)
        return np.clip(rho_vol,a_min=self.env['rho_amb'],a_max=None)
        
    def stream(self, th0=pi/2.,p0=0,shock=True): #calculates the streamline position and properties (density and velocity)
        thstart = np.radians(self.env['theta_min'])
        th = self.theta[self.theta>=th0]
        Rc = self.env['Rc']
        Min = self.env['Min']
        rho0 = self.env['rho_0']
        Rin = Rc*np.sin(np.radians(self.env['theta_min']))**2 #inner landing radius of the stream
        Mfac = np.sqrt(1. - sqrt(Rin/Rc)) #adjusts the total mass to fit within the chosen streams
        streamline = {}

        vwr=0.
        vwphi = 0.
        vwth = 0.

        v_r = np.array([])
        v_phi = np.array([])
        v_theta = np.array([])

        if th0 != pi/2. and th0 != 0.:
            cosa = (np.cos(th)/np.cos(th0))
            r = Rc*(np.sin(th0)**2)/np.clip((1. - cosa),a_min=1e-18,a_max=None)
            r[0] = self.env['Rmax']
            r[-1] = Rc*(np.sin(th0)**2)
        elif th0 <= 1e-5:
            r = np.ones_like(th)*self.r[0]
        else:
            r= np.ones_like(th)*Rc

        om = np.sqrt(Gconv*self.star['Ms']*Msun*(r*AU)**3)
        t1 = ((Min/Mfac)*Msun/yr)/(8*pi*om)
        t2 = np.sqrt(1+(np.cos(th)/np.cos(th0)))
        t3 = ((np.cos(th)/(2*np.cos(th0))) + (Rc/r)*(np.cos(th0)**2))
        rho1 = t1/(t2*t3)

        cosa = np.cos(th)/np.cos(th0)
        sina = np.sqrt(1-cosa**2)
        tanphi = (sina/cosa)/np.sin(th0)
        delphi = np.arctan(tanphi)
        ph = p0+delphi

        vkep = self.vk(r)
        if th0 <= thstart: # outside of opening angle
            v_r = np.append(v_r,vwr*np.ones_like(r))
            v_phi = np.append(v_phi,vwphi*np.ones_like(r))
            v_theta = np.append(v_theta,vwth*np.ones_like(r))
            rho1 = rho0*r**(-self.env['exf'])
        else:
            v_r = np.append(v_r, -vkep*np.sqrt(1+cosa))
            v_phi = np.append(v_phi,vkep*np.sqrt(1-cosa)*(np.sin(th0)/np.sin(th)))
            v_theta = np.append(v_theta,vkep*np.sqrt((1+cosa))*(np.cos(th0)-np.cos(th))/np.sin(th))

        if shock == True and th0 > thstart: 
            rhod = self.rho_disk_profile(r=r*np.sin(th),z=r*np.cos(th))
            Tstream = np.zeros_like(r)
            rho_m = 0.5*(rho1+rhod)
            surface = (rhod*rho_m)/(rho1**2) #calculate where the envelope should hit the disk surface to form a shock
            if np.amin(surface) < 1.0: 
                s = np.amax(np.where(surface <= 1.))
                rupper = r[s] #top of shock surface
                vs = (rho1[s]/(2*rho_m[s]))*np.sqrt(v_r**2 + v_theta**2)[s]
                cs0 = self.cs(rupper)
                Ts = (3/16.)*mu*mh*(vs)**2 
                rhos = 2*rho_m[s]
                Q0 = 1e-31*(rhos/(mu*mh))**2 * Ts**(1.5) #optically thin cooling limit for CO
                s0 = (1./AU)*(2.5*rhos*vs**2*((vs/4) - cs0))/Q0 #cooling decay coefficient = amt of cooling/peak cooling
                Q = Q0*np.exp(-(r[s:]/s0)) #cooling over rest of streamline
                v_cool = (Q-Q0*s0*AU)/(2.5*rhos*vs**2) + (vs/4) # deceleration as energy is radiated away
                rho_cool = rhos*vs/v_cool #compression
                rho1[s:] = np.clip(rho_cool - rhod[s:],a_min=self.env['rho_amb'],a_max=None)
                T_cool = (mu*mh/kb)*rhos*(vs**2)/(rho1[s:] + rhod[s:]) 
                v_theta[s:] *= (v_cool/vs)
                v_r[s:] *= (v_cool/vs)
                Tstream[s:] = T_cool

            streamline['Tg'] = Tstream
        streamline['path'] = [r,th,ph]
        streamline['rho'] = rho1
        streamline['v'] = [v_r,v_theta,v_phi]
        return streamline
    
    
    def rho_env(self,fluid=0):
        shock = self.env['shock']
        thstart = np.radians(self.env['theta_min'])
        rho0 = self.env['rho_0']
        R = np.array([])
        TH = np.array([])
        rho_tot = np.array([])
        for j in self.theta:
            streamline = self.stream(th0=j,shock=shock)
            R = np.append(R,streamline['path'][0])
            TH = np.append(TH,streamline['path'][1])
            rho_tot = np.append(rho_tot,streamline['rho'])
        R_CYL = R*np.sin(TH)
        Z_CYL = R*np.cos(TH)
        r_cyl, z_cyl = self.make_rz()
        RHO = griddata((R_CYL, Z_CYL), rho_tot, (r_cyl[:,:,0],z_cyl[:,:,0]), method='nearest')
        RHO = np.reshape(RHO,(np.shape(RHO)[0],np.shape(RHO)[1],1))
        rho_vol = np.repeat(RHO,len(self.phi),axis=2) #assumes axisymmetry
        if fluid == 0:
            return rho_vol + self.env['rho_amb']
        elif fluid == 1:
            return (rho_vol + self.env['rho_amb'])*self.env['d2g']
        else:
            return rho_vol*0 + self.env['rho_amb']
        
    def v_env(self):
        shock = self.env['shock']
        thstart = np.radians(self.env['theta_min'])
        rho0 = self.env['rho_0']
        R = np.array([])
        TH = np.array([])
        V = {}
        vr = np.array([])
        vth = np.array([])
        vphi = np.array([])
        for j in self.theta:
            streamline = self.stream(th0=j,shock=shock)
            R = np.append(R,streamline['path'][0])
            TH = np.append(TH,streamline['path'][1])
            vr = np.append(vr,streamline['v'][0])
            vth = np.append(vth,streamline['v'][1])
            vphi = np.append(vphi,streamline['v'][2])
        R_CYL = R*np.sin(TH)
        Z_CYL = R*np.cos(TH)
        r_cyl, z_cyl = self.make_rz()
        for stream,key in zip([vr,vth,vphi],['r','theta','phi']):
            V[key] = griddata((R_CYL, Z_CYL), stream, (r_cyl[:,:,0],z_cyl[:,:,0]), method='nearest')
            V[key] = np.reshape(V[key],(np.shape(V[key])[0],np.shape(V[key])[1],1))
        return V
    
    def v_disk(self):
        R_CYL,Z_CYL = self.make_rz()
        R = np.sqrt(R_CYL**2 + Z_CYL**2)
        V = {}
        V['phi'] = self.vk(R_CYL)
        vin = -self.star['accrate']/(2*pi*R_CYL*np.sqrt(2*pi)*self.H(R_CYL)*AU*self.rho_disk(fluid=0))
        V['r'] = vin*R_CYL/R
        V['theta'] = vin*Z_CYL/R
        return V
    
    def v_embedded(self):
        V = {}
        vdisk = self.v_disk()
        venv = self.v_env()
        for key in ['r','theta','phi']:
            V[key] = (self.rho_disk(fluid=0)*vdisk[key] + self.rho_env(fluid=0)*venv[key])/self.rho_embedded(fluid=0)
        return V
            
        
           
    def rho_embedded(self,fluid=0):
        rho_d = self.rho_disk(fluid=fluid)
        if fluid < 2: # for gas and small dust
            rho_e = self.rho_env(fluid=fluid)
        else: # for large dust only
            rho_e = np.zeros_like(rho_d)
        floor = self.env['rho_amb'] 
        factor = min(fluid,1)*self.env['d2g']
        if factor != 0: 
            floor *= factor
        if fluid < 2: #for gas and small dust
            rho_both = np.clip(rho_d + rho_e,a_min=floor,a_max=None)
        else: # for large dust only
            rho_both = rho_d + rho_e
        if fluid != 0: #no dust inside sublimation radius
            R,THETA,PHI = self.make_grid()
            rsub = self.Rsub()
            rho_both[R <= rsub] = floor
        return rho_both
    
    def v_transform(self,v): #spherical to cartesian
        R,THETA,PHI = self.make_grid()
        V = {}
        for key in ['r','theta','phi']:
            v[key] = np.repeat(v[key],len(self.phi),axis=2) #assumes axisymmetry
        V['x'] = v['r']*np.sin(THETA)*np.cos(PHI) + v['theta']*np.cos(THETA)*np.cos(PHI) - v['phi']*np.sin(PHI)
        V['y'] = v['r']*np.sin(THETA)*np.sin(PHI) + v['theta']*np.cos(THETA)*np.sin(PHI) + v['phi']*np.cos(PHI)
        V['z'] = v['r']*np.cos(THETA) + v['theta']*np.sin(THETA)
        return V

    
    

    
    
    
    
    
    