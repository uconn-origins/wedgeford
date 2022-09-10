from pylab import *
import radmc3dPy as rpy
import numpy as np
import os
from scipy import interpolate
from . units import *


class model_:
    def __init__(self,params, outdir = '/m1_test/'):
        stellar_params = {'Ms': 1, 'Rs': 2.0, 'Ts': 5000, 'accrate':1e-7,'f':0.01}
        disk_params = {'Mdisk': 0.06, 'Mfrac': [0.01,0.01],'R0':[5,5], 'H0':[1,0.3], 'p':[-1,-1], 'Rdisk':[125,125],
                      'Tfac':1, 'q':0.5, 'hydro':[0,0,0]}
        envelope_params = {'Min': 1e-6, 'Rc':125, 'rho_amb':1e-25, 'rho_0': 3e-22,'theta_min': 25,'exf':0.25,'Rmax':1.5e4, 'd2g': 0.01, 'shock':False, 'nstreams': 1, 'stream_frac':1}
        grid_params = {'N':[180,90,48], 'min':[0.1,pi/16.,0], 'max':[400,pi/2.,2*pi], 'spacing':['log','lin','lin']}
        dust_params = {'rho_si':3.1518, 'amin_chem':0.06, 'amax_ism': 1.0, 'amin': [0.005,0.005], 'amax': [1,1e3], 'apow': [3.5,3.5]}
        RT_params = {'cr_model': 'ssx','zetacr': 1.3e-17, 'G0':1, 'viscous_heating':True, 'fLya': 1e-4}
        for key in params.keys():
            if key in stellar_params.keys():
                stellar_params[key] = params[key]
            elif key in disk_params.keys():
                disk_params[key] = params[key]
            elif key in envelope_params.keys():
                envelope_params[key] = params[key]
            elif key in dust_params.keys():
                dust_params[key] = params[key]
            elif key in RT_params.keys():
                RT_params[key] = params[key]
            elif key in grid_params.keys():
                grid_params[key] = params[key]
            
        self.star = stellar_params
        self.disk = disk_params
        self.env = envelope_params
        self.grid = grid_params
        self.dust = dust_params
        self.rad = RT_params
        self.print_params()
        from pathlib import Path
        p = Path().resolve()
        self.parent_dir = str(p).split('wedgeford')[0] + 'wedgeford/'
        self.models_dir = self.parent_dir+'models/'
        self.outdir = self.parent_dir+outdir
        self.coords = [0,0,0]
        self.rhomin = int(np.log10(self.env['rho_amb']))
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
        
        tgas = np.load(self.models_dir+'templates/tgas.npy')
        log_n = np.cumsum(np.ones(49)*0.125) - 0.125 #second axis
        log_F = np.linspace(-3.29588079,-3.29588079+ 57*0.125, 57) - np.log10(G0) #first axis
        self.T_heat = interpolate.interp2d(log_n,log_F,tgas)
        if self.env['stream_frac'] < 1:
            self.isolate_stream(dphi_frac = self.env['stream_frac'] , nstream = self.env['nstreams'])
        
    
    def print_params(self):
        print('stellar_params:',self.star)
        print('disk_params:',self.disk)
        print('envelope_params:',self.env)
        print('grid_params:',self.grid)
        print('dust_params:',self.dust)
        print('RT_params:',self.rad)
    
    def T(self,R):
        Tfac = self.disk['Tfac']
        q = self.disk['q']
        return Tfac*self.star['Ts']*((0.5*self.star['Rs']*Rsun)/(R*AU))**(q)
        
    def cs(self,R):
        return np.sqrt(kb*self.T(R)/(mu*mh))
    
    def vk(self,R):
        return sqrt((Gconv*self.star['Ms']*Msun)/(R*AU))
    
    def H(self,R,fluid=0):
        def H_gas(R): #computes the pressure scale height for temp. distribution.
            return self.cs(R)*R/(self.vk(R))
        H0 = self.disk['H0'][fluid-1]
        if fluid == 0:
            return H_gas(R)
        else:
            return H_gas(R)*H0
        
    def disk_transport(self,setnew=False):
        Rc = self.env['Rc']
        Min = self.env['Min']
        Rin = Rc*np.sin(np.radians(self.env['theta_min']))**2
        tauc = Rc**(3/2.)
        sigc = self.sig_profile(self.r)[self.r <= Rc][-1]
        sigin = (Min*Msun)/(4*pi*((Rc*AU)**2-(Rin*AU)**2)) * tauc
        logalpha = min(2*(np.log10(sigin) - np.log10(sigc)) - 0.5, -1.75)
        self.disk['alpha'] = 10**logalpha
        Ms = self.star['Ms']
        R0 = self.disk['R0'][0]
        vc = np.sqrt(Gconv*Ms*Msun/(R0*AU))
        h0 = self.H(R0)
        rho0 = self.rho_disk_profile(r=np.array([R0,R0]),z=np.array([0,0]))[0]
        Macc = sqrt(18*pi**3)*(10**logalpha)*vc*rho0*(h0*AU)**3/(R0*AU)/(Msun/yr)
        print('Accretion rate from infall: {} Msun/yr'.format(Macc))
        if setnew == True:
            self.star['accrate'] = Macc
            self.radmcpar['accrate'] = self.star['accrate']*Msun/yr
        
    def make_grid(self,quadrant=False):
        if quadrant == True:
            return np.meshgrid(self.r,self.theta,self.phi[:int(len(self.phi)/4)])
        else:
            return np.meshgrid(self.r,self.theta,self.phi)
        
        
    def dvol(self):
        R,THETA,PHI = self.make_grid()
        DR = np.gradient(R,axis=1)
        DTH = np.gradient(THETA,axis=0)
        DPH = np.gradient(PHI,axis=-1)
        dvol = (R*AU)**2 * (DR*AU) * DTH *np.sin(THETA) * DPH
        return dvol
    
    def make_rz(self): #non-uniform spherical grid mapped onto cylindrical coordinates
        R,THETA,PHI = self.make_grid()
        R_CYL = R*np.sin(THETA)
        Z_CYL = R*np.cos(THETA)
        return R_CYL, Z_CYL
    
    def make_xyz_quadrant(self):
        R,THETA,PHI = self.make_grid(quadrant=True)
        X_CYL = R*np.sin(THETA)*np.cos(PHI)
        Y_CYL = R*np.sin(THETA)*np.sin(PHI)
        Z_CYL = R*np.cos(THETA)
        return X_CYL, Y_CYL, Z_CYL
        
    def sig_profile(self,R,fluid=0): #fluid = 0 for gas, 1 for small dust (follows gas), 2 for large dust (settled)
        i = max(fluid,1)-1 #small dust follows gas
        Rd = self.disk['Rdisk'][i]
        p = self.disk['p'][i]
        R0 = self.disk['R0'][i]
        frac = [1.-self.disk['Mfrac'][0]-self.disk['Mfrac'][1],self.disk['Mfrac'][0],self.disk['Mfrac'][1]]
        Mtot = self.disk['Mdisk']*frac[fluid]*Msun #assigns mass of disk component
        def sig_r(R):
            sig = (R/R0)**(p)*np.exp(-(R/Rd)**(2.+p))
            sig[R<R0] = 0.
            #sig[R>Rout] *= (Rd/R[R>Rout])**2
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
        if len(self.disk['hydro'][fluid]) < 2:
            R_CYL,Z_CYL = self.make_rz()
            rho_mid = self.rho_midplane(R_CYL,fluid=fluid)
            self.rhomax = np.log10(np.amax(rho_mid)/5)
            h_mid = self.H(R_CYL,fluid=fluid)
            rho_vol = rho_mid*np.exp(-0.5*(Z_CYL/h_mid)**2)
            return np.clip(rho_vol,a_min=self.env['rho_amb'],a_max=None)
        else:
            return self.disk['hydro'][fluid]
        
    def stream(self, th0=pi/4.,p0=0,shock=True): #calculates the streamline position and properties (density and velocity)
        thstart = np.radians(self.env['theta_min'])
        thmin = self.grid['min'][1]
        #th = self.theta[self.theta>=th0]
        th = np.arccos(np.linspace(np.cos(th0),1e-6, 360))
        Rc = self.env['Rc']
        Min = self.env['Min']
        rho0 = self.env['rho_0']
        Rin = Rc*np.sin(thstart)**2 #inner landing radius of the stream
        Mfac = np.sqrt(1. - sqrt(Rin/Rc)) #adjusts the total mass to fit within the chosen streams
        streamline = {}

        vwr=0.
        vwphi = 0.
        vwth = 0.

        v_r = np.array([])
        v_phi = np.array([])
        v_theta = np.array([])

        cosa = (np.cos(th)/np.cos(th0))
        amax = 1.0 - (Rc*(np.sin(th0)**2))/self.env['Rmax']
        cosa = np.clip(cosa, a_min = None, a_max = amax)
        sina = np.sqrt(1-cosa**2)
        r = Rc*(np.sin(th0)**2)/(1. - cosa)
        r[0] = self.env['Rmax']
        
        
        tanphi = (sina/cosa)/np.sin(th0)
        delphi = np.arctan(tanphi)
        ph = p0+delphi
        
        #r = np.append(r,Rc*(np.sin(th0)**2))
        #th = np.append(th,pi/2.)
        #cosa = np.append(cosa,0)
        #ph = np.append(ph, p0+pi/2.)
        
        om = np.sqrt(Gconv*self.star['Ms']*Msun*(r*AU)**3)
        t1 = ((Min/Mfac)*Msun/yr)/(8*pi*om)
        #t2 = np.sqrt(1+(np.cos(th)/np.cos(th0)))
        t2 = np.sqrt(1 + cosa)
        #t3 = ((np.cos(th)/(2*np.cos(th0))) + (Rc/r)*(np.cos(th0)**2))
        t3 = 0.5*cosa + (Rc/r)*(np.cos(th0)**2)
        rho1 = t1/(t2*t3) * (1./self.env['stream_frac'])
        
        vkep = self.vk(r)
        cs = self.cs(r)
        if th0 <= thstart: # outside of opening angle
            v_r = np.append(v_r,vwr*np.ones_like(r))
            v_phi = np.append(v_phi,vwphi*np.ones_like(r))
            v_theta = np.append(v_theta,vwth*np.ones_like(r))
            rho1 = rho0*r**(-self.env['exf'])
        else:
            v_r = np.append(v_r, -vkep*np.sqrt(1+cosa))
            v_phi = np.append(v_phi,vkep*np.sqrt(1-cosa)*(np.sin(th0)/np.sin(th)))
            v_theta = np.append(v_theta,vkep*np.sqrt((1+cosa))*(np.cos(th0)-np.cos(th))/np.sin(th))
        
        rhod = self.rho_disk_profile(r=r*np.sin(th),z=r*np.cos(th))
        rho_m = rho1+rhod
        u_r = (v_r*rho1)/rho_m
        u_phi = (v_phi*rho1 + vkep*rhod)/rho_m
        u_theta = (v_theta*rho1)/rho_m
        Tstream = self.T(R=r)
        
        if shock == True and th0 > thstart: 
            #ds = np.abs(np.gradient(r) + np.gradient(th)*r)#+ r*np.sin(th)*np.gradient(ph))
            #calculate where shock should start with highest velocity transition
            vs = np.sqrt(u_r**2 + u_theta**2)
            ds = np.sqrt(np.gradient(r)**2 + (np.gradient(th)*r)**2)
            dv = np.gradient(vs)
            s = np.argmax(-dv/cs)
            if np.any(-dv > 0) and (vs[s]/cs[s] >= 1.0) : # if deceleration + supersonic
                vs0 = np.sqrt(u_r[s]**2 + u_theta[s]**2)
                #cs0 = self.cs(r[s]*np.sin(th[s]))
                cs0 = self.cs(r[s])
                Ts = (3/16.)*(mu*mh/kb)*(vs0)**2
                rhos = rho_m[s]
                Q0 = 1e-31*(rhos/(mu*mh))**2 * Ts**(1.5) #optically thin cooling limit for CO #ergs/cm^3/s
                s0 = (1./AU)*(2.5*rhos*(vs0**2)*(vs0/4. - cs0))/Q0 #cooling decay coefficient = amt of cooling/peak 
                Q = Q0*np.exp(-(r[s:]-r[s])/s0) #cooling over the streamline (radially)
                v_cool = ((Q-Q0)*s0*AU)/(2.5*rhos*vs0**2) + (vs0/4.) # deceleration as energy is radiated away
                v_cool = np.clip(v_cool, a_max = vs0/4., a_min = cs0)
                rho_cool = rhos*(vs0/v_cool) #compression of inward flow
                rho_m[s:] = (rho_m[s:] - rhos) + rho_cool #compression of inwards flow
                T_cool = (mu*mh/kb)*rhos*(vs0**2)/(rho_m[s:]) #cooling region as disk gets denser
                Tstream[s:] = (T_cool*rho_cool + Tstream[s:]*(rho_m[s:]-rho_cool))/rho_m[s:]
                u_theta[s:] = vs0/(rho_m[s:]/rhos)
                u_r[s:] = vs0/(rho_m[s:]/rhos)
            v_r = u_r
            v_phi = u_phi
            v_theta = u_theta
            rho1 = rho_m - rhod
        streamline['Tg'] = Tstream
        streamline['path'] = [r,th,ph]
        streamline['rho'] = rho1
        streamline['v'] = [v_r,v_theta,v_phi]
        return streamline
    
    def solve_envelope(self,prop='rho',shock=None):
        if shock == None:
            shock = self.env['shock']
        thstart = np.radians(self.env['theta_min'])
        axisymmetric = True
        prop_tot = np.array([])
        if prop == 'v':
            PROP = {}
            vr = np.array([])
            vth = np.array([])
            vphi = np.array([])
        R = np.array([])
        TH = np.array([])
        PHI = np.array([])
        if axisymmetric == True:
            np0 = 1
        acosj = np.arccos(np.linspace(np.cos(self.grid['min'][1]),1e-6, 180))
        for j in acosj:
            streamline = self.stream(th0=j,p0=0,shock=shock)
            R = np.append(R,np.repeat(streamline['path'][0],np0,axis=-1))
            TH = np.append(TH,np.repeat(streamline['path'][1],np0,axis=-1))
            phi_ =  np.array([streamline['path'][2] + j for j in np.linspace(0,pi/2,np0)]).flatten()
            PHI = np.append(PHI,phi_)
            if np.ndim(streamline[prop]) < 2:
                prop_tot = np.append(prop_tot,np.repeat(streamline[prop],np0,axis=-1))
            else:
                vr = np.append(vr,np.repeat(streamline['v'][0],np0,axis=-1))
                vth = np.append(vth,np.repeat(streamline['v'][1],np0,axis=-1))
                vphi = np.append(vphi,np.repeat(streamline['v'][2],np0,axis=-1))
                prop_tot = [vr,vth,vphi]
        X = R*np.sin(TH)*np.cos(PHI)
        Y = R*np.sin(TH)*np.sin(PHI)
        R_CYL = np.sqrt(X**2 + Y**2)
        Z_CYL = R*np.cos(TH)
        r_cyl, z_cyl = self.make_rz()
        #x, y, z = self.make_xyz_quadrant()
        if np.ndim(prop_tot) < 2:
            #PROP = interpolate.griddata((X,Y,Z_CYL), prop_tot, (x,y,z), method='nearest',rescale=True)
            PROP = interpolate.griddata((R_CYL,Z_CYL), prop_tot, (r_cyl[:,:,0],z_cyl[:,:,0]), method='linear',fill_value=0,rescale=True)
            #PROP = np.repeat(PROP,4,axis=-1)
            PROP = np.repeat(np.expand_dims(PROP,axis=-1),len(self.phi),axis=-1)
        else:
            for vec,key in zip(prop_tot,['r','theta','phi']):
                #PROP[key] = interpolate.griddata((X,Y,Z_CYL), vec, (x,y,z), method='nearest',rescale=True)
                PROP[key] = interpolate.griddata((R_CYL,Z_CYL), vec, (r_cyl[:,:,0],z_cyl[:,:,0]), method='linear',fill_value=0,rescale=True)
                #PROP[key] = np.repeat(PROP[key],4,axis=-1)
                PROP[key] = np.repeat(np.expand_dims(PROP[key],axis=-1),len(self.phi),axis=-1)
        return PROP
    
    def isolate_stream(self, dphi_frac = 0.3, nstream = 1):
        shock = False
        np0 = len(self.phi)
        npstream = int(dphi_frac*np0)/nstream #how many phi0 per stream
        stream_phi = np.array([(np.arange(0,npstream,1) + i*np0/nstream).astype(int) for i in range(nstream)]).flatten()
        
        r_s = np.array([])
        th_s = np.array([])
        phi_s = np.array([])
        stream_s = np.array([])
        thstart = np.radians(self.env['theta_min'])
        acosj = np.arccos(np.linspace(np.cos(thstart), 1e-6, 180))
        
        dphi = np.gradient(self.coords[2])[0]
        subset_phi = self.phi[stream_phi]
        index_phi = stream_phi


        for j in acosj:
            streamline = self.stream(th0=j,p0=0,shock=shock)
            r_s = np.append(r_s,streamline['path'][0])
            th_s = np.append(th_s,streamline['path'][1])
            phi_s = np.append(phi_s,streamline['path'][2])
            stream_s = np.append(stream_s,np.ones_like(streamline['path'][0]))

        for ip in range(1,np0):
            if ip in index_phi:
                stream_ = np.ones_like(r_s)
            else:
                stream_ = np.zeros_like(r_s)
            stream_s = np.append(stream_s,stream_)

        r_s = np.repeat(r_s,np0,axis=-1).flatten()
        th_s = np.repeat(th_s,np0,axis=-1).flatten()
        phi_s = np.mod(np.array([phi_s + dphi for dphi in self.phi]).flatten(),2*pi)
        
        zcyl_s = r_s*np.cos(th_s)
        rcyl_s = r_s*np.sin(th_s)

        R_CYL,Z_CYL = self.make_rz()
        stream_mask = np.zeros_like(R_CYL[:,:,0])
        R = R_CYL[:,:,0]
        Z = Z_CYL[:,:,0]
        for iphi in self.phi:
            rcyl = rcyl_s[np.where(np.abs(phi_s-iphi)<=dphi/2.)]
            zcyl = zcyl_s[np.where(np.abs(phi_s-iphi)<=dphi/2.)]
            stream_val = stream_s[np.where(np.abs(phi_s-iphi)<=dphi/2.)]
            if len(rcyl) > 0:
                mask = np.clip(interpolate.griddata((rcyl,zcyl),stream_val,(R,Z),method='nearest',fill_value=0,rescale=True),a_min = 0.0, a_max=None)
            else:
                mask = np.zeros_like(R)
            stream_mask = np.dstack([stream_mask,mask])
        self.stream_mask = stream_mask[:,:,1:]

        
    def rho_env(self,fluid=0):
        rho0 = self.env['rho_0']
        rho_vol = self.solve_envelope(prop='rho')
        if fluid == 0:
            return rho_vol + self.env['rho_amb']
        elif fluid == 1:
            return (rho_vol + self.env['rho_amb'])*self.env['d2g']
        else:
            return rho_vol*0 
        
    def rho_streamvelope(self,fluid=0):
        #stream_mask = self.isolate_stream(dphi_frac=dphi_frac,nstream=nstream)
        rho0 = self.env['rho_0']
        rho_vol = self.solve_envelope(prop='rho')*self.stream_mask
        if fluid == 0:
            return rho_vol + self.env['rho_amb']
        elif fluid == 1:
            return (rho_vol + self.env['rho_amb'])*self.env['d2g']
        else:
            return rho_vol*0 
        
    def v_env(self):
        V = self.solve_envelope(prop='v',shock=False)
        return V
    
    def v_disk(self):
        R_CYL,Z_CYL = self.make_rz()
        R = np.sqrt(R_CYL**2 + Z_CYL**2)
        V = {}
        V['phi'] = self.vk(R)
        vin = -self.star['accrate']/(2*pi*R_CYL*np.sqrt(2*pi)*self.H(R_CYL)*AU*self.rho_disk(fluid=0))
        V['r'] = vin*R_CYL/R
        V['theta'] = vin*Z_CYL/R
        return V
    
    def v_embedded(self):
        shock = self.env['shock']
        if self.env['nstreams'] == 1 and self.env['stream_frac'] == 1:
            rho_e = self.rho_env(fluid=0)
        else:
            rho_e = self.rho_streamvelope(fluid=0)
        if shock == True:
            V = self.solve_envelope(prop='v')
        if shock == False:
            V = {}
            for key in ['r','theta','phi']:
                V[key] = (self.v_disk()[key]*self.rho_disk() + self.v_env()[key]*rho_e)/self.rho_embedded()
        return V
            
    def rho_embedded(self,fluid=0,rhod=''):
        shock = self.env['shock']
        rho_d = self.rho_disk(fluid=fluid) #calc disk from model
        if self.env['nstreams'] == 1 and self.env['stream_frac'] == 1:
            rho_e = self.rho_env(fluid=fluid)
        else:
            rho_e = self.rho_streamvelope(fluid=fluid)
        floor = self.env['rho_amb']*np.array([1,self.env['d2g'],0])
        rho_both = np.clip(rho_d + rho_e,a_min=floor[fluid],a_max=None)
        if fluid != 0: #no dust inside sublimation radius
            R,THETA,PHI = self.make_grid()
            rho_both[self.T(R) > 1500.] = floor[fluid]
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

    
    

    
    
    
    
    
    