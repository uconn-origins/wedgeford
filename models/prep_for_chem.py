#Converts radmc3d model input/output into format needed for UV/X-ray radiative transfer and chemistry
dir = './example_gmaur/'
outdir = '../highenergy_RT/'
import sys
sys.path.append(dir)
import final_setup as p
import numpy as np
import os

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
c = 2.99792458e10 # speed of light
mu = 2.36
h = 6.6260755e-27 #planck
mh = 1.67e-24 # mass of hydrogen cgs
Habing = 1.83590e+08  # Photons/cm^2/s between 930-2000
G0 = 2.7e-3 # erg/cm^2/s^-1


au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]


def write_out(model,outname='model',ndust=2): #write the .out file for RT and chemistry, pass the model class into this
    m = model
    if os.getcwd() != m.outdir:
        os.chdir(m.outdir)
    r,z = m.make_rz_H()
    T = m.get_Tdust()
    rho_g = m.rho_embedded(fluid=0)
    rho_d = np.zeros_like(T)
    for fluid in np.arange(ndust):
        rho_d += m.get_rhodust(fluid=fluid+1)
    rhod_2d = m.make_quadrant(rho_d, fill_value = np.median(rho_d)/2.,order='F').T.flatten()
    T_2d = m.make_quadrant(T,fill_value = np.median(T)/2.,order='F').T.flatten()
    rhog_2d = m.make_quadrant(rho_g, fill_value = np.median(rho_g)/2., order='C').T.flatten()
    
    header1 = 'R(AU)    z(AU)   rhogas(g/cm3)  rhodust(g/cm^3)    T(K)   [    4 types of dust]      f_H\n'
    header2 = '-ForChem+Plots--------------------------------------------------------\n'
    outfile = open(outname+'.out',"w")
    outfile.write(header1)
    outfile.write(header2)
    for n in np.arange(len(r.flatten())):
        outfile.write('%10.4f %10.4f   %5.2e   %5.2e %7.1f %12.7f %12.7f %12.7f %12.7f   %7.2e\n'%(r.T.flatten()[n],z.T.flatten()[n],rhog_2d[n], rhod_2d[n] ,T_2d[n],1.0,0.0,0.0,0.0,1.0))
    outfile.close()
    return outname + '.out'

def read_out(outname='model.out'): #convenience function for reading the .out file and returning the bits accessible as dictionary
    data = np.loadtxt(outname,skiprows=2,usecols=(0,1,2,3,4))
    r = data[:,0]
    z = data[:,1]
    shape = (len(np.unique(r)),50)
    rhog = data[:,2].reshape(shape)
    rhod = data[:,3].reshape(shape)
    T = data[:,4].reshape(shape)
    out_dict = {'r':r.reshape(shape), 'z':z.reshape(shape),'rhog':rhog,'rhod':rhod,'Td':T}
    return out_dict

def plot_out(prop='rhog',outname='model.out',log=True,method=contourf,**pk): #plotting convenience function to inspect .out file contents
    data = read_out(outname)
    if log == True:
        method(data['r'],data['z'],np.log10(data[prop]),**pk)
    else:
        method(data['r'],data['z'],(data[prop]),**pk) 




def main():
	dust2gas = 0.01
	rho_si = 3.0 #cancels out
	amin_chem = 0.06
	amax_ism = 1.0
	ndust = 3
	zetacr = 1.3e-17
	amin = [0.005e-6,0.005e-6,0.005e-6]#meters
	amax = [1e-6,1e-6,1e-3]#meters

	infile = open(dir+'amr_grid.inp',"r")
	n = 6
	line = infile.readlines()
	infile.close()
	grid_sty = int(line[1]) #0 is regular, 1 is oct tree amr, 10 is layer
	print('grid_sty ',grid_sty)
	coord_sys = int(line[2]) #cartesian<100, 100<=spherical<200, 200<=cylindrical<300
	print("coord_sys", coord_sys)
	nx = int(line[5].split(' ')[0])
	ny = int(line[5].split(' ')[1])
	nz = int(line[5].split(' ')[2])
	print('nx ',nx,'ny ',ny,'nz ',nz)

	infile = open(dir+'dust_density.inp','r')
	line = infile.readlines()
	infile.close()
	ndust = int(line[2])


	if grid_sty == 0:
		amr = np.loadtxt(dir+'amr_grid.inp',skiprows=6)
		rhod = np.loadtxt(dir+'dust_density.inp',skiprows=3)# ‘F’ means to index the elements in column-major, Fortran-style order, with the first index changing fastest, and the last index changing slowest
		rhod = np.reshape(rhod,(ndust,int(len(rhod)/ndust)),order='C')
		Td = np.loadtxt(dir+'dust_temperature.dat',skiprows=3)
		Td = np.reshape(Td,(ndust,int(len(Td)/ndust)),order='C')
		print("file shapes",amr.shape,rhod.shape,Td.shape)
		if coord_sys >= 100 and coord_sys < 200:
            		rr = p.rr
            		tt = p.tt
            		pp = rr*np.sin(tt)
            		zz = rr*np.cos(tt)
            		pp = pp.ravel(order='F')/au
            		zz = zz.ravel(order='F')/au
            		print(np.amin(pp),np.amax(pp))
            		print(np.amin(zz),np.amax(zz))
            		print("ravled arrays",pp.shape,zz.shape,rhod.shape)
		else:
            		sys.exit("need to write function for this coordinate system")

	elif grid_sty == 1:
        	lvlmax = int(line[6].split(' ')[0])
        	nleafmax = int(line[6].split(' ')[1])
        	nbranchmax = int(line[6].split(' ')[2])
        	sys.exit("need to write function for this file type")
	else:
		sys.exit("need to write function for this file type")

	try:
		rhog = p.rhog.ravel(order='F')
	except:
		rhog = np.sum(rhod,axis=0)/dust2gas #this needs fixed for more than 1 dust type

	#make grid with range of z for each r
	rin = p.rin/au
	rout = p.rout/au
	thetaup = p.thetaup
	nr = int(nx/2)
	#rhofac = p.rhofac
	#rvec = np.logspace(np.log10(rin),np.log10(rout/rhofac),nr)
	rvec = np.logspace(np.log10(rin),np.log10(rout),nr)
	#rvec = 0.5*(rvec[0:nr]+rvec[1:nr+1]) #center of cell
	rarr = np.zeros((nr,ny))
	zarr = np.zeros((nr,ny))

 	#flip along midplane for interpolation
	pp = np.concatenate((pp,np.flip(pp,axis=-1)))
	zz = np.concatenate((zz,np.flip(-1*zz,axis=-1)))
	rhod = np.concatenate((rhod,np.flip(rhod,axis=-1)),axis=-1)
	rhog = np.concatenate((rhog,np.flip(rhog,axis=-1)))
	Td = np.concatenate((Td,np.flip(Td,axis=-1)),axis = -1)
	for x in range(nr):
        	rarr[x,:] = rvec[x]
        	#zmax = rvec[x]
        	zmax = rvec[x]*np.cos(thetaup)
        	zmin = 0.01*rvec[x]
		#zvec = np.logspace(np.log10(zmax),np.log10(zmin),ny)
        	zvec = np.linspace(zmax,zmin,ny-1)
        	zvec = np.concatenate((zvec,[0.0]))
        	#zvec = 0.5*(zvec[0:ny]+zvec[1:ny+1]) #center of cell
        	zarr[x,:] = zvec
	#now interpolate to find values in new r z grid
	rhog_grid = griddata((pp,zz),rhog,(rarr,zarr),method='linear')
	rhod_grid = np.zeros((ndust,len(rhog_grid[:,0]),len(rhog_grid[0,:])))
	Td_grid = np.zeros_like(rhod_grid)
	for d in range(ndust):
		rhod_grid[d,:,:] = griddata((pp,zz),rhod[d,:],(rarr,zarr),method='linear')
		Td_grid[d,:,:] = griddata((pp,zz),Td[d,:],(rarr,zarr),method='linear')
	rhod_grid[np.where(rhod_grid < 1e-99)] = 0.0
	rhog_grid[np.where(rhog_grid < 1e-99)] = 0.0

	r = rarr.ravel(order='C')
	z = zarr.ravel(order='C')
	rhog = rhog_grid.ravel(order='C')
	rhod = np.zeros((ndust,len(rhog)))
	Td = np.zeros_like(rhod)
	for d in range(ndust):
		rhod[d,:] = rhod_grid[d,:,:].ravel(order='C')
		Td[d,:] = Td_grid[d,:,:].ravel(order='C')


# Writing 1environ files

	if not os.path.exists(outdir):
		os.makedirs(outdir)

	for n in range(len(r)):
		#outfile = open(outdir+'1environ.inp.e1.'+np.round(r[n],decimals=4))
		#outfile.write('  R             rho           Tgas          Tdust        zAU           zcm           Nrz           ZetaCR        DustFrac\n')
		#outfile.write('1\n')
		#outfile.write(str(nz)+'\n')
		if n !=0 and r[n] == r[n-1]:
			nsig_unset = -3.0 *rhog[n] * dust2gas/rho_si * (amax_ism**-0.5 - amin_chem**-0.5)/ \
			(amax_ism**-0.5 - amin_chem**-0.5)
			nsig_tot = 0.0
			for d in range(ndust):
				nsig = -3.0 *rhod[d,n]/rho_si * (amax[d]**-0.5 - amin[d]**-0.5)/ \
				(amax[d]**-0.5 - amin[d]**-0.5)
				nsig_tot = nsig_tot + nsig
			dustfrac = nsig_tot/nsig_unset
			zcm = zcm + np.abs(z[n]-z[n-1])
		else:
			zcm = 1e-10
			print("new radius:",r[n],z[n])
		#outfile.write(' %6.2E   %6.2E   %6.2E   %6.2E   %6.2E   %6.2E   %6.2E   %6.2E   %6.2E\n'%(r[n],rhod[n],Td[n],Td[n],z[n],zcm,1e10,zetacr,dustfrac))

main()





