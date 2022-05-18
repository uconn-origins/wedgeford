#Converts radmc3d model input/output into format needed for UV/X-ray radiative transfer and chemistry


dir = '/home/kschwarz/radmc-3d/version_0.41/examples/gmaur/'
outdir = '../highenergy_RT/'
import sys
sys.path.append(dir)
import final_setup as p
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
import os

au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]

def main():
	dust2gas = 0.01
	rho_si = 3.0 #cancels out
	amin_chem = 0.06
	amax_ism = 1.0
	ndust = 3
	zetacr = 1.3e-17
	amin = [0.005e-6,0.005e-6,0.005e-6]
	amax = [1e-6,1e-6,1e-3]

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
	nz =int( len(z)/len(np.unique(r)))

# Writing .out file
	#outfile = open('gas_radmctest.out',"w")
	#outfile.write('R(AU)    z(AU)   rhogas(g/cm3)  rhodust(g/cm^3)    T(K)   [    4 types of dust]      f_H\n')
	#outfile.write('-ForChem+Plots--------------------------------------------------------\n')
	#    outfile.write('%10.4f %10.4f   %5.2e   %5.2e %7.1f %12.7f %12.7f %12.7f %12.7f   %7.2e\n'%(r[n],z[n],rhog[n],np.sum(rhod[:,n]),Td[n],1.0,0.0,0.0,0.0,1.0))
    #outfile.close()

#zcm is height from surface where surface is 1e-10cm
#Nrz is radial column number density, deprecated 

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





