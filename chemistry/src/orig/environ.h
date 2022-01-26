C..............................................................................
C
C Global parameter(s):
C
C nzone       == maximal number of zones
C 
C..............................................................................
C
C Common block(s):
C
C..............................................................................
C
C BL2:
C
C Input Parameters (from 1environ.inp)
C
C Nr          == number of radial grid points,
C
C Nz          == number of vertical grid points,
C
C Rs(Nz)      == radius (AU),
C
C Tg(Nz)       == gas temperature [K],
C
C Td(Nz)       == dust temperature [K],
C
C rho(Nz)     == density [g/cm^3],
C
C ZetaCR(Nz)  == cosmic ray ionization rate.
C
C Nrz(Nz)	  == Radial column density [1/cm^2], used for x-ray ionziation rate
C
C zAU(Nz)	  == Height above the midplane [AU]
C
C zcm(Nz)	  == Height from disk surface, surface = 1.0e10 cm [cm]
C
C zone		  == variable that holds the number of the current zone
C
C CRatten(Nz) == Cosmic Ray attenuation factor
C
C totflux(Nz) == integrated flux [photons / cm^2 / s]
C..............................................................................
C Global parameter(s):
      INTEGER nzone, nwavl, nheight, ZMAX, RMAX
      PARAMETER (nzone=200, nwavl=120, nheight=100, ZMAX=4000, RMAX=4000)

C Common block 2 (input parameters)
	  INTEGER Nr, Nz, zone, UVmaxzone
      DOUBLE PRECISION Rs, Tg, Td, ZetaCR, rho, Nrz
	  DOUBLE PRECISION uvfield, lambda
	  DOUBLE PRECISION zAU, zcm, CRatten, totflux, ndust
	  CHARACTER*10 dust, tmpdust
	  LOGICAL firstzone, xraydust
	  DIMENSION Rs(nzone), Tg(nzone), Td(nzone), ZetaCR(nzone), 
     1 rho(nzone), Nrz(nzone), uvfield(nheight,nwavl), 
     2 lambda(nwavl), zAU(nzone), zcm(nzone), CRatten(nzone),
     3 totflux(nzone), UVmaxzone(nwavl)
	  
	  COMMON /BL2/ uvfield, Nr, Nz, UVmaxzone, Rs, Tg, Td, ZetaCR, 
     1  rho, Nrz, zAU, zcm, lambda, CRatten, totflux, ndust, zone, 
     2  firstzone, xraydust, dust, tmpdust
	 
