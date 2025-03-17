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
C ngr(Nz)     == dust number density [1/cm^3],
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
C
C nxray_photons(Nz) == # x-ray photons, from Tom Bethell scattering calculation [photons]
C
C xrayrate(Nz) == xray ionization rate at each zone
C..............................................................................
C Global parameter(s):
      INTEGER nzone, nwavl, nheight, ZMAX, RMAX
      PARAMETER (nzone=200, nwavl=120, nheight=100, ZMAX=4000,
     1			 RMAX=4000)

C Common block 2 (input parameters)
	  INTEGER Nr, Nz, zone, UVmaxzone
      INTEGER isavzone
      DOUBLE PRECISION Rs, Tg, Td, ZetaCR, rho, ngr, Nrz
	  DOUBLE PRECISION uvfield, xrayfield, lambda
	  DOUBLE PRECISION zAU, zcm, CRatten, totflux, ndust
	  DOUBLE PRECISION xraylevels, xrayrate, xrayratesimon
	  DOUBLE PRECISION isrffield, RNatten, RNrate,locdust
      DOUBLE PRECISION rsavzone
	  CHARACTER*10 dust, tmpdust
	  LOGICAL firstzone,xraydust,incl_radionuc,incl_isrf,incl_2dabun
	  LOGICAL write_2dabun, incl_locdust
      LOGICAL runzone
	  DIMENSION isavzone(nzone,100), Rs(nzone), Tg(nzone), Td(nzone), ZetaCR(nzone),
     1 rho(nzone), ngr(nzone), Nrz(nzone), uvfield(nheight,nwavl),
     2 xrayfield(nheight,nwavl), lambda(nwavl), zAU(nzone),
     3 zcm(nzone), CRatten(nzone), totflux(nzone),
     4 UVmaxzone(nwavl), xraylevels(nwavl), xrayrate(nzone),
     5 xrayratesimon(nzone), isrffield(nheight,nwavl),
     6 RNatten(nzone,2), RNrate(nzone), locdust(nzone), rsavzone(nzone,100), runzone(nzone)

	  COMMON /BL2/ uvfield, xrayfield, Nr, Nz, UVmaxzone, isavzone, Rs, Tg,
     1  Td, ZetaCR, rho, ngr, zAU, zcm, lambda, CRatten, totflux,
     2  xraylevels, xrayrate, ndust, zone, firstzone, xraydust,
     3  dust, tmpdust, xrayratesimon, isrffield, RNatten, RNrate, rsavzone,
     4  incl_radionuc,incl_isrf,incl_2dabun, runzone
