C..............................................................................
C
C Common Variables:
C
C..............................................................................
C
C Tcr		== CR temperature (CR desorption)
C
C frac		== time spent at 70 K (CR desorption)
C
C g0_100	== G0 value, normalized at 100 AU
C
C g0_CO 	== G0 value for CO (ignore lyman alpha)
C
C Lxray		== X-ray luminosity of the central star (erg/s)
C
C RadNuc	  == Radionuclide Ionization Rate [1/s]
C
C fg        == dust to gas mass ratio (normal is 100)
C
C CRdesorp	== turn on/off CR desorption
C
C CRionization	== turn on/off CR ionization
C
C photodesorp == turn on/off photodesorption
C
C include_LyA == remove lya from UV field or not
C..............................................................................
	INTEGER NBANDS
	PARAMETER (NBANDS=150)

	DOUBLE PRECISION Tcr, frac
	PARAMETER(Tcr=70.0, frac=3.16D-19)
	DOUBLE PRECISION g0_100, g0_CO, Lxray, fg, freezeeffic
	LOGICAL CRdesorp, CRionization, photodesorp, LyAphotodesorp
	LOGICAL thermaldesorp, include_LyA
	DOUBLE PRECISION RadNuc

	COMMON /FLAGS/ g0_100, g0_CO, Lxray, fg, freezeeffic, CRdesorp,
     &		CRionization, photodesorp, LyAphotodesorp, thermaldesorp,
     &		include_LyA, RadNuc

C DOUBLE PRECISION xrayrate
C
C COMMON /FLAGS/ g0_100, g0_CO, Lxray, fg, xrayrate, CRdesorp,
C
C..............................................................................
C
C XSECT:
C
C siglam_l     == lower lambda value from xsect file
C
C siglam_u     == upper lambda value from xsect file
C
C sigma        == cross section from xsect file
C
C..............................................................................
	INTEGER nband_read, nlam, xnlam
	DOUBLE PRECISION siglam_l, siglam_u, sigma, e_sil, xs_sil
	DIMENSION siglam_l(NBANDS), siglam_u(NBANDS), sigma(NBANDS),
     &		e_sil(2000), xs_sil(2000)

	COMMON /XSECT/ siglam_l, siglam_u, sigma, e_sil, xs_sil,
     &		nband_read, nlam, xnlam

C..............................................................................
C
C SELFSHIELD:
C
C ih2		==
C
C ico		==
C
C..............................................................................
	INTEGER ih2, ico
	PARAMETER (ih2=6, ico=8)
	CHARACTER*80 fnameSS

	COMMON /SSTEST/ fnameSS

C..............................................................................
C
C Monolayers
C
C n_ice_init == Total abundance of species on grains initially
C
C n_ice	== Total abundance of species on grains
C
C numlayers	== number of layers on grain
C
C PDadjust	== monolayer adjustment for photodesorption
C
C Madjust	== monolayer adjustment for other reactions
C
C..............................................................................
	DOUBLE PRECISION ngr_init, n_ice_init, n_ice, numlayers, PDadjust, 
     1  Madjust

	COMMON /MLAYER/ ngr_init, n_ice_init, n_ice, numlayers, PDadjust,
     1  Madjust
