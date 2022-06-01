C..............................................................................
C
C This function checks to see if the on grain abundance is larger than
C the grain abundance and, if so, reduces the rate of the reaction
C
C..............................................................................
	DOUBLE PRECISION FUNCTION grain_abun_check(ralpha,rbeta,rgamma,species)
	implicit none

C Initialization of common blocks:
	INCLUDE "Fcn.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	DOUBLE PRECISION ralpha, rbeta, rgamma
	DOUBLE PRECISION nr1, ngr
	DOUBLE PRECISION calcrate
	DOUBLE PRECISION vel
	DOUBLE PRECISION EHc,EHp,Es,nuHc,Flu,squiggle,effh2
	CHARACTER*13 species
	INTEGER r1index, grindex
        DOUBLE PRECISION sigadjust

        sigadjust = 1.0
        if (incl_locdust) sigadjust = locdust(zone)

c	calcrate = ralpha*(Tg(zone)/300.0D0)**rbeta
	vel = ((8.0 * kbol * Tg(zone))/(pi * rbeta * mH))**0.5

C   Now adding efficiency term from Cazaux and Tielens 2002, 2004

	EHc  = 1e4
	EHp  = 600.0
	Es   = 200.0
	nuHc = 1.3e13
	Flu = 1e-10

C	squiggle = (1+nuHc/(2.0*Flu)*exp(-1.5*EHc/(Tg(zone)))
C     &		 *(1+sqrt((EHc - Es)/(EHp-Es)))**2)**(-1)


C	effh2 = (1.0+0.25*(1+sqrt((EHc - Es)/(EHp-Es)))**2
C     &		* exp(-Es/(Tg(zone))))**(-1)*squiggle




	effh2 = 1.0

	calcrate = ralpha * vel * rgamma * sigadjust * effh2

	CALL ispecies(species, ns, s, r1index)
	CALL ispecies('GRAIN        ', ns, s, grindex)
	IF (timestep .EQ. 1) THEN
		nr1 = abundances(r1index, timestep)
		ngr = abundances(grindex, timestep)
	ELSE
		nr1 = abundances(r1index, timestep-1)
		ngr = abundances(grindex, timestep-1)
	END IF

	IF (ngr .LT. 1.0e-20) THEN
		ngr = ngr_init
	ENDIF

	IF (nr1 .GT. ngr) THEN
		calcrate = calcrate * ngr / nr1
	ENDIF

	if (nr1.lt.1e-30) calcrate = 0.0
	if (calcrate .lt. 1e-30) calcrate = 0.0

	grain_abun_check = calcrate
c	print *, 'H freeze: ', calcrate, 7.5e-5*(Tg(zone)/300.0)**0.5

	RETURN
	END

C..............................................................................
C
C binding_surface
C Binding energies should be all reflective of bare silicate (H2O and CO included).
C 855.0*1180.0/960.0
C
C..............................................................................
	DOUBLE PRECISION FUNCTION binding_surface(Ebind,massSp,rgamma)
	implicit none

C Initialization of common blocks:
	INCLUDE "Fcn.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	INTEGER i_co,i_co2,i_h2o
	DOUBLE PRECISION nh2o,nco,nco2
	DOUBLE PRECISION f_h2o,f_co,f_co2
	DOUBLE PRECISION nu0,Ebind,massSp,rgamma,icetot
	DOUBLE PRECISION calcrate
	DOUBLE PRECISION h2o_v_sil, co2_v_sil, co_v_sil
	DOUBLE PRECISION Ebind_adj,nHtot, coicetot

	CALL ispecies('CO(gr)       ', ns, s, i_co)
	CALL ispecies('H2O(gr)      ', ns, s, i_h2o)
	CALL ispecies('CO2(gr)      ', ns, s, i_co2)
	IF (timestep .EQ. 1) THEN
		nco  = abundances(i_co,  timestep)
		nco2 = abundances(i_co2, timestep)
		nh2o = abundances(i_h2o, timestep)
	ELSE
		nco  = abundances(i_co,  timestep-1)
		nco2 = abundances(i_co2, timestep-1)
		nh2o = abundances(i_h2o, timestep-1)
	END IF

	icetot = nco + nco2 + nh2o
	coicetot = nco + nco2
	nHtot = rho(zone) / aMp / amu

	f_co = nco/coicetot
	f_h2o = nh2o/coicetot
	f_co2 = nco2/coicetot

C Evaluate fractional coverage of ice, as well as total amount of ice coveraged.  For example, 1 monolayer of ice is 1e6 sites.
C The abundance of a single monolayer of ice is: icetot/6e-12 = 1e6.  Or, icetot = 6e-6 is the minimum total
C abundance of ice to have a monolayer. If the abundance is below x(ice) < 1e-6 then ice is not fully coated and
C bare silicate binding energies are assumed.

	h2o_v_sil = 1323.0/857.0 * 855.0
	co2_v_sil = 1113.0/857.0 * 855.0
	co_v_sil  = 855.0


c	if (icetot/nHtot.ge.6e-6) then
c		IF (Td(zone) .le. 17.1) THEN      ! Below 17 K the ice mantles will be dominated by CO
c			Ebind_adj = Ebind*co_v_sil
c		ELSE IF (nco/nHtot.ge.6e-6) THEN  ! If you have a monolayer of CO alone:
c			Ebind_adj = Ebind*co_v_sil
c		ELSE
c			Ebind_adj = Ebind*(f_h2o*h2o_v_sil + f_co2*co2_v_sil +
c     &         f_co*co_v_sil)
cc		ENDIF
c	endif

C   17.1 for CO
C	22.16 for CO2
C   26.46 for H2O
C   96 K for H2O in gas (silicate surf)

C  If dust > 26.46

C CO2 freezes onto H2O at XX Temp.  Between 17.1 - CO2 on H2O
C For Td < 17.1, CO surface 855 K.
C For Td > 17.1 and < Eb(CO2/H2O) = 48 K, have a CO2 surface or 1113.0 K.
C For Td > 48 K and < Eb(H2O/SiO2) = 96 K, have an H2O surface or 1323.0 K.
C For Td > Eb(H2O/SiO2) = 96 K, have a SiO2 surface surface or 1180.0 K.

	Ebind_adj = 1180.0 ! For bare silicate.

	if (icetot/nHtot.ge.6e-6) then
		IF (Td(zone) .le. 25.0) THEN      ! Below 17 K the ice mantles will be dominated by CO
			Ebind_adj = co_v_sil
		ELSE IF ((Td(zone) .gt. 25.0).and.(Td(zone) .le. 50.0)) THEN  ! If you have a monolayer of CO alone:
c		ELSE IF (Td(zone) .gt. 30.0) THEN  ! If you have a monolayer of CO alone:
			Ebind_adj = co2_v_sil
		ELSE IF ((Td(zone) .gt. 50.0).and.(Td(zone) .le. 100.0)) THEN  ! If you have a monolayer of CO alone:
			Ebind_adj = h2o_v_sil
		ELSE
			Ebind_adj = 1180.0
		ENDIF
	endif
c	if (icetot/nHtot.ge.6e-6) then
c		Ebind_adj = co_v_sil
c	else
c		Ebind_adj = 1180.0
c	endif

c	if (icetot/nHtot.ge.6e-6) then
c		IF (Td(zone) .le. 17.1) THEN      ! Below 17 K the ice mantles will be dominated by CO
c			Ebind_adj = Ebind*co_v_sil
c        ELSE IF ((Td(zone) .gt. 48.0).and.(Td(zone) .le. 96.0)) THEN
c			Ebind_adj = Ebind*h2o_v_sil
c		ELSE
c			Ebind_adj = Ebind*(f_co2*co2_v_sil + f_co*co_v_sil)
C		ELSE IF ((Td(zone) .gt. 48.0).and.(Td(zone) .le. 96.0)) THEN  ! If you have a monolayer of CO alone:
C			Ebind_adj = Ebind*h2o_v_sil
C		ELSE
C			Ebind_adj = Ebind
c		ENDIF
c	endif

c	print *,
c	write(*,'(a16,F10.3,a15,F10.3,a15,E10.3,a15,E10.3,a15,f10.3,f10.3,
c     &		 f10.3)')
c     &		 'Ebind silicate: ',Ebind,' and adjusted: ',Ebind_adj,
c     &		 ' Ice abundance:',icetot/nHtot,
c     &       ' CO ice abund:',nco/nHtot,'fractions:',f_co,f_co2
c	write(*, '()') calcrate, calcrate * ngr / nr1


	nu0 = ((2.0 * rgamma * kbol * Ebind_adj)/(pi*pi*mH*massSp))**0.5
	calcrate = nu0 * dexp(-Ebind_adj / Td(zone))
	if (calcrate .lt. 1e-30) calcrate = 0.0
	binding_surface = calcrate

	RETURN
	END


C..............................................................................
C
C Cazaux and Tielens grain surf H2 formation
C
C..............................................................................
	DOUBLE PRECISION FUNCTION h2form(ralpha,rbeta,rgamma,species)
	implicit none

C Initialization of common blocks:
	INCLUDE "Fcn.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	DOUBLE PRECISION ralpha, rbeta, rgamma
	DOUBLE PRECISION nr1, ngr
	DOUBLE PRECISION calcrate
	DOUBLE PRECISION vel
	DOUBLE PRECISION EHc,EHp,Es,nuHc,Flu,squiggle,effh2
	CHARACTER*13 species
	INTEGER r1index, grindex
        DOUBLE PRECISION sigadjust

        sigadjust = 1.0
        if (incl_locdust) sigadjust = locdust(zone)

c	calcrate = ralpha*(Tg(zone)/300.0D0)**rbeta
	vel = ((8.0 * kbol * Tg(zone))/(pi * rbeta * mH))**0.5


C   Now adding efficiency term from Cazaux and Tielens 2002, 2004

	EHc  = 1e4
	EHp  = 600.0
	Es   = 200.0
	nuHc = 1.3e13
	Flu = 1e-10

	squiggle = (1+nuHc/(2.0*Flu)*exp(-1.5*EHc/(Tg(zone))) *
     &		 (1+sqrt((EHc - Es)/(EHp-Es)))**2)**(-1)


	effh2 = (1.0+0.25*(1+sqrt((EHc - Es)/(EHp-Es)))**2
     &		 * exp(-Es/(Tg(zone))))**(-1)*squiggle
c	print *, 'H2 form',Tg(zone),effh2,squiggle
c      print *, ralpha,rbeta,rgamma,species

	calcrate = ralpha * vel * rgamma * sigadjust * effh2

	CALL ispecies(species, ns, s, r1index)
	CALL ispecies('GRAIN        ', ns, s, grindex)
	IF (timestep .EQ. 1) THEN
		nr1 = abundances(r1index, timestep)
		ngr = abundances(grindex, timestep)
	ELSE
		nr1 = abundances(r1index, timestep-1)
		ngr = abundances(grindex, timestep-1)
	END IF

	IF (ngr .LT. 1.0e-20) THEN
		ngr = ngr_init
	ENDIF

	calcrate = calcrate * ngr / nr1

	if (nr1.lt.1e-20) calcrate = 0.0

c	print *,'calcrate h2form:',calcrate,species

	h2form = calcrate

	RETURN
	END

C..............................................................................
	DOUBLE PRECISION FUNCTION compmono()
	IMPLICIT NONE

C Include the environmental variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"
	INCLUDE "Fcn.h"
	INCLUDE "constants.h"

	CHARACTER*13 specname, nogrspec
	INTEGER i, grspecindex, grindex,Hindex,Dindex
	DOUBLE PRECISION  Nsites, ngrspec
	DOUBLE PRECISION ngr, Mlayers
	DOUBLE PRECISION fr,molarea

C This part of the code currently assumes all grains are small but there are fewer of them.

	Nsites = 1.0e6

	n_ice = 0.0
	CALL ispecies('H(gr)        ', ns, s, Hindex)
	CALL ispecies('D(gr)        ', ns, s, Dindex)

C Calculate the abundance of on-grain species
	i = 0
	do while (grspecs(i) .NE. -1)
		IF (timestep .EQ. 1) THEN
            if ((grspecs(i) .ne. Hindex).and.(grspecs(i) .ne. Dindex))
     &		 n_ice = n_ice + abundances(grspecs(i), timestep)
		ELSE
            if ((grspecs(i) .ne. Hindex).and.(grspecs(i) .ne. Dindex))
     &		 n_ice = n_ice + abundances(grspecs(i), timestep-1)
		ENDIF
		i = i+1
	end do

c Get the current abundance of grains
	CALL ispecies('GRAIN        ', ns, s, grindex)
	IF (timestep .EQ. 1) THEN
		ngr = abundances(grindex, timestep)
	ELSE
		ngr = abundances(grindex, timestep-1)
	END IF


c set the abundance for the first zone (when abundances = 0)
	IF (n_ice .LT. MINABUN) THEN
		n_ice = n_ice_init
	ENDIF
	IF (ngr .LT. MINABUN) THEN
		ngr = ngr_init
	ENDIF

c Number of monolayers is nice / ngr / Nsites
c Abundance of monolayers is then nice / Mlayers = ngr * Nsites
	compmono = n_ice / ngr / Nsites

	RETURN
	END


C..............................................................................
	DOUBLE PRECISION FUNCTION occupyice(species)
	IMPLICIT NONE

C Include the environmental variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"
	INCLUDE "Fcn.h"
	INCLUDE "constants.h"

	CHARACTER*13 specname, nogrspec, species
	INTEGER i, grspecindex, grindex,r1index
	DOUBLE PRECISION  ngrspec
	DOUBLE PRECISION ngr, Mlayers
	DOUBLE PRECISION fr,nr1

C This part of the code currently assumes all grains are small but there are fewer of them.

	nr1 = 0.0


c Get the current abundance of grains
	CALL ispecies('GRAIN        ', ns, s, grindex)
	CALL ispecies(species, ns, s, r1index)

	IF (timestep .EQ. 1) THEN
		ngr = abundances(grindex, timestep)
        nr1 = abundances(r1index, timestep)
	ELSE
		ngr = abundances(grindex, timestep-1)
        nr1 = abundances(r1index, timestep-1)
	END IF

	occupyice = nr1/ngr

	if (nr1 .lt. 1e-30) occupyice = 0.0


	RETURN
	END



C..............................................................................
C
C This function calculates the H2 photodissociation rate and takes into account
C the H2 self-shielding.
C Uses rates and k_Av0 from Lee et al. 1996, A&A, 311, 690
C..............................................................................
	DOUBLE PRECISION FUNCTION self_shield_H2()
	implicit none

C Initialization of common blocks:
	INCLUDE "Fcn.h"
	INCLUDE "rates.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	DOUBLE PRECISION avshock, k_Av0, sigmad_1000, taud_1000
	DOUBLE PRECISION colden, colden_h2, colden_h2_tot(ntime)
	DOUBLE PRECISION H2abun, maxlam, etau
	INTEGER n_lee, i
	PARAMETER (n_lee = 105)
	DOUBLE PRECISION lee_nh2(n_lee), lee_theta(n_lee), theta
	DOUBLE PRECISION lnh2(n_lee), interp
	DOUBLE PRECISION itheta

	double PRECISION G0eval

	DATA colden/0.0/, colden_h2/0.0/, colden_h2_tot/ntime * 0.0/
	DATA lee_nh2/n_lee * 0.0/, lee_theta/n_lee * 0.0/

	IF (zone .EQ. 1) then
		zcm(zone-1) = 0.0
	endif

	if (timestep .EQ. 1) then
		call lee_factors(lee_nh2, lee_theta)
	endif
c	write(*,*) '0. ColdenH2',colden_h2_tot(timestep)
	H2abun = rho(zone) / mH / amu
	colden_h2 = H2abun * (zcm(zone) - zcm(zone-1))
	colden_h2_tot(timestep) = colden_h2_tot(timestep) + colden_h2

c take log10 of column density table to prepare for interpolation
	do i=1,n_lee
		if (lee_nh2(i) .EQ. 0.0) then
			lnh2(i) = 0.0
		else
			lnh2(i) = dlog10(lee_nh2(i))
		endif
	enddo

c Linear interpolation
	theta = interp(dlog10(colden_h2_tot(timestep)),
     &		 lnh2, lee_theta, n_lee)

C Old method, calculating G0 on the fly now, LIC 9/18/2013
c Unshielded value
	k_Av0 = 2.54e-11 * g0_100 * (100.0 / Rs(zone))**2.0

c changed to I/I0 at 1000 Ang
	maxlam = 1000.0

C	call find_etau(maxlam, etau, zone)
c	self_shield_H2 = k_Av0 * theta * etau

	call find_G0(maxlam, G0eval, zone, 1)
	self_shield_H2 = 2.54e-11 * G0eval * theta

c write values to selfshield.test
	IF (shieldtest) THEN
		OPEN(unit=02, file=fnameSS, status='old',
     &		 ACCESS='APPEND')
		write(2,110), Rs(zone), self_shield_H2, k_Av0, theta, etau,
     &		 colden_h2_tot(timestep), timestep
		CLOSE(02)
	END IF
 110	format(1x, 'H2 ', 1x, f6.2, 5(1x,e10.3), i4)

	RETURN
	END

C..............................................................................
C Returns
C..............................................................................
	DOUBLE PRECISION FUNCTION interp(xval, x, y, npts)
	DOUBLE PRECISION xval
	INTEGER npts,i,ipt
	DOUBLE PRECISION x(npts), y(npts)


c Return first element if less than the first element in the array
	if (xval .LE. x(1)) then
		interp = y(1)
		return
	endif

c If beyond the last point, return the last value (no extrapolation)
	if (xval .GE. x(npts)) then
		interp = y(npts)
		return
	endif

c find points to interpolate around
	do i=1,npts
		if (xval .LT. x(i)) then
			ipt = i-1
			exit
		endif
	enddo

c 1D Linear Interpolation
	interp = y(ipt) + (xval-x(ipt)) * (y(ipt+1)-y(ipt))
     &		 / (x(ipt+1)-x(ipt))

	return
	end

C..............................................................................
C This function takes a column density as input and returns an H2 shielding
C factor based on Lee et al. 1996
C..............................................................................
	SUBROUTINE lee_factors(nh2, theta)
	implicit none
	INTEGER n_lee
	PARAMETER (n_lee = 105)
	DOUBLE PRECISION nh2(n_lee), theta(n_lee)

	nh2(1)  = 0.000
	nh2(2)  = 3.690e11
	nh2(3)  = 3.715e12
	nh2(4)  = 3.948e13
	nh2(5)  = 1.233e14
	nh2(6)  = 2.536e14
	nh2(7)  = 4.342e14
	nh2(8)  = 6.653e14
	nh2(9)  = 6.689e14
	nh2(10) = 9.075e14
	nh2(11) = 1.234e15
	nh2(12) = 1.631e15
	nh2(13) = 2.105e15
	nh2(14) = 2.363e15
	nh2(15) = 2.899e15
	nh2(16) = 3.207e15
	nh2(17) = 3.848e15
	nh2(18) = 4.636e15
	nh2(19) = 5.547e15
	nh2(20) = 6.604e15
	nh2(21) = 7.855e15
	nh2(22) = 9.368e15
	nh2(23) = 1.122e16
	nh2(24) = 1.352e16
	nh2(25) = 1.643e16
	nh2(26) = 2.017e16
	nh2(27) = 2.515e16
	nh2(28) = 3.190e16
	nh2(29) = 4.128e16
	nh2(30) = 5.439e16
	nh2(31) = 7.315e16
	nh2(32) = 1.009e17
	nh2(33) = 1.432e17
	nh2(34) = 2.092e17
	nh2(35) = 3.123e17
	nh2(36) = 4.738e17

	nh2(37) = 5.388e17
	nh2(38) = 8.935e17
	nh2(39) = 1.381e18
	nh2(40) = 2.164e18
	nh2(41) = 3.330e18
	nh2(42) = 5.024e18
	nh2(43) = 7.404e18
	nh2(44) = 9.029e18
	nh2(45) = 1.316e19
	nh2(46) = 1.813e19
	nh2(47) = 2.453e19
	nh2(48) = 3.248e19
	nh2(49) = 4.216e19
	nh2(50) = 5.370e19
	nh2(51) = 6.722e19
	nh2(52) = 8.277e19
	nh2(53) = 9.894e19
	nh2(54) = 1.186e20
	nh2(55) = 1.404e20
	nh2(56) = 1.644e20
	nh2(57) = 1.908e20
	nh2(58) = 2.197e20
	nh2(59) = 2.510e20
	nh2(60) = 2.849e20
	nh2(61) = 3.214e20
	nh2(62) = 3.604e20
	nh2(63) = 4.019e20
	nh2(64) = 4.456e20
	nh2(65) = 4.915e20
	nh2(66) = 5.393e20
	nh2(67) = 5.886e20
	nh2(68) = 6.392e20
	nh2(69) = 6.909e20
	nh2(70) = 7.433e20
	nh2(71) = 7.965e20

	nh2(72) = 8.505e20
	nh2(73) = 9.056e20
	nh2(74) = 9.627e20
	nh2(75) = 1.011e21
	nh2(76) = 1.068e21
	nh2(77) = 1.125e21
	nh2(78) = 1.185e21
	nh2(79) = 1.250e21
	nh2(80) = 1.327e21
	nh2(81) = 1.428e21
	nh2(82) = 1.578e21
	nh2(83) = 1.851e21
	nh2(84) = 2.128e21
	nh2(85) = 2.298e21
	nh2(86) = 2.389e21
	nh2(87) = 2.459e21
	nh2(88) = 2.519e21
	nh2(89) = 2.571e21
	nh2(90) = 2.618e21
	nh2(91) = 2.707e21
	nh2(92) = 2.790e21
	nh2(93) = 2.887e21
	nh2(94) = 3.001e21
	nh2(95) = 3.139e21
	nh2(96) = 3.303e21
	nh2(97) = 3.497e21
	nh2(98) = 3.722e21
	nh2(99) = 3.983e21
	nh2(100) = 4.283e21
	nh2(101) = 4.644e21
	nh2(102) = 5.127e21
	nh2(103) = 5.945e21
	nh2(104) = 8.205e21
	nh2(105) = 1.015e22


	theta(1)  = 1.000
	theta(2)  = 9.983e-1
	theta(3)  = 9.853e-1
	theta(4)  = 8.761e-1
	theta(5)  = 7.199e-1
	theta(6)  = 5.728e-1
	theta(7)  = 4.455e-1
	theta(8)  = 3.431e-1
	theta(9)  = 3.418e-1
	theta(10) = 2.732e-1
	theta(11) = 2.110e-1
	theta(12) = 1.619e-1
	theta(13) = 1.236e-1
	theta(14) = 1.084e-1
	theta(15) = 8.447e-2
	theta(16) = 7.410e-2
	theta(17) = 5.774e-2
	theta(18) = 4.416e-2
	theta(19) = 3.390e-2
	theta(20) = 2.625e-2
	theta(21) = 2.048e-2
	theta(22) = 1.606e-2
	theta(23) = 1.264e-2
	theta(24) = 9.987e-3
	theta(25) = 7.937e-3
	theta(26) = 6.343e-3
	theta(27) = 5.088e-3
	theta(28) = 4.089e-3
	theta(29) = 3.283e-3
	theta(30) = 2.640e-3
	theta(31) = 2.130e-3
	theta(32) = 1.725e-3
	theta(33) = 1.397e-3
	theta(34) = 1.129e-3
	theta(35) = 9.097e-4
	theta(36) = 7.340e-4

	theta(37) = 6.883e-4
	theta(38) = 5.377e-4
	theta(39) = 4.352e-4
	theta(40) = 3.475e-4
	theta(41) = 2.771e-4
	theta(42) = 2.205e-4
	theta(43) = 1.753e-4
	theta(44) = 1.549e-4
	theta(45) = 1.210e-4
	theta(46) = 9.666e-5
	theta(47) = 7.705e-5
	theta(48) = 6.148e-5
	theta(49) = 4.904e-5
	theta(50) = 3.909e-5
	theta(51) = 3.112e-5
	theta(52) = 2.473e-5
	theta(53) = 1.997e-5
	theta(54) = 1.578e-5
	theta(55) = 1.244e-5
	theta(56) = 9.769e-6
	theta(57) = 7.634e-6
	theta(58) = 5.932e-6
	theta(59) = 4.581e-6
	theta(60) = 3.515e-6
	theta(61) = 2.679e-6
	theta(62) = 2.029e-6
	theta(63) = 1.527e-6
	theta(64) = 1.144e-6
	theta(65) = 8.523e-7
	theta(66) = 6.332e-7
	theta(67) = 4.693e-7
	theta(68) = 3.475e-7
	theta(69) = 2.574e-7
	theta(70) = 1.907e-7
	theta(71) = 1.413e-7

	theta(72) = 1.047e-7
	theta(73) = 7.739e-8
	theta(74) = 5.677e-8
	theta(75) = 4.386e-8
	theta(76) = 3.227e-8
	theta(77) = 2.385e-8
	theta(78) = 1.750e-8
	theta(79) = 1.248e-8
	theta(80) = 8.389e-9
	theta(81) = 5.026e-9
	theta(82) = 2.382e-9
	theta(83) = 6.259e-10
	theta(84) = 1.653e-10
	theta(85) = 7.399e-11
	theta(86) = 4.824e-11
	theta(87) = 3.474e-11
	theta(88) = 2.633e-11
	theta(89) = 2.069e-11
	theta(90) = 1.663e-11
	theta(91) = 1.099e-11
	theta(92) = 7.506e-12
	theta(93) = 4.825e-12
	theta(94) = 2.864e-12
	theta(95) = 1.534e-12
	theta(96) = 7.324e-13
	theta(97) = 3.087e-13
	theta(98) = 1.135e-13
	theta(99) = 3.591e-14
	theta(100) = 9.689e-15
	theta(101) = 2.045e-15
	theta(102) = 2.618e-16
	theta(103) = 8.918e-18
	theta(104) = 3.041e-21
	theta(105) = 1.739e-23

	return
	end

C..............................................................................
C
C This function calculates the CO photodissociation rate and takes into account
C the CO self-shielding.  Isotopes of C (Y = 13C) and O (Z = 18O) are included.
C Based heavily on Ted Bergin's self_shield_v2.f code.
C
C Uses rates and k_Av0 from Visser et al. 2009, A&A 503, 323
C
C..............................................................................
C
C Input parameter(s):
C
C abundances		== array storing the chemical abundances
C
C..............................................................................
	DOUBLE PRECISION FUNCTION self_shield_CO(iso)
	implicit none

C Initialization of common blocks:
	INCLUDE "Fcn.h"
	INCLUDE "rates.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"

	DOUBLE PRECISION iso
	double precision xnh2(ih2),xnco(ico)
	double precision theta(ih2,ico), dtheta(ih2,ico)
	double precision theta_highco(ih2), dtheta_highco(ih2)
	double precision theta_highh2(ico), dtheta_highh2(ico)
	double precision thetaco(ih2,ico), dthetaco(ih2,ico)

	double precision colden_co, colden_13co, colden_c18o, colden_13c18o
	double precision colden_h2_12, colden_h2_13, colden_h2_18, colden_h2_1318
	double precision colden_co_tot(ntime), colden_13co_tot(ntime)
	double precision colden_c18o_tot(ntime), colden_13c18o_tot(ntime)
	double precision colden_h2_tot_12(ntime), colden_h2_tot_13(ntime)
	double precision colden_h2_tot_18(ntime), colden_h2_tot_1318(ntime)
	DATA colden_co_tot/ntime * 0.0/, colden_13co_tot/ntime * 0.0/
	DATA colden_c18o_tot/ntime * 0.0/, colden_13c18o_tot/ntime * 0.0/
	DATA colden_h2_tot_12/ntime * 0.0/, colden_h2_tot_13/ntime * 0.0/
	DATA colden_h2_tot_18/ntime * 0.0/, colden_h2_tot_1318/ntime * 0.0/

	double precision thetatot
	double precision nCO, H2abun
	double precision k_Av0, etau, maxlam
	double precision lnh2, lnco, t, u
	double precision y1,y2,y3,y4

	double precision lnc18o

	double PRECISION G0eval
	double precision test1,test2

	integer i, j, i_h2, i_co
	integer flag
	integer COindex
	integer smallcol

	IF (shieldtest) THEN
		OPEN(unit=02, file=fnameSS, status='old', ACCESS ='APPEND')
		IF (timestep.EQ.1) THEN
			WRITE(02,*)
			WRITE(02,*) 'Zone: ', zone
		END IF
	END IF

c Zero out theta and dtheta.  ico and ih2 are parameters set in rates.h that determine the table size.
	do i = 1,ih2
	  do j = 1, ico
		theta(i,j) = 0.0
		dtheta(i,j) = 0.0
	  end do
	end do

c  Read in the visser self-shielding parameter tables, which will then be spline interpolated
c change this
	if (iso .eq. 12)  call visser_coshield(theta, xnco, xnh2)
	if (iso .eq. 13)  call visser_13coshield(theta, xnco, xnh2)
	if (iso .eq. 18)  call visser_c18oshield(theta, xnco, xnh2)
	if (iso .eq. 38)  call visser_13c18oshield(theta, xnco, xnh2)

c/***
c               setup the table for the spline by calling the splie2 routine
c               returns the array dtheta containing derivitives
c***/
	call splie2(xnh2,xnco,theta,ih2,ico,dtheta)
c/***
c               also have to have it set up to handle limits such as
c               h2_colden > 10^23 and co_colden < 10^19 so interpolate
c               only in one direction.  of course we have to have it set
c               up for both ways.
c***/


c/ First handle theta for high nH2 values:
	do i = 1, ico
		theta_highh2(i) = theta(ih2,i)
	end do

	call spline(xnco,theta_highh2,ico,1.0e30,1.0e30,dtheta_highh2)

c Then handle theta for high nCO values:
	do i = 1, ih2
		theta_highco(i) = theta(i,ico)
	end do

	call spline(xnh2,theta_highco,ih2,1.0e30,1.0e30,dtheta_highco)

C Get the current CO abundance
	if (iso .eq. 12) then
		CALL ispecies('CO           ', ns, s, COindex)
		IF (timestep .EQ. 1) THEN
			nCO = abundances(COindex, timestep)
		ELSE
			nCO = abundances(COindex, timestep-1)
		END IF
	end if
        if (iso .eq. 13) then
                CALL ispecies('CO           ', ns, s, COindex)
                IF (timestep .EQ. 1) THEN
                        nCO = abundances(COindex, timestep)
                ELSE
                        nCO = abundances(COindex, timestep-1)
                END IF
        end if
        if (iso .eq. 18) then
                CALL ispecies('CO           ', ns, s, COindex)
                IF (timestep .EQ. 1) THEN
                        nCO = abundances(COindex, timestep)
                ELSE
                        nCO = abundances(COindex, timestep-1)
                END IF
        end if
        if (iso .eq. 38) then
                CALL ispecies('CO           ', ns, s, COindex)
                IF (timestep .EQ. 1) THEN
                        nCO = abundances(COindex, timestep)
                ELSE
                        nCO = abundances(COindex, timestep-1)
                END IF
        end if

C Change abundance for different isotopes
c$$$	IF (iso .eq. 13) THEN
c$$$	   nCO=nCO/75.0
c$$$	ELSE IF (iso .eq. 18) THEN
c$$$	   nCO=nCO/500.0
c$$$	END IF

c/***
c   Compute the shielding due to the dust continuum
C   this section is no longer used:
c***/
c	print *,'g0_CO Rs(zone): ',g0_CO,Rs(zone)
c	if (iso .eq. 12) k_Av0 = 2.590e-10 * g0_CO * (100.0 / Rs(zone))**2.0
c	if (iso .eq. 13) k_Av0 = 2.597e-10 * g0_CO * (100.0 / Rs(zone))**2.0
c	if (iso .eq. 18) k_Av0 = 2.396e-10 * g0_CO * (100.0 / Rs(zone))**2.0
c	if (iso .eq. 38) k_Av0 = 2.510e-10 * g0_CO * (100.0 / Rs(zone))**2.0

C Updated to calculate G0_CO on the fly as well 9/18/13.  These *do not* need etau -- include extinction already!
	maxlam = 1100.0 ! updated from 960
	call find_G0(maxlam, G0eval, zone, 0)

	if (iso .eq. 12) k_Av0 = 2.590e-10 * G0eval
	if (iso .eq. 13) k_Av0 = 2.597e-10 * G0eval
	if (iso .eq. 18) k_Av0 = 2.396e-10 * G0eval
	if (iso .eq. 38) k_Av0 = 2.510e-10 * G0eval


c/***
c   Here we have to save the column density of H2 and CO at each time step
c   and add it to the column density in the previous zone.
c   This requires the time steps to be the same in each zone!!!! which
c   should be the case - but put in a check just in case.  the check is
c   located in the main program.
c***/

	if (zone.eq.1) zcm(zone-1) = 0.0
	if (iso .eq. 12) then

c Column density in current zone:
		colden_co = nCO * (zcm(zone) - zcm(zone-1))
		if (nCO .LE. 0.0) colden_co = 1.0e4
		colden_co_tot(timestep) = colden_co_tot(timestep) + colden_co
		lnco = dlog10(colden_co_tot(timestep))
		H2abun = rho(zone) / mH / amu
		colden_h2_12 = H2abun * (zcm(zone) - zcm(zone-1))
		colden_h2_tot_12(timestep) = colden_h2_tot_12(timestep)
     &		+ colden_h2_12
		lnh2 = dlog10(colden_h2_tot_12(timestep))
	end if

	if (iso .eq. 13) then
		colden_13co = nCO * (zcm(zone) - zcm(zone-1))
		if (nCO .LE. 0.0) colden_13co = 1.0e4
		colden_13co_tot(timestep) = colden_13co_tot(timestep)
     &		+ colden_13co
		lnco = dlog10(colden_13co_tot(timestep))
		H2abun = rho(zone) / mH / amu
		colden_h2_13 = H2abun * (zcm(zone) - zcm(zone-1))
		colden_h2_tot_13(timestep) = colden_h2_tot_13(timestep)
     &		+ colden_h2_13
		lnh2 = dlog10(colden_h2_tot_13(timestep))
	end if

	if (iso .eq. 18) then
		colden_c18o = nCO * (zcm(zone) - zcm(zone-1))
		if (nCO .LE. 0.0) colden_c18o = 1.0e4
		colden_c18o_tot(timestep) = colden_c18o_tot(timestep)
     &		+ colden_c18o
		lnco = dlog10(colden_c18o_tot(timestep))
		H2abun = rho(zone) / mH / amu
		colden_h2_18 = H2abun * (zcm(zone) - zcm(zone-1))
		colden_h2_tot_18(timestep) = colden_h2_tot_18(timestep)
     &		+ colden_h2_18
		lnh2 = dlog10(colden_h2_tot_18(timestep))
	end if

C  Temporary testing to find C18O self-shielding factor:
	lnc18o = dlog10((10**lnco)/500.0)
c	write(*,*) 'lnc18o: ',lnco,lnc18o,lnh2,k_Av0

	if (iso .eq. 38) then
		colden_13c18o = nCO * (zcm(zone) - zcm(zone-1))
		if (nCO .LE. 0.0) colden_13c18o = 1.0e4
		colden_13c18o_tot(timestep) = colden_13c18o_tot(timestep)
     &		+ colden_13c18o
		lnco = dlog10(colden_13c18o_tot(timestep))
		H2abun = rho(zone) / mH
		colden_h2_1318 = H2abun * (zcm(zone) - zcm(zone-1))
		colden_h2_tot_1318(timestep) = colden_h2_tot_1318(timestep)
     &		+ colden_h2_1318
		lnh2 = dlog10(colden_h2_tot_1318(timestep))
	end if

C  "flag" specifies whether the data is within the Visser table range (flag=0).
C  The following sequence of if-statements handles the special case of large
c  or small column.  When flag=0, the value is taken from the table.
	flag = 0
	if ((lnco .le. 13.0) .and. (lnh2 .le. 13.0)) then
		thetatot = theta(1,1)
		self_shield_CO = theta(1,1) * k_Av0
		flag = 1
	end if

	if ((lnco .ge. 19.0) .and. (lnh2 .ge. 23.0)) then
		thetatot = theta(6,8)
		self_shield_CO = theta(6,8) * k_Av0
		flag = 1
	else if (lnco .ge. 19.0) then
		call splint(xnh2,theta_highco,dtheta_highco,ih2,lnh2,thetatot)

c spline between 22 and 23 is negative, so linearly interpolate in that range
        if (thetatot .LT. 0) then
            do i = 1,ih2
                if (lnh2 .GT. xnh2(i)) i_h2 = i
            end do
            thetatot = theta(i_h2,8) + (lnh2 - floor(lnh2)) *
     &		(theta(i_h2+1,8) - theta(i_h2,8)) / (ceiling(lnh2) - floor(lnh2))
        endif

        self_shield_CO = thetatot * k_Av0
        flag = 1

	else if (lnh2 .ge. 23.0) then
		call splint(xnco,theta_highh2,dtheta_highh2,
     +          ico, lnco, thetatot)
		self_shield_CO = thetatot * k_Av0
		flag = 1
	end if

	smallcol = 0

	if (flag .eq. 0) then

		if (lnco .LT. 0.0) then
			print *, 'FYI: lnco is less than zero - setting it to 0.0.'
			lnco = 0.0
			smallcol = 1
		end if

		do i = 1,ih2
			if (lnh2 .GE. xnh2(i)) i_h2 = i
		end do

		do i = 1,ico
			if (lnco .GE. xnco(i)) i_co = i
		end do

		t = (lnh2-xnh2(i_h2))/(xnh2(i_h2+1)-xnh2(i_h2))
		u = (lnco-xnco(i_co))/(xnco(i_co+1)-xnco(i_co))
		thetatot = 0

		y1 = theta(i_h2,i_co)
		y2 = theta(i_h2+1,i_co)
		y3 = theta(i_h2+1,i_co+1)
		y4 = theta(i_h2,i_co+1)

		thetatot = (1-t)*(1-u)*y1
		thetatot = thetatot + t*(1-u)*y2
		thetatot = thetatot + t*u*y3
		thetatot = thetatot + (1-t)*u*y4

c		thetatot = 1.0
c theta = 1 equivalent to turning off self shielding
		self_shield_CO = k_Av0 * thetatot
C		print *, 'k_Av0 thetatot: ',k_Av0,thetatot
	end if
	flag = 0

c   ETAU (Outdated, do not use).
c	maxlam = 960.0
c	call find_etau(maxlam, etau, zone)

C	print *, 'columns lnco, lnh2:',lnco,lnh2,thetatot,self_shield_CO,self_shield_CO * etau
c	self_shield_CO = self_shield_CO * etau


	IF (shieldtest) THEN
		if (iso .eq. 12) then
			write(2, 112) iso,zcm(zone),self_shield_CO,
     &			colden_co_tot(timestep), colden_h2_tot_12(timestep),
     &			thetatot, etau, nCO, timestep
		end if
		if (iso .eq. 13) then
			write(2, 112) iso,zcm(zone), self_shield_CO,
     &			colden_13co_tot(timestep), colden_h2_tot_13(timestep),
     &			thetatot, etau, nCO, timestep
		end if
		if (iso .eq. 18) then
			write(2, 112) iso,zcm(zone), self_shield_CO,
     &			colden_c18o_tot(timestep), colden_h2_tot_18(timestep),
     &			thetatot, etau, nCO, timestep
		end if
		if (iso .eq. 38) then
			write(2, 112) iso,zcm(zone), self_shield_CO,
     &			colden_13c18o_tot(timestep), colden_h2_tot_1318(timestep),
     &			thetatot, etau, nCO, timestep
		end if

		CLOSE(02)
	END IF

 112     format(1x,f3.0,2x,7(1x,e8.3),2x, i4)

	if (self_shield_CO .LT. 0) then
		print *, 'You broke the code!! selfshield_CO: ',self_shield_CO
		if (smallcol .eq. 1) then
		     self_shield_CO = 0.0 !! LIC?
		else
			if (abs(self_shield_CO) .lt. 1e-5) then
				self_shield_CO = 0.0
			else
			   print *,'Exiting now.'
			   stop
			endif
		end if
	end if

	return
	end


C.......................................................................
C FIND_ETAU - Interpolate the value of the UV field at current z
C
C
C.......................................................................
	subroutine find_etau(x, etau, jzone)

	INCLUDE "environ.h"

	DOUBLE PRECISION x, etau
	INTEGER jzone

	double precision dlam, dphotuv, photuv
	double precision dphotuv_max, photuv_max
	integer i, lambdai

	do i = 1, nwavl
		if(lambda(i) .gt. x) then
			ilam = i
			exit
		end if
	end do

	dlam = lambda(i) - lambda(i-1)
	dphotuv = ((x-lambda(i-1))/dlam) *
     &		  (uvfield(jzone,i)-uvfield(jzone,i-1))
	photuv = uvfield(jzone,i-1) + dphotuv

C zone 1 = top of disk = I_0
C	dphotuv0 = ((x-lambda(i-1))/dlam) *
C     &		   (uvfield(1,i)-uvfield(1,i-1))
C	photuv0 = uvfield(1,i-1) + dphotuv0

C UVmaxzone(i) = zone with maximum flux at wavelength x
c	print *, lambda(i-1), UVmaxzone(i-1)
	dphotuv_max = ((x-lambda(i-1))/dlam) *
     &		   (uvfield(UVmaxzone(i-1),i)-uvfield(UVmaxzone(i-1),i-1))
	photuv_max = uvfield(UVmaxzone(i-1),i-1) + dphotuv_max

C etau = I / I_max
c	if (x .EQ. 960) then
c		print *, 'IN ETAU: ', jzone, photuv, photuv_max
c	endif
	etau = photuv/photuv_max

	return
	end


C.......................................................................
C FIND_G0 - Evaluate G0 in habings at current z
C
C  lyaflag == 1 means calculate g0_100, 0 means g0_CO
C.......................................................................
        subroutine find_G0(x, G0eval, jzone, lyaflag)

        INCLUDE "environ.h"

        DOUBLE PRECISION x, etau
		DOUBLE PRECISION G0eval,totuv,En_i,G0_Habing
        INTEGER jzone, lyaflag

        double precision dlam, dphotuv, photuv
        double precision dphotuv_max, photuv_max
        integer i, lambdai

        totuv = 0.0
        do i = 1, nwavl
			if(lambda(i) .le. x) then
			   En_i = 6.6260755e-27*2.99792458e10/(lambda(i)*1e-8)
			   if (lambda(i) .EQ. 1210) then
				   if (lyaflag .EQ. 1) then
					   totuv = uvfield(jzone,i) * 10.0 * En_i + totuv
				   else	 ! Interpolate over Lyman alpha.
					   totuv = (uvfield(jzone,i-1) +
     &		                uvfield(jzone,i+1))/2 * 10.0 * En_i + totuv
				   end if
			   else
				   totuv = uvfield(jzone,i) * 10.0 * En_i + totuv
			   end if
			else
c			if(lambda(i) .gt. x) then
c				ilam = i
c				totuv=totuv-(uvfield(jzone,i)*(lambda(i)-x))*En_i ! Subtract the extra over limit x.
c						En_i = 6.6260755e-27*2.99792458e10/(lambda(i)*1e-8)
c						totuv = uvfield(jzone,i) * 10.0 * En_i + totuv  ! erg/cm^2/s units.
				exit
			end if
        end do

        G0_Habing = 1.6e-3 ! erg/cm^2/s from Herczeg et al. 2004
		G0eval = totuv/G0_Habing

        return
        end

C..............................................................................
C VISSER_COSHIELD - returns xnco, xnh2 and theta values for CO
C
C From Visser et al. 2009, A&A, 503, 323
C b(CO) = 0.3 km/s, Tex(CO) = 50 K, N(12CO)/N(13CO) = 69
C..............................................................................
	subroutine visser_coshield(theta, xnco, xnh2)
	INCLUDE "rates.h"
	DOUBLE PRECISION theta(ih2,ico), xnco(ico), xnh2(ih2)

	xnco(1) = 0.0
	xnco(2) = 13.0
	xnco(3) = 14.0
	xnco(4) = 15.0
 	xnco(5) = 16.0
	xnco(6) = 17.0
	xnco(7) = 18.0
	xnco(8) = 19.0

	xnh2(1) = 0.0
	xnh2(2) = 19.0
	xnh2(3) = 20.0
	xnh2(4) = 21.0
	xnh2(5) = 22.0
	xnh2(6) = 23.0

	theta(1,1) = 1.00000
	theta(1,2) = 9.405e-1
	theta(1,3) = 7.046e-1
	theta(1,4) = 4.015e-1
	theta(1,5) = 9.964e-2
	theta(1,6) = 1.567e-2
	theta(1,7) = 3.162e-3
	theta(1,8) = 4.839e-4

	theta(2,1) = 7.546e-1
	theta(2,2) = 6.979e-1
	theta(2,3) = 4.817e-1
	theta(2,4) = 2.577e-1
	theta(2,5) = 6.505e-2
	theta(2,6) = 1.135e-2
	theta(2,7) = 2.369e-3
	theta(2,8) = 3.924e-4

	theta(3,1) = 5.752e-1
	theta(3,2) = 5.228e-1
	theta(3,3) = 3.279e-1
	theta(3,4) = 1.559e-1
	theta(3,5) = 3.559e-2
	theta(3,6) = 6.443e-3
	theta(3,7) = 1.526e-3
	theta(3,8) = 2.751e-4

	theta(4,1) = 2.493e-1
	theta(4,2) = 2.196e-1
	theta(4,3) = 1.135e-1
	theta(4,4) = 4.062e-2
	theta(4,5) = 7.864e-3
	theta(4,6) = 1.516e-3
	theta(4,7) = 4.448e-4
	theta(4,8) = 9.367e-5

	theta(5,1) = 1.550e-3
	theta(5,2) = 1.370e-3
	theta(5,3) = 6.801e-4
	theta(5,4) = 2.127e-4
	theta(5,5) = 5.051e-5
	theta(5,6) = 1.198e-5
	theta(5,7) = 6.553e-6
	theta(5,8) = 3.937e-6

	theta(6,1) = 8.492e-8
	theta(6,2) = 8.492e-8
	theta(6,3) = 8.492e-8
	theta(6,4) = 8.492e-8
	theta(6,5) = 8.492e-8
	theta(6,6) = 8.492e-8
	theta(6,7) = 8.488e-8
	theta(6,8) = 8.453e-8

	return
	end

C..............................................................................
C VISSER_13COSHIELD - returns xnco, xnh2 and theta values for 13CO
C
C From Visser et al. 2009, A&A, 503, 323
C b(CO) = 0.3 km/s, Tex(CO) = 50 K, N(12CO)/N(13CO) = 69
C..............................................................................
	subroutine visser_13coshield(theta, xnco, xnh2)
	INCLUDE "rates.h"
	DOUBLE PRECISION theta(ih2,ico), xnco(ico), xnh2(ih2)
c	print *, 'In 13coshield'
	xnco(1) = 0.0
	xnco(2) = 13.0
	xnco(3) = 14.0
	xnco(4) = 15.0
 	xnco(5) = 16.0
	xnco(6) = 17.0
	xnco(7) = 18.0
	xnco(8) = 19.0

	xnh2(1) = 0.0
	xnh2(2) = 19.0
	xnh2(3) = 20.0
	xnh2(4) = 21.0
	xnh2(5) = 22.0
	xnh2(6) = 23.0

	theta(1,1) = 1.00000
	theta(1,2) = 9.765e-1
	theta(1,3) = 8.965e-1
	theta(1,4) = 7.701e-1
	theta(1,5) = 4.459e-1
	theta(1,6) = 1.415e-1
	theta(1,7) = 1.748e-2
	theta(1,8) = 8.346e-4

	theta(2,1) = 7.780e-1
	theta(2,2) = 7.550e-1
	theta(2,3) = 6.789e-1
	theta(2,4) = 5.753e-1
	theta(2,5) = 3.262e-1
	theta(2,6) = 1.032e-1
	theta(2,7) = 1.183e-2
	theta(2,8) = 6.569e-4

	theta(3,1) = 5.523e-1
	theta(3,2) = 5.309e-1
	theta(3,3) = 4.617e-1
	theta(3,4) = 3.793e-1
	theta(3,5) = 1.958e-1
	theta(3,6) = 5.888e-2
	theta(3,7) = 7.515e-3
	theta(3,8) = 4.382e-4

	theta(4,1) = 2.100e-1
	theta(4,2) = 1.979e-1
	theta(4,3) = 1.601e-1
	theta(4,4) = 1.244e-1
	theta(4,5) = 5.373e-2
	theta(4,6) = 1.573e-2
	theta(4,7) = 2.535e-3
	theta(4,8) = 1.361e-4

	theta(5,1) = 1.318e-3
	theta(5,2) = 1.268e-3
	theta(5,3) = 1.107e-3
	theta(5,4) = 9.017e-4
	theta(5,5) = 3.658e-4
	theta(5,6) = 1.114e-4
	theta(5,7) = 2.414e-5
	theta(5,8) = 2.608e-6

	theta(6,1) = 4.511e-8
	theta(6,2) = 4.511e-8
	theta(6,3) = 4.511e-8
	theta(6,4) = 4.511e-8
	theta(6,5) = 4.511e-8
	theta(6,6) = 4.511e-8
	theta(6,7) = 4.509e-8
	theta(6,8) = 4.490e-8

	return
	end

C..............................................................................
C VISSER_C18OSHIELD - returns xnco, xnh2 and theta values for C18O
C
C From Visser et al. 2009, A&A, 503, 323
C b(CO) = 0.3 km/s, Tex(CO) = 50 K, N(12CO)/N(13CO) = 69
C..............................................................................
	subroutine visser_c18oshield(theta, xnco, xnh2)
	INCLUDE "rates.h"
	DOUBLE PRECISION theta(ih2,ico), xnco(ico), xnh2(ih2)
c	print *, 'In c18oshield'
	xnco(1) = 0.0
	xnco(2) = 13.0
	xnco(3) = 14.0
	xnco(4) = 15.0
 	xnco(5) = 16.0
	xnco(6) = 17.0
	xnco(7) = 18.0
	xnco(8) = 19.0

	xnh2(1) = 0.0
	xnh2(2) = 19.0
	xnh2(3) = 20.0
	xnh2(4) = 21.0
	xnh2(5) = 22.0
	xnh2(6) = 23.0

	theta(1,1) = 1.00000
	theta(1,2) = 9.822e-1
	theta(1,3) = 9.163e-1
	theta(1,4) = 8.067e-1
	theta(1,5) = 5.498e-1
	theta(1,6) = 2.188e-1
	theta(1,7) = 3.412e-2
	theta(1,8) = 1.992e-3

	theta(2,1) = 8.007e-1
	theta(2,2) = 7.833e-1
	theta(2,3) = 7.206e-1
	theta(2,4) = 6.272e-1
	theta(2,5) = 4.226e-1
	theta(2,6) = 1.597e-1
	theta(2,7) = 2.469e-2
	theta(2,8) = 1.630e-3

	theta(3,1) = 5.848e-1
	theta(3,2) = 5.688e-1
	theta(3,3) = 5.123e-1
	theta(3,4) = 4.363e-1
	theta(3,5) = 2.816e-1
	theta(3,6) = 9.826e-2
	theta(3,7) = 1.653e-2
	theta(3,8) = 1.158e-3

	theta(4,1) = 2.277e-1
	theta(4,2) = 2.187e-1
	theta(4,3) = 1.881e-1
	theta(4,4) = 1.546e-1
	theta(4,5) = 9.811e-2
	theta(4,6) = 3.114e-2
	theta(4,7) = 5.506e-3
	theta(4,8) = 3.655e-4

	theta(5,1) = 1.411e-3
	theta(5,2) = 1.375e-3
	theta(5,3) = 1.256e-3
	theta(5,4) = 1.126e-3
	theta(5,5) = 8.200e-4
	theta(5,6) = 2.912e-4
	theta(5,7) = 5.310e-5
	theta(5,8) = 5.241e-6

	theta(6,1) = 8.937e-8
	theta(6,2) = 8.937e-8
	theta(6,3) = 8.937e-8
	theta(6,4) = 8.937e-8
	theta(6,5) = 8.937e-8
	theta(6,6) = 8.936e-8
	theta(6,7) = 8.933e-8
	theta(6,8) = 8.896e-8

	return
	end

C..............................................................................
C VISSER_13C18OSHIELD - returns xnco, xnh2 and theta values for 13CO
C added 06/25/14 - KS
C From Visser et al. 2009, A&A, 503, 323, Table 6
C b(CO) = 0.3 km/s, Tex(CO) = 50 K, N(12CO)/N(13CO) = 69
C..............................................................................
	subroutine visser_13c18oshield(theta, xnco, xnh2)
	INCLUDE "rates.h"
	DOUBLE PRECISION theta(ih2,ico), xnco(ico), xnh2(ih2)
c	print *, 'In 13c18oshield'
	xnco(1) = 0.0
	xnco(2) = 13.0
	xnco(3) = 14.0
	xnco(4) = 15.0
 	xnco(5) = 16.0
	xnco(6) = 17.0
	xnco(7) = 18.0
	xnco(8) = 19.0

	xnh2(1) = 0.0
	xnh2(2) = 19.0
	xnh2(3) = 20.0
	xnh2(4) = 21.0
	xnh2(5) = 22.0
	xnh2(6) = 23.0

	theta(1,1) = 1.00000
	theta(1,2) = 9.873e-1
	theta(1,3) = 9.424e-1
	theta(1,4) = 8.744e-1
	theta(1,5) = 6.740e-1
	theta(1,6) = 3.885e-1
	theta(1,7) = 1.268e-1
	theta(1,8) = 2.980e-2

	theta(2,1) = 7.975e-1
	theta(2,2) = 7.852e-1
	theta(2,3) = 7.424e-1
	theta(2,4) = 6.894e-1
	theta(2,5) = 5.447e-1
	theta(2,6) = 3.236e-1
	theta(2,7) = 1.036e-1
	theta(2,8) = 2.722e-2

	theta(3,1) = 5.869e-1
	theta(3,2) = 5.745e-1
	theta(3,3) = 5.364e-1
	theta(3,4) = 4.953e-1
	theta(3,5) = 3.931e-1
	theta(3,6) = 2.445e-1
	theta(3,7) = 8.070e-2
	theta(3,8) = 1.894e-2

	theta(4,1) = 2.465e-1
	theta(4,2) = 2.401e-1
	theta(4,3) = 2.187e-1
	theta(4,4) = 2.011e-1
	theta(4,5) = 1.628e-1
	theta(4,6) = 1.079e-1
	theta(4,7) = 3.267e-2
	theta(4,8) = 4.406e-3

	theta(5,1) = 2.164e-3
	theta(5,2) = 2.138e-3
	theta(5,3) = 2.050e-3
	theta(5,4) = 1.971e-3
	theta(5,5) = 1.772e-3
	theta(5,6) = 1.421e-3
	theta(5,7) = 5.147e-4
	theta(5,8) = 4.066e-5

	theta(6,1) = 4.454e-8
	theta(6,2) = 4.454e-8
	theta(6,3) = 4.454e-8
	theta(6,4) = 4.454e-8
	theta(6,5) = 4.454e-8
	theta(6,6) = 4.452e-8
	theta(6,7) = 4.434e-8
	theta(6,8) = 4.490e-8

	return
	end


C..............................................................................
C SPLIN2 - Given X1A, X2A, YA, M, N as described in SPLIE2 and Y2A as produced
C by that routine, and given a desired interpolating point X1, X2, this routine
C returns an interpolated function value Y by bicubic spline interpolation.
C
C "two-dimensional spline interpolation"
C
C Description from IDL function
C Code from Numerical Recipes
C..............................................................................
      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
      implicit double precision (a-h, o-z)
      PARAMETER (NN=100)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN),YYTMP(
     *NN)
      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
          Y2TMP(K)=Y2A(J,K)
11      CONTINUE
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J))
12    CONTINUE
      CALL SPLINE(X1A,YYTMP,M,1.E30,1.E30,Y2TMP)
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)
      RETURN
      END

C..............................................................................
C SPLIE2 - Given an M by N tabulated function YA, and tabulated independent
C variables X1A (M values) and X2A (N values), this routine constructs
C one-dimensional natural cubic splines of the rows of YA and returns the
C second derivatives in the array Y2A.
C
C "construct two-dimensional spline"
C
C Description from IDL function
C Code from Numerical Recipes
C..............................................................................
      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
      implicit double precision (a-h,o-z)
      PARAMETER (NN=100)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)
      DO 13 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
        CALL SPLINE(X2A,YTMP,N,1.E30,1.E30,Y2TMP)
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
13    CONTINUE
      RETURN
      END

C..............................................................................
C SPLINE - "construct a cubic spline"
C
C Code from Numerical Recipes
C..............................................................................
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      implicit double precision (a-h, o-z)
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

C..............................................................................
C SPLINT - "cubic spline interpolation"
C
C Code from Numerical Recipes
C..............................................................................
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      implicit double precision (a-h, o-z)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) THEN
		write(*,*) 'ERROR: Bad XA input in splint'
		STOP
	  END IF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END

C..............................................................................
C
C This function calculates the needed values for monolayer adjustments
C TIME DEPENDENT due to the ngr term and n_ice.
C used for types -23, -24 (photodesorption), 21 (classical evaporation)
C and 22 (CR desorption)
C
C NOT USED AS OF 6/22/09 due to different calculation done in
C photodesorp_rate
C
C..............................................................................
	SUBROUTINE monolayer_calc()
	implicit none

C Initialization of common blocks:
	INCLUDE "Fcn.h"
	INCLUDE "rates.h"
	INCLUDE "constants.h"
        INCLUDE "environ.h"

	DOUBLE PRECISION a_gr, Mlayer, ngr
	INTEGER i, grindex
        DOUBLE PRECISION sigadjust, rgr

        sigadjust = 1.0
        if (incl_locdust) sigadjust = locdust(zone)


c sum the abundance of on-grain species
	i = 0
	do while (grspecs(i) .NE. -1)
		IF (timestep .EQ. 1) THEN
			n_ice = n_ice + abundances(grspecs(i), timestep)
		ELSE
			n_ice = n_ice + abundances(grspecs(i), timestep-1)
		ENDIF
		i = i+1
	end do

	CALL ispecies('GRAIN        ', ns, s, grindex)
	IF (timestep .EQ. 1) THEN
		ngr = abundances(grindex, timestep)
	ELSE
		ngr = abundances(grindex, timestep-1)
	END IF

C catch for first run, before any grain species have been calculated
C (before first solver call).
	IF (n_ice .LT. 1.0e-20) THEN
		n_ice = n_ice_init
	ENDIF
	IF (ngr .LT. 1.0e-20) THEN
		ngr = ngr_init
	ENDIF

	!a_gr = rgr*freezeeffic**(-1.0/1.5)  ! The freezeeffic term accounts for grain growth (small freeze effic, bigger grains, but less surf area).
	Mlayer = 4 * pi * (rgr*1e-4)**2 * sigadjust / ((3e-8)**2)
	numlayers = n_ice / ngr / Mlayer
	print *, '# of monolayers = ', numlayers

C perc_ice = n(i) / n_ice
C fac = ngr * perc_ice / n(i), perc_ice = n(i) / n_ice
C fac = PDadjust
	IF (numlayers .GT. 1.0) THEN
		PDadjust = ngr / n_ice
c		if (perc_ice(rspec(i,1)) .lt. 1e-3) fac = 1.0e-3
	ELSE
		PDadjust = ngr * numlayers / n_ice
c		if (perc_ice(rspec(i,1)) .lt. 5.0e-5) fac = 5.0e-5
	ENDIF

	IF (numlayers .LT. 1.0) THEN
		Madjust = 1.0
	ELSE
		Madjust = 1.0 / numlayers
	ENDIF

c	print *, 'Madjust = ', Madjust
c	print *, 'PDadjust = ', PDadjust

	RETURN
	END
C..............................................................................
C
C This function calculates the HD photodissociation rate and takes into account
C the H2 self-shielding.
C Uses rates and k_Av0 from Lee et al. 1996, A&A, 311, 690
C..............................................................................
	DOUBLE PRECISION FUNCTION self_shield_HD()
	implicit none

C Initialization of common blocks:
	INCLUDE "Fcn.h"
	INCLUDE "rates.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"

	DOUBLE PRECISION colden_h2_hd, colden_h2_tot_hd(ntime)
	DOUBLE PRECISION colden_h_hd, colden_h_tot_hd(ntime)
	DOUBLE PRECISION colden_hd, colden_hd_tot(ntime)
	DOUBLE PRECISION H2abun, maxlam, etau
	DOUBLE PRECISION nHD, HDabun, nHI
	INTEGER n_lee_hd, i
	PARAMETER (n_lee_hd = 105)
	DOUBLE PRECISION lee_nh2_hd(n_lee_hd),lee_theta_hd(n_lee_hd),theta
	DOUBLE PRECISION lnh2(n_lee_hd), interp
	DOUBLE PRECISION itheta
	INTEGER HDindex, HIindex
	DOUBLE PRECISION G0eval
	DOUBLE PRECISION zet_2,x_2,del_2,alph_2,f_2,zet_1,x_1,del_1,alph_1,f_1

	DATA colden_hd/0.0/,colden_h2_hd/0.0/,colden_h2_tot_hd/ntime * 0.0/
	DATA colden_hd/0.0/, colden_hd_tot/ntime * 0.0/
	DATA lee_nh2_hd/n_lee_hd * 0.0/, lee_theta_hd/n_lee_hd * 0.0/
	DATA colden_h_hd/0.0/,colden_h_tot_hd/ntime * 0.0/

	IF (zone .EQ. 1) then
		zcm(zone-1) = 0.0
	endif

	if (timestep .EQ. 1) then
		call lee_factors(lee_nh2_hd, lee_theta_hd)
	endif

c	H2abun = rho(zone) / mH / amu ! LIC fix
c	colden_h2 = H2abun * (zcm(zone) - zcm(zone-1))
c	colden_h2_tot(timestep) = colden_h2_tot(timestep) + colden_h2

C Get the current CO abundance
	CALL ispecies('HD           ', ns, s, HDindex)
	IF (timestep .EQ. 1) THEN
		nHD = abundances(HDindex, timestep)
	ELSE
		nHD = abundances(HDindex, timestep-1)
	END IF


	CALL ispecies('H            ', ns, s, HIindex)
	IF (timestep .EQ. 1) THEN
		nHI = abundances(HIindex, timestep)
	ELSE
		nHI = abundances(HIindex, timestep-1)
	END IF

	colden_hd = nHD * (zcm(zone) - zcm(zone-1))
	colden_hd_tot(timestep) = colden_hd_tot(timestep) + colden_hd

	colden_h2_hd = rho(zone) / mH / amu * (zcm(zone) - zcm(zone-1))
	colden_h2_tot_hd(timestep) = colden_h2_tot_hd(timestep)
     &		+ colden_h2_hd

c	colden_h2_hd = Nrz(zone)  ! Trying radial column density

	colden_h_hd = nHI *  (zcm(zone) - zcm(zone-1))
	colden_h_tot_hd(timestep) = colden_h_tot_hd(timestep) + colden_h_hd




c take log10 of column density table to prepare for interpolation
	do i=1,n_lee_hd
		if (lee_nh2_hd(i) .EQ. 0.0) then
			lnh2(i) = 0.0
		else
			lnh2(i) = dlog10(lee_nh2_hd(i))
		endif
	enddo

c changed to I/I0 at 1000 Ang
	maxlam = 1000.0
	call find_G0(maxlam, G0eval, zone, 1)

c Linear interpolation
	theta = interp(dlog10(colden_hd_tot(timestep)),
     &		 lnh2, lee_theta_hd, n_lee_hd)

C Shielding of HD by H2 (Wolcott-Green and Haimann 2011)

	zet_2 = 2.34e+19
	x_2 = colden_h2_tot_hd(timestep)/zet_2
	del_2 =  0.238
	alph_2 = 5.2e-3

	f_2 = 1.0/(1.0+x_2)**del_2 * exp(-alph_2*x_2)

	zet_1 = 2.85e23
	x_1 = colden_h_tot_hd(timestep)/zet_1
	del_1 = 1.62
	alph_1 = 0.149

	f_1 = 1.0/(1.0+x_1)**del_1 * exp(-alph_1*x_1)

	if (colden_h2_tot_hd(timestep).ge.1e22) f_2 = 0.0
	if (colden_h2_tot_hd(timestep).ge.1e24) f_1 = 0.0

c	write(*,'(A5,E10.3,A5,E10.3,1x,E10.3,1x,E10.3,1x,A5,E10.3)') 'f_2',f_2,'f_1',
c     &		 f_1,colden_hd_tot(timestep),colden_h2_tot_hd(timestep),'th:',theta

c	IF ((timestep.EQ.1).and.(zone.eq.25)) THEN
c		OPEN(unit=15, file='shieldinginfoHD.dat', status='new')
c		WRITE(15,*) '# HD shielding'
c		WRITE(15,*) '# N(HD), N(H2), f1(HI), f2(H2), theta(HD)'
c		close(15)
c	END IF

c	OPEN(unit=15, file='shieldinginfoHD.dat', status='old', ACCESS ='APPEND')
c	WRITE(15,'(E10.3,E10.3,1x,E10.3,1x,E10.3,1x,E10.3)')
c     &		 colden_hd_tot(timestep),colden_h2_tot_hd(timestep),f_1,
c     &		 f_2,theta
c	close(15)

c	write(*,*) 'ColdenH2',colden_h2_tot_hd(timestep)
c	write(*,*) 'ColdenHD',colden_hd_tot(timestep)

	self_shield_HD = 2.54e-11 * G0eval * theta * f_2 * f_1
c	write(*,*) 'HD factor:',theta * f_2 * f_1

	RETURN
	END


C..............................................................................
C
C This function calculates the D2 photodissociation rate and takes into account
C the H2 self-shielding.
C Uses rates and k_Av0 from Lee et al. 1996, A&A, 311, 690
C..............................................................................
	DOUBLE PRECISION FUNCTION self_shield_D2()
	implicit none

C Initialization of common blocks:
	INCLUDE "Fcn.h"
	INCLUDE "rates.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"

	DOUBLE PRECISION colden_h2_d2, colden_h2_tot_d2(ntime)
	DOUBLE PRECISION colden_h_d2, colden_h_tot_d2(ntime)
	DOUBLE PRECISION colden_d2, colden_d2_tot(ntime)
	DOUBLE PRECISION maxlam, etau
	DOUBLE PRECISION nD2, D2abun, nHI
	INTEGER n_lee_d2, i
	PARAMETER (n_lee_d2 = 105)
	DOUBLE PRECISION lee_nh2_d2(n_lee_d2),lee_theta_d2(n_lee_d2),theta
	DOUBLE PRECISION lnh2(n_lee_d2), interp
	DOUBLE PRECISION itheta
	INTEGER D2index, HIindex
	DOUBLE PRECISION G0eval
	DOUBLE PRECISION zet_2,x_2,del_2,alph_2,f_2,zet_1,x_1,del_1,alph_1,f_1

	DATA colden_d2/0.0/,colden_h2_d2/0.0/,colden_h2_tot_d2/ntime * 0.0/

	DATA colden_d2_tot/ntime * 0.0/
	DATA lee_nh2_d2/n_lee_d2 * 0.0/, lee_theta_d2/n_lee_d2 * 0.0/
	DATA colden_h_d2/0.0/,colden_h_tot_d2/ntime * 0.0/

	IF (zone .EQ. 1) then
		zcm(zone-1) = 0.0
	endif

	if (timestep .EQ. 1) then
		call lee_factors(lee_nh2_d2, lee_theta_d2)
	endif

c	H2abun = rho(zone) / mH / amu ! LIC fix
c	colden_h2 = H2abun * (zcm(zone) - zcm(zone-1))
c	colden_h2_tot(timestep) = colden_h2_tot(timestep) + colden_h2

C Get the current CO abundance
	CALL ispecies('D2           ', ns, s, D2index)
	IF (timestep .EQ. 1) THEN
		nD2 = abundances(D2index, timestep)
	ELSE
		nD2 = abundances(D2index, timestep-1)
	END IF

	colden_d2 = nD2 * (zcm(zone) - zcm(zone-1))
	colden_d2_tot(timestep) = colden_d2_tot(timestep) + colden_d2

	colden_h2_d2 = rho(zone) / mH / amu * (zcm(zone) - zcm(zone-1))
	colden_h2_tot_d2(timestep) = colden_h2_tot_d2(timestep)
     &		+ colden_h2_d2

c	colden_h2_d2 = Nrz(zone)  ! Trying radial column density


c	write(*,*) 'ColdenH2',colden_h2_tot_hd(timestep)
c	write(*,*) 'ColdenHD',colden_hd_tot(timestep)

	CALL ispecies('H            ', ns, s, HIindex)
	IF (timestep .EQ. 1) THEN
		nHI = abundances(HIindex, timestep)
	ELSE
		nHI = abundances(HIindex, timestep-1)
	END IF
	colden_h_d2 = nHI *  (zcm(zone) - zcm(zone-1))
	colden_h_tot_d2(timestep) = colden_h_tot_d2(timestep) + colden_h_d2


c take log10 of column density table to prepare for interpolation
	do i=1,n_lee_d2
		if (lee_nh2_d2(i) .EQ. 0.0) then
			lnh2(i) = 0.0
		else
			lnh2(i) = dlog10(lee_nh2_d2(i))
		endif
	enddo

c changed to I/I0 at 1000 Ang
	maxlam = 1000.0
	call find_G0(maxlam, G0eval, zone, 1)

c Linear interpolation
	theta = interp(dlog10(colden_d2_tot(timestep)),
     &		 lnh2, lee_theta_d2, n_lee_d2)

C Shielding of HD by H2 (Wolcott-Green and Haimann 2011)

	zet_2 = 2.34e+19
	x_2 = colden_h2_tot_d2(timestep)/zet_2
	del_2 =  0.238
	alph_2 = 5.2e-3

	f_2 = 1.0/(1.0+x_2)**del_2 * exp(-alph_2*x_2)

	zet_1 = 2.85e23
	x_1 = colden_h_tot_d2(timestep)/zet_1
	del_1 = 1.62
	alph_1 = 0.149

	f_1 = 1.0/(1.0+x_1)**del_1 * exp(-alph_1*x_1)

	if (colden_h2_tot_d2(timestep).ge.1e22) f_2 = 0.0
	if (colden_h2_tot_d2(timestep).ge.1e24) f_1 = 0.0

c	print *, 'f_2',f_2,colden_h2_tot_hd(timestep)

	self_shield_D2 = 2.54e-11 * G0eval * theta * f_2 * f_1

	RETURN
	END
