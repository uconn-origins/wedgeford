C..............................................................................
C
C This function calculates a reaction's rate coefficient
C
C..............................................................................
C
C Input parameter(s):
C
C alpha			== first component of rate coeffs
C
C beta			== second component of rate coeffs
C
C gamma			== third component of rate coeffs
C
C rtype			== type of reaction
C
C r1			== Reactant 1 (used for photo reactions)
C
C timedep		== Flag to calculate time dependent reactions or not
C
C..............................................................................
C
C RATES:
C
C < 0 - time-dependent (self shielding, grain)
C -1 - CO self-shielding
C -2 - H2 self-shielding
C -3 - grain reactions (check to see how many grains are used)
C -23 - Photodesorption 
C -24 - Lyman-alpha photodesorption - same as -23, but with LYAPHOTON
C
C 0 - Grain reaction (r1 + Grain -> r1(gr) + Grain)
C 1 - Cosmic ray ionization (r2 = CRP):
C 2-12 - Two-body reaction
C 13 - Photoreaction (r2 = PHOTON):
C 14 - X-ray ionization (r2 = XRAY):
C 15 - UV Photolysis induced by X-rays
C		Taken from Aikawa & Herbst, 2001, eqn 8.
C		psi = number of Lyman-Warner photons produced per H2 ionization
C		omega = grain albedo
C 20 - Depletion terms (e = 100), basically same as 0 but from Bergin's
C		gasgr.f code
C 21 - Classical evaporation terms (e = 101), from Bergin's gasgr.f code
C 22 - Cosmic-ray desorption  from Hasegawa et al. - just evaporation
C		rate at 70 K times the time spent at 70 K when hit by cosmic
C		ray.  so t_cosmic_ray = 70 K.  time spent at 70 K = frac (e = 102)
C		beta = 0 is the old e = -2 (modified rates from Bringa & Johnson 04)
C		From Bergin's gasgr.f code
C 23 - Depletion terms for ions
C > 30 - Manually calculated photodissociation rates.  Otherwise rtype = 13
C		numbers > 30 are for multiple paths.  Path # = rtype - 29
C 99 - alpha = k
C
C..............................................................................
	DOUBLE PRECISION FUNCTION calcrate(alpha, beta, gamma, rtype, r1,
     +								   r2,grnfrac,timedep)
	IMPLICIT NONE

C Common Blocks
c	INCLUDE "Fcn.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

C Input variables
	DOUBLE PRECISION alpha, beta, gamma
	INTEGER rtype, timedep, i
	CHARACTER*13 r1,r2, strg1,strg2
	REAL grnfrac

C Local variables
	DOUBLE PRECISION xrayion_rate, radius
	DOUBLE PRECISION vel, stick, nu0
	DOUBLE PRECISION psi, xh2, omega, etau, maxlam
	DOUBLE PRECISION tempval, Cfac
	DOUBLE PRECISION a_gr, Mlayer, ndens
	DOUBLE PRECISION heliorad
	DOUBLE PRECISION mu_i_j, kappa_ij
	REAL mass_i,mass_j	
	REAL masscalculator
	PARAMETER (heliorad = 0.0)  ! this should also be an external input LIC

C External Functions
	DOUBLE PRECISION self_shield_CO, self_shield_H2
	double precision grain_abun_check, photorate_nocalc
	DOUBLE PRECISION photorate, photodesorp_rate

C   for type 41
	REAL Nsites, nsfsites, a_spacing, v_zeri, v_zerj, t_hop_i
	DOUBLE PRECISION t_hop_j, R_diff_i, R_diff_j 
	Double precision R_diff_i_q, ndustdens,t_quant
	double PRECISION calcrateo,calcrateq

c Get required info
	radius = Rs(zone)
	ndens = rho(zone) / mH

c Initialize calcrate to 0
	calcrate = 0

C..............................................................................
C Time-dependent reactions
C..............................................................................
	IF (timedep .EQ. 1) THEN
C -1) CO Self Shielding
		IF (rtype .EQ. -1) THEN
			calcrate = self_shield_CO(alpha)
C -2) H2 Self Shielding
		ELSE IF (rtype .EQ. -2) THEN
			calcrate = self_shield_H2()
C -3) Grain reactions / Surface Chemistry
		ELSE IF (rtype .EQ. -3) THEN
			calcrate = grain_abun_check(alpha, beta, r1)
C -23) Photodesorption (e = 103), from Bergin's gasgr.f code, r2 = PHOTON
		ELSE IF (rtype .EQ. -23) THEN
			calcrate = photodesorp_rate(r1, rtype, beta)
			IF (.NOT. photodesorp) calcrate = 0
C -24) Lyman alpha Photodesorption, r2 = LYAPHOTON
		ELSE IF (rtype .EQ. -24) THEN
			calcrate = photodesorp_rate(r1, rtype, beta)
			IF (.NOT. LyAphotodesorp) calcrate = 0
		END IF

		RETURN
	END IF

C..............................................................................
C Time-independent reactions
C..............................................................................
C 0) Grain Reaction (r2 = GRAIN*):
	IF (rtype.EQ.0) THEN
		calcrate = alpha*(Tg(zone)/300.0D0)**beta

C 1) Cosmic ray ionization (r2 = CRP):
	ELSE IF (rtype.EQ.1) THEN
		calcrate = alpha * (zetaCR(zone)/1.3D-17) * CRatten(zone)
		if (.NOT. CRionization) calcrate = 0
		if (Rs(zone) .LE. heliorad) calcrate = 0

C 2-12) Two-body reaction:
	ELSE IF (rtype .GE. 2 .AND. rtype .LE. 12) THEN
		calcrate = alpha*(Tg(zone)/300.0D0)**beta*dexp(-gamma/Tg(zone))
		
C 13) Photoreaction (r2 = PHOTON):
	ELSE IF (rtype.EQ.13) THEN
c  changed to I/Imax at 1500 Ang
		calcrate = photorate_nocalc(r1, alpha)

C 14) X-ray ionization (r2 = XRAY):
	ELSE IF (rtype.EQ.14) THEN
		calcrate = xrayrate(zone)
		print *, 'X-ray ionization rate at zone ', zone, ': ', calcrate

C 15) UV Photolysis induced by X-rays
C	  Taken from Aikawa & Herbst, 2001, eqn 8.
C	  psi = number of Lyman-Warner photons produced per H2 ionization
C     omega = grain albedo
	ELSE IF (rtype .EQ. 15) THEN
		psi = 1.4
		xh2 = 0.5
		omega = 0.5
		calcrate = 2.6 * psi * alpha * xh2 * xrayrate(zone) * 
     &			1.0/(1-omega)
C 20) Depletion terms (e = 100), basically same as 0 but from Bergin's
C	  gasgr.f code
	ELSE IF (rtype.EQ.20) THEN
		vel = ((8.0 * kbol * Tg(zone))/(pi * beta * mH))**0.5
		stick = gamma
		calcrate = alpha * vel * stick

C 21) Classical evaporation terms / thermal desorption (e = 101), from Bergin's gasgr.f code
	ELSE IF (rtype.EQ.21) THEN
		nu0 = ((2.0 * gamma * kbol * alpha)/(pi*pi*mH*beta))**0.5
		calcrate = nu0 * dexp(-alpha / Td(zone))
		IF (.NOT. thermaldesorp) calcrate = 0
		
C 22) Cosmic-ray desorption  from Hasegawa et al. - just evaporation
C     rate at 70 K times the time spent at 70 K when hit by cosmic
C     ray.  so Tcr = 70 K.  time spent at 70 K = frac (e = 102)
C     beta = 0 is the old e = -2 (modified rates from Bringa & Johnson 04)
C	  From Bergin's gasgr.f code
	ELSE IF (rtype.EQ.22) THEN
		IF (beta.EQ.0) THEN
			calcrate = (alpha + beta * Td(zone) ** (-1 * gamma)) * 
     &					CRatten(zone)
		ELSE
			nu0 = ((2.0 * gamma * kbol * alpha)/(pi*pi*mH*beta))**0.5
			calcrate = nu0 * frac * dexp(-alpha / Tcr) * CRatten(zone)
		END IF
		IF (.NOT. CRdesorp) calcrate = 0
		if (Rs(zone) .LE. heliorad) calcrate = 0
		
C 23) Depletion terms for ions.  From Bergin's gasgr.f code (mostly)
C	  Cfac = factor which takes into account the charge of the accreting particle.  
C	  Eqn. 9 from Willacy et al. 1998 (a = 1000 A)
	ELSE IF (rtype.EQ.23) THEN
		vel = ((8.0 * kbol * Tg(zone))/(pi * beta * mH))**0.5
		stick = gamma
		Cfac = 1 + 16.71D-4 / (1000D-8 * Tg(zone))
		calcrate = Cfac * alpha * vel * stick

C 30) Manually calculated photodissociation rates.  Otherwise rtype = 13
C     numbers > 30 are for multiple paths.  Path # = rtype - 29
	ELSE IF (rtype.GE.30 .AND. rtype.LT.40) THEN
C		calcrate = photorate_nocalc(r1, alpha)
		calcrate = photorate(r1, rtype)
		
C 41) Adding in grain surface reactions. Grain properties need to be calculated.
C     Give E_des and E_act. 
C     alpha = E_desi
C     beta = E_desj
C	  gamma = E_act

 	ELSE IF (rtype.EQ.41) THEN
                write(*,*) 'Type 41'
		Nsites = 1E6  ! This should really be an input parameter, also used for uvfield.f
		nsfsites = 1.5e15 ! cm^-2 (TA)
		a_spacing = 1E-8 ! barrier between adjacent sitesf
		
		strg1=trim(r1)
		strg1=adjustl(strg1)
		strg2 = trim(r2)
		strg2=adjustl(strg2)
				
		! Take grain off the end
		mass_i = masscalculator(strg1)
		mass_j = masscalculator(strg2)
		
		!write(*,*) '-----------------------------'
		!write(*,*) '1 ',strg1, '2 ',strg2,'1 ', r1,'2 ', r2

		!write(*,*) mass_i, mass_j
		!write(*,*) '-----------------------------'
		
		mu_i_j = mass_i*mass_j/(mass_i+mass_j)*mH
	    kappa_ij = exp(-2*(a_spacing/(hp/(2*pi)))*   
     &			 (2*mu_i_j*kbol*gamma)**0.5)
			 
		v_zeri = 1.0/((2*nsfsites*kbol*alpha/(pi*pi*mH*mass_i))**0.5)
		v_zerj = 1.0/((2*nsfsites*kbol*beta/(pi*pi*mH*mass_j))**0.5)

		! E_b = 0.3*E_d, HHL
		t_hop_i  = v_zeri*exp(alpha*0.3/Td(zone))
		t_hop_j  = v_zerj*exp(beta*0.3/Td(zone))
		R_diff_i = 1.0/(Nsites*t_hop_i)
		R_diff_j = 1.0/(Nsites*t_hop_j)
        ndustdens = grnfrac*ndens*ndust
		!write(*,*) ndustdens
		calcrateq = 0.0
		! Quantum tunneling
c		write(*,*) 'kappa: ',kappa_ij,'R_diff_i: ',R_diff_i, R_diff_j
		select case(strg1)
         	case('H(gr)')
			!case('jn')
			t_quant = v_zeri*exp((2*a_spacing/(hp/(2*pi)))*
     &			 (2*mass_i*mH*kbol*alpha*0.3)**0.5)
			R_diff_i_q = 1.0/(Nsites*t_quant)
			if (strg2 .eq. 'H(gr)') then 
			calcrateq = kappa_ij*(2.0*R_diff_i_q)/ndustdens
			else
			calcrateq = kappa_ij*(R_diff_i_q)/ndustdens 
			endif

c			write(*,*) 'kappa: ',kappa_ij,'R_diff_i_q: ',R_diff_i_q,
c     &			 ' Tdust',Td(zone)
c			case('H2(gr)')
			t_quant = v_zeri*exp((2*a_spacing/(hp/(2*pi)))*
     &			 (2*mass_i*mH*kbol*alpha*0.3)**0.5)
			R_diff_i_q = 1.0/(Nsites*t_quant)
			calcrateq = kappa_ij*(R_diff_i_q)/ndustdens			
         	case default
			calcrateo = kappa_ij*(R_diff_i+R_diff_j)/ndustdens	
        end select
        calcrate = max(calcrateo,calcrateq)
		write(*,*) 'grain rate: ',calcrate
		
C 42) Grain surface reactions. Temperature dependent barrier (H + CO...)
C     Give E_des and E_act. 
C     alpha = E_desi
C     beta = E_desj
C	  gamma = E_act

 	ELSE IF (rtype.EQ.42) THEN
		Nsites = 1E6  ! This should really be an input parameter, also used for uvfield.f
		nsfsites = 1.5e15 ! cm^-2 (TA)
		a_spacing = 1E-8 ! barrier between adjacent sitesf
		
		strg1=trim(r1)
		strg1=adjustl(strg1)
		strg2 = trim(r2)
		strg2=adjustl(strg2)
				
		! Take grain off the end
		mass_i = masscalculator(strg1)
		mass_j = masscalculator(strg2)
		
		!write(*,*) '-----------------------------'
		!write(*,*) '1 ',strg1, '2 ',strg2,'1 ', r1,'2 ', r2

		!write(*,*) mass_i, mass_j
		!write(*,*) '-----------------------------'
		
		mu_i_j = mass_i*mass_j/(mass_i+mass_j)*mH
	    kappa_ij = exp(-2*(a_spacing/(hp/(2*pi)))*
     &			(2*mu_i_j*gamma)**0.5)  ! 
		
		v_zeri = 1.0/((2*nsfsites*alpha/(pi*pi*mH*mass_i))**0.5)
		v_zerj = 1.0/((2*nsfsites*beta/(pi*pi*mH*mass_j))**0.5)
		! Just include hopping from the lighter reactant. Make sure*** lighter reactant
		! always first?  Could add a flag here with the mass_i,j calculations above to be sure.
		! E_b = 0.3*E_d, HHL
		t_hop_i  = v_zeri*exp(alpha*0.3/(kbol*Td(zone)))
		t_hop_j  = v_zerj*exp(beta*0.3/(kbol*Td(zone)))
		R_diff_i = 1.0/(Nsites*t_hop_i)
		R_diff_j = 1.0/(Nsites*t_hop_j)
        ndustdens = grnfrac*ndens*ndust
		!write(*,*) ndustdens
		
		! Quantum tunneling
		
		select case(strg1)
         	case('H(gr)')
			t_quant = v_zeri*exp((2*a_spacing/(hp/(2*pi)))*
     &			(2*mass_i*alpha*0.3)**0.5)
			R_diff_i_q = 1.0/(Nsites*t_quant)
			calcrate = kappa_ij*(R_diff_i_q)/ndustdens
			case('H2(gr)')
			t_quant = v_zeri*exp((2*a_spacing/(hp/(2*pi)))*
     &			(2*mass_i*alpha*0.3)**0.5)
			R_diff_i_q = 1.0/(Nsites*t_quant)
			calcrate = kappa_ij*(R_diff_i_q)/ndustdens			
         	case default
			calcrate = kappa_ij*(R_diff_i+R_diff_j)/ndustdens	
        end select

C 99) alpha = k, no other calculations needed
	ELSE IF (rtype.EQ.99) THEN
		calcrate = alpha

C Ignore Time-dependent Reaction for now
 	ELSE IF (rtype .LT. 0) THEN
		calcrate = 0
		
	ELSE
		write(*,*) 'ERROR: No Rate calculation for given reaction'
		write(*,*) 'rtype = ', rtype
		stop

	END IF 

	RETURN
	END

C..............................................................................
C  
C CR Attenuation
C
C Calculate the attenuation factor for cosmic rays, from Bergin PPV paper.
C
C..............................................................................
	SUBROUTINE calcCRatten()
	
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	
	INTEGER i
	DOUBLE PRECISION colden_tot, colden1(Nz), colden2(Nz)
	DOUBLE PRECISION sigmaCR
	PARAMETER (sigmaCR = 96)
		
	colden_tot = 0.0
	DO i=1,Nz
		if (i .eq. 1) then
			colden_tot = colden_tot + rho(i) * 2 * zcm(i)
		else
			colden_tot = colden_tot + rho(i) * 2 * (zcm(i) - zcm(i-1))
		endif
	END DO
	 
	DO i=1,Nz
		if (i .eq. 1) then
			colden1(i) = rho(i) * 2 * zcm(i)
		else
			colden1(i) = colden1(i-1) + rho(i) * 2 * (zcm(i) - zcm(i-1))
		end if


		colden2(i) = 2 * colden_tot - colden1(i)
        write(*,*) 'colden1: ',colden1(i),'colden2: ',colden2(i)
		CRatten(i) = 0.5 * (exp(-1 * colden1(i)/sigmaCR) + 
     &				 exp(-1 * colden2(i)/sigmaCR))
     
		print *, 'CR Atten for zone ',i,' = ', CRatten(i)
	END DO
	
	RETURN
	END


C..............................................................................
C
C X-Ray Ionization
C
C Rate calculated from Aikawa & Herbst 1999 (1999A&A...351..233A) eqn 5 
C and Glassgold et al. 1997
C
C Assumptions:
C kTxr = 2 keV
C E1 = 0.8 keV, E2 = 20 keV
C Nsec (number of secondary ionizations per unit photoelectron energy) = 30
C E range 1 - 30 keV
C..............................................................................
	DOUBLE PRECISION FUNCTION xrayion_rate()
	IMPLICIT NONE

	DOUBLE PRECISION xrayinteg
	EXTERNAL xrayinteg
	DOUBLE PRECISION kTmin, kTmax
	PARAMETER (kTmin = 1.0, kTmax = 20.0)
		
	CALL qromb(xrayinteg, ktmin, kTmax, xrayion_rate)
	RETURN
	END

C..............................................................................
C Integrand for X-ray ionization rate calculation
C
C Assumptions:
C kTxr = 2 keV
C E1 = 0.8 keV, E2 = 10 keV
C Nsec (number of secondary ionizations per unit photoelectron energy) = 30
C E range 1 - 10 keV
C..............................................................................
	DOUBLE PRECISION FUNCTION xrayinteg(E)
	
C Common Blocks
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"
	
	DOUBLE PRECISION E

	DOUBLE PRECISION Lv0, Lv, F0, F
	DOUBLE PRECISION KTXR, E1, E2, Nsec, c
	PARAMETER (KTXR = 2, E1 = 0.8, E2 = 10.0, Nsec = 30)
	PARAMETER (c = 2.99792458D10)

C External Functions
	DOUBLE PRECISION sigmaE, Habssigma, xray_flux

	
c	Lv0 = Lxray/KTXR * (dexp(-E1/KTXR) - dexp(-E2/KTXR))
c	Lv = Lv0 * dexp(-E/KTXR)
c	F0 = Lv/(8*pi*Rs(zone)**2*AU**2 * E * keV)
c	xrayinteg = Nsec * Habssigma(E) * F0 * dexp(-sigmaE(E)*Nrz(zone))

c	xrayinteg = Nsec * Habssigma(E) * nxray_photons(zone) * c

	xrayinteg = Nsec * Habssigma(E) * xray_flux(E)

	if (xrayinteg .LE. 1.0e-30) then
		xrayinteg = 1.0e-30
	endif

c	print *, E, xrayinteg, sigmaE(E), Nrz(zone)

	RETURN 
	END

C..............................................................................
C
C SigmaE
C
C Calculates the sigma value for a given E, using the work by Bethell
C 2010 in prep  
C
C Assume fb = self-blanketing = 1
C For now, assume no dust settling dependence for sigma
C sigmatot = sigmagas + eps * fb * sigmadust
C
C..............................................................................
	DOUBLE PRECISION FUNCTION sigmaE(Einp)
	
	INCLUDE "environ.h"
	
	DOUBLE PRECISION Evals(14), Einp
	DOUBLE PRECISION cgas(14,3), cdust(14,3)
	DOUBLE PRECISION fb, sigmagas, sigmadust
	PARAMETER(fb = 1.0)
	INTEGER i, j

	Evals(1) = 0.030
	Evals(2) = 0.100
	Evals(3) = 0.284
	Evals(4) = 0.400
	Evals(5) = 0.532
	Evals(6) = 0.708
	Evals(7) = 0.867
	Evals(8) = 1.303
	Evals(9) = 1.840
	Evals(10) = 2.471
	Evals(11) = 3.210
	Evals(12) = 4.038
	Evals(13) = 7.111
	Evals(14) = 8.331
	
	cgas = reshape ( (/ 17.9, 37.0, 49.3, 58.4, 49.2, 78.0,  
     &				79.8, 117.0, 109.0, 108.0, 138.0, 142.0, 139.0, 95.9,
     &				560.0, 175.0, 84.4, 38.0, 127.0, 44.7, 70.1, 7.52, 
     &				14.5, 12.7, -2.16, -4.71, -3.59, 6.61, 
     &				-2320.0, -309.0, -100.0, -41.0, -79.2, -21.0, -28.3, 
     &				-1.89, -3.4, -2.49, -0.153, 0.24, 0.147, 
     &				-0.461 /), shape(cgas))

	cdust = reshape( (/ -0.0359, -1.77, -0.507, -3.63, -6.53, -73.5, 
     &						9.67, -14.2, 30.4, 28.5, 118.0, 186.0, 1040.0,31.0,
     &					0.873, 27.9, 29.0, 55.8, 84.2, 261.0, 71.8, 113.0, 
     &						76.3, 82.0, 27.0, -1.06, -133.0, 123.0, 
     &					71.4, -22.8, -8.1, -38.1, -43.8, -132.0, -20.5, 
     &						-27.8, -12.3, -10.8, -2.41, 0.918, 10.2, 
     &						-5.61 /), shape(cdust))


	DO i = 1,14 
		IF (Einp .LT. Evals(1)) THEN
			sigmaE = 0.0
			RETURN
		ELSE IF (Einp .LT. Evals(i)) THEN
			IF (Einp .GE. Evals(14)) THEN
				j = 14
			ELSE
				j = i-1
			END IF
			EXIT
		ENDIF
	END DO

		
	sigmagas = (cgas(j,1) + cgas(j,2) * Einp + 
     &				cgas(j,3) * Einp**2)
	sigmadust = (cdust(j,1) + cdust(j,2) * 
     &				Einp + cdust(j,3) * Einp**2)

	IF (xraydust) THEN
		sigmadust = sigmadust * fb * ndust
	ELSE
		sigmadust = sigmadust * fb * 1.0
	ENDIF

	sigmaE = (sigmagas + sigmadust) * 1d-24 * Einp**(-3)

	RETURN
	END

C..............................................................................
C
C HabsSigma
C
C Calculates the sigma value for a given E, using the data from 
C http://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
C for photoelectric absorption of H 
C
C values here in 1e-24 cm**2, so multiply the values from that website by
C 1.6733 to get the values below
C
C..............................................................................
	DOUBLE PRECISION FUNCTION Habssigma(Einp)
	
	INTEGER nvals
	PARAMETER (nvals=11)
	DOUBLE PRECISION sigmavals(nvals), Evals(nvals), Einp
	INTEGER i, j
	
	Evals(1) = 1.0
	Evals(2) = 1.5
	Evals(3) = 2.0
	Evals(4) = 3.0
	Evals(5) = 4.0
	Evals(6) = 5.0
	Evals(7) = 6.0
	Evals(8) = 8.0
	Evals(9) = 10.0
	Evals(10) = 15.0
	Evals(11) = 20.0
	
	sigmavals(1) = 1.14D+01
	sigmavals(2) = 2.93D+00
	sigmavals(3) = 1.11D+00
	sigmavals(4) = 2.81D-01
	sigmavals(5) = 1.05D-01
	sigmavals(6) = 4.91D-02
	sigmavals(7) = 2.63D-02
	sigmavals(8) = 9.82D-03
	sigmavals(9) = 4.56D-03
	sigmavals(10) = 1.13D-03
	sigmavals(11) = 4.18D-04
	
	
	DO i=1,nvals
		IF (Evals(i) .EQ. Einp) THEN
			Habssigma = 1D-24 * sigmavals(i)
			RETURN
		ELSE IF (Einp < Evals(i)) THEN
			j = i-1
			EXIT
		ENDIF
	ENDDO

	Habssigma = ((Einp-Evals(j))*sigmavals(j+1) + 
     &			 (Evals(j+1)-Einp)*sigmavals(j)) / 
     &			 (Evals(j+1)-Evals(j))
	
	Habssigma = Habssigma * 1d-24

	RETURN
	END
	
C..............................................................................
C
C QROMB
C
C Returns the integral of the function func from a to b.  
C Integration is performed by Romberg's method of order 2K, where 
C K = 2 is Simpson's rule.
C Taken from Numerical Recipes
C..............................................................................
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=25, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd2(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss) .or. abs(dss).le.1.0e-25) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
	  write(*,*) 'ERROR: too many steps in qromb'
	  STOP
      END

C..............................................................................
C
C TRAPZD2
C
C trapzd routine taken from NR
C "Routine implementing the extended trapezoidal rule"
C Version in uv_field.f is modified, hence the second copy.
C..............................................................................
      SUBROUTINE trapzd2(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END

C..............................................................................
C
C Polint
C
C polint routine taken from NR
C..............................................................................
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
			write(*,*) 'failure in polint'
			stop
		  end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END


C..............................................................................
C
C Read Bethell X-ray File
C
C Read in scattering photon values from Tom Bethell's file (filename given as 
C input parameter in 0io file)
C..............................................................................
	SUBROUTINE read_bethell_xrayscatter(xrayfile)
	IMPLICIT NONE
	
	INCLUDE "environ.h"

	CHARACTER*80 xrayfile
	INTEGER fend, i
	LOGICAL foundrad
	DOUBLE PRECISION Rfile(nzone), zfile(nzone), nphfile(nzone)
	DOUBLE PRECISION Rtmp, ztmp, nphtmp
	DATA Rfile/nzone*0.0/, zfile/nzone*0.0/, nphfile/nzone*0.0/

C Read in xray scattering data -> store in nxray_photons(zone)
	OPEN(unit=01, file=trim(xrayfile))
	DO i=1,7
		read(01,*)
	ENDDO
	
	fend = 0
	foundrad = .FALSE.
	DO WHILE (fend .EQ. 0)
		read(01, *, iostat=fend) Rtmp, ztmp, nphtmp
		if (fend .NE. 0) exit

		IF ( Rtmp .NE. Rfile(1) ) THEN
			if ( foundrad ) exit

			i = 1

			IF ( Rtmp .GT. Rs(zone) ) THEN
				if ( (Rtmp - Rs(zone)) .LT. (Rs(zone) - Rfile(1)) ) then
					foundrad = .TRUE.
				else
					exit
				endif
			ENDIF
		ENDIF

		Rfile(i) = Rtmp
		zfile(i) = ztmp
		nphfile(i) = nphtmp
		i = i + 1
	ENDDO

	CLOSE(01)
	
c	DO i = i,Nr*Nz
c		nxray_photons(i) = nphfile(i)
c	ENDDO

	RETURN
	END
	
C..............................................................................
C
C masscalculator
C lic 1/14
C Calculate mass of a given molecule
C..............................................................................

	REAL FUNCTION masscalculator(molstring)

	implicit none
	CHARACTER*10 molstring, d,dn,f, fis
	DOUBLE PRECISION atomarr(5), masses(5)
	INTEGER ct, x,xold, oldatomindex, numberofatoms, oldnumberofatoms
	DOUBLE PRECISION C_n, O_n, N_n, He_n, H_n, Mg_n, Na_n, Ca_n, S_n
	DOUBLE PRECISION massval, valnum
	CHARACTER*15 is_digit
	character*80 formt
	LOGICAL ios, res, wasdigit
	CHARACTER*20 str,fmt
	real rnum

	H_n = 1.0
	C_n = 12.0
	O_n = 16.0
	N_n = 14.0
	He_n = 4.0
	H_n = 1.0
	Mg_n = 24.31
	Na_n = 23.0
	Ca_n = 40.08
	S_n = 32.0
	
	masses = [H_n, C_n, O_n, N_n, S_n]
	
	!write(*,*) 'function sees: ',molstring
	d = 'wow'
	wasdigit = .false.

c	          [H,C,O,N,S]
	atomarr = [0,0,0,0,0]
	ct = 1
	
	do while ((d .ne. '').and.(d .ne. '('))
	     d = molstring(ct:ct)
         select case(d)
         	case('0':'9')
			if (wasdigit) then
			    xold = x
			    read( d, '(i10)' )  x
                formt='('//trim('i1,i1')//')'
                write(str,formt) xold,x
                str=adjustl(str)
				read(str, '(i20)') numberofatoms				
			    atomarr(oldatomindex) = atomarr(oldatomindex)-
     &				 xold+numberofatoms
		        oldnumberofatoms = numberofatoms
				wasdigit = .true.			    
			else
			    read( d, '(i10)' )  x
			    numberofatoms = x
			    atomarr(oldatomindex) = atomarr(oldatomindex) 
     &				 + numberofatoms-1
		        oldnumberofatoms = numberofatoms
				wasdigit = .true.
			endif
			
         	case default		
		    wasdigit = .false.  
	          SELECT CASE (d)       ! number of type integer
	          CASE ('H')                 ! all values below 0
				   dn = molstring(ct+1:ct+1)
				   if (dn .eq. 'e') then 
				   !write(*,*) 'He!'
				   ct = ct + 1
			       else
	               atomarr(1) = atomarr(1) + 1
				   oldatomindex = 1
				   endif
				   
	          CASE ('C')                   
	               atomarr(2) = atomarr(2) + 1
				   oldatomindex = 2
	          CASE ('O')                  
	               atomarr(3) = atomarr(3) + 1
				   oldatomindex = 3
			  CASE ('N')                  
	               atomarr(4) = atomarr(4) + 1
				   oldatomindex = 4   
			  CASE ('S')                  
				   dn = molstring(ct+1:ct+1)
				   if (dn .eq. 'i') then 
				   !write(*,*) 'Si!'
				   ct = ct + 1	
			       else
	               atomarr(5) = atomarr(5) + 1
				   oldatomindex = 5 
				   endif
				   	
			  !CASE ('+')                  
			  ! write(*,*) d		
			  !CASE ('E')                    	   
		      CASE DEFAULT
			 
	          END SELECT			
		 end select 
c		 write(*,*) d
		 ct = ct + 1
		 
	enddo

	masscalculator = sum(atomarr*masses)

	return
	end




CC    SUBROUTINE dustdensity(numdens,numdust)
C  C  INCLUDE "Fcn.h"
C	CINCLUDE "environ.h"
C	INCLUDE "constants.h"
C	INCLUDE "rates.h"
	
C	real numdens, numdust

C   numdust = frc(12)*numdens
C  RETURN
C    END
