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
     +								   r2,grnfrac,timedep,dtaudt,timeval)
	IMPLICIT NONE

C Common Blocks
c	INCLUDE "Fcn.h"
	INCLUDE "BindingE.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

C Input variables
	DOUBLE PRECISION alpha, beta, gamma
	DOUBLE PRECISION dtaudt(64,3), temperaturevec(64), dervec(64),v1_n(63)
	DOUBLE PRECISION v2_n(63),numer_v1,denom_v2,del_h(62)
	DOUBLE PRECISION timeval
	INTEGER rtype, timedep, i, qn
	CHARACTER*13 r1,r2, strg1,strg2
	REAL grnfrac

C Local variables
	DOUBLE PRECISION xrayion_rate, radius
	DOUBLE PRECISION xrayionHe_rate
	DOUBLE PRECISION vel, stick, nu0
	DOUBLE PRECISION psi, xh2, omega, etau, maxlam
	DOUBLE PRECISION tempval, Cfac
	DOUBLE PRECISION a_gr, Mlayer, ndens
	DOUBLE PRECISION heliorad, crayreduc
	DOUBLE PRECISION mu_i_j, kappa_ij
	DOUBLE PRECISION grnfoc
	DOUBLE PRECISION stickingco,sigmagr,rgr
	REAL mass_i,mass_j
	REAL masscalculator

C External Functions
	DOUBLE PRECISION self_shield_CO, self_shield_H2
	DOUBLE PRECISION self_shield_HD, self_shield_D2
	DOUBLE precision grain_abun_check, photorate_nocalc
	DOUBLE PRECISION photorate, photodesorp_rate
	DOUBLE PRECISION h2form
	DOUBLE PRECISION binding_surface
	DOUBLE PRECISION compmono
	DOUBLE PRECISION occupyice

C   for type 41
	REAL Nsites, nsfsites, a_spacing, v_zeri, v_zerj, t_hop_i
	DOUBLE PRECISION t_hop_j, R_diff_i, R_diff_j
	Double precision R_diff_i_q, ndustdens,t_quant
	double PRECISION calcrateo,calcrateq,rate_xray
	double precision time_xray_hit,time_xray_cooldown
	DOUBLE PRECISION xrayheat_rate, normrate
	Double PRECISION opratH2,f_orthoH2,f_paraH2,opratH2Dp
	Double Precision f_orthoH2Dp,f_paraH2Dp,fro
	Double Precision E_Mg,E_Bt,E_An,sig_Mg,sig_Bt,sig_An,sigH
	Double Precision X26AL,Eatten1,Eatten2,f26Al1,f26Al2
	Double Precision meangrnsz
	Double Precision bind1, bind2
	INTEGER beind1, beind2
	DOUBLE PRECISION sigappr,zetanew
	DOUBLE PRECISION physlim,chemtemp
	DOUBLE PRECISION monoly
	INTEGER r1index,r2index
	DOUBLE PRECISION occupation1,occupation2
	DOUBLE PRECISION sigadjust

	sigadjust = 1.0
        if (incl_locdust) then
          sigadjust = locdust(zone)
         ! print *, 'Local dust fraction adjusted by: ',sigadjust
        end if


c Get required info
	radius = Rs(zone)
	ndens = rho(zone) / mH

c Initialize calcrate to 0
	calcrate = 0.0

C   At high dust temperatures (set by physlim) change the H-binding energy to the chemisorbed value (chemtemp).
	physlim = 125.0       ! Td (Kelvin)
	chemtemp = 10000.0    ! Eb (Kelvin)

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
C -3) Grain reactions / Surface Chemistry  (Obsolete).
		ELSE IF (rtype .EQ. -3) THEN
			calcrate = grain_abun_check(alpha, beta, gamma, r1)
C -4) Grain reactions / Surface Chemistry of H2 with Cazaux and Tielens formalism   (Obsolete).
		ELSE IF (rtype .EQ. -4) THEN
			calcrate = h2form(alpha, beta, gamma, r1)
C -5) HD Self Shielding
		ELSE IF (rtype .EQ. -5) THEN
			calcrate = self_shield_HD()
C -6) D2 Self Shielding
		ELSE IF (rtype .EQ. -6) THEN
			calcrate = self_shield_D2()

C -21) Classical evaporation terms for H and D only.
	ELSE IF (rtype.EQ.-21) THEN

        if ((r1 .ne. 'H(gr)        ').and.((r1 .ne. 'D(gr)        '))) then
			write(*,*) "Wrong application of binding energy.",r1,rtype
			stop
		end if

        ! Check how any monolayers exist on grains to see which binding energy is appropriate
		monoly = compmono()

		! Increasing H(gr) and D(gr) B.E. to simulate energy between chemi and physisorption.
        bind1 = 0.0
		if (monoly .lt. 1.0) then ! If no ice coating, use chemisorption value.
        	bind1 = chemtemp
        else ! Otherwise read in from 6grain.be.
            CALL ispecies(r1, size(bindspec), bindspec, beind1)
        	bind1 = bindenerg(beind1)
        end if

		nu0 = ((2.0 * gamma * kbol * bind1)/(pi*pi*mH*beta))**0.5
		calcrate = nu0 * dexp(-bind1 / Td(zone))
		IF (.NOT. thermaldesorp) calcrate = 0

C -23) Photodesorption (e = 103), from Bergin's gasgr.f code, r2 = PHOTON
C      <*> Depends on grain size.
		ELSE IF (rtype .EQ. -23) THEN
			calcrate = photodesorp_rate(r1, rtype, beta)
			IF (.NOT. photodesorp) calcrate = 0

C -24) Lyman alpha Photodesorption, r2 = LYAPHOTON
C      <*> Depends on grain size.
		ELSE IF (rtype .EQ. -24) THEN
			calcrate = photodesorp_rate(r1, rtype, beta)
			IF (.NOT. LyAphotodesorp) calcrate = 0

C -25) Surface dependent Binding Energies for CO and N2 ice.
C      Calls time_dependent.f
		ELSE IF (rtype .EQ. -25) THEN
            if ((r1 .ne. 'CO(gr)       ').and.((r1 .ne. 'N2(gr)       '))) then
                write(*,*) "Wrong application of binding energy.",r1,rtype
                stop
            end if
		    calcrate = binding_surface(alpha,beta,gamma)

C -41) Grain surface chemistry reduced efficiency (restricted to surface layer)
C      <*> Depends on grain size.

		ELSE IF (rtype .EQ. -41) THEN
            calcrate = 0.0
            if (zone .gt. 3) then

                Nsites = 1E6		! Number of sites per grain.
                if (freezeeffic .ne. 1.0) then
                    sigmagr = pi*(rgr*freezeeffic**(-1.0/1.5))**2
                    Nsites = Nsites * sigmagr/(pi*rgr**2)  ! scale number of sites by increase/decrease in surface area per grain (grain growth)
                end if

                nsfsites = 1.5e15   ! number density of surface sites. cm^-2 (TA) (do not change).
                a_spacing = 1E-8    ! Barrier between adjacent sites (Ang).

                strg1=trim(r1)
                strg1=adjustl(strg1)
                strg2=trim(r2)
                strg2=adjustl(strg2)

                ! Take grain off the end
                mass_i = masscalculator(strg1)
                mass_j = masscalculator(strg2)

                mu_i_j = mass_i*mass_j/(mass_i+mass_j)*mH
                kappa_ij = exp(-2*(a_spacing/(hp/(2*pi)))*
     &			 (2*mu_i_j*kbol*gamma)**0.5)

                beind1 = -1
                beind2 = -1
                bind1 = 0.0
                bind2 = 0.0

                CALL ispecies(r1, size(bindspec), bindspec, beind1)
                CALL ispecies(r2, size(bindspec), bindspec, beind2)
                bind1 = bindenerg(beind1)
                bind2 = bindenerg(beind2)

                if (bind1 .EQ. 0 .or. bind2 .EQ. 0) then
				write(*,*) "Can't find B.E. in 6grainbe.inp -gr.",r1,r2
				stop
                end if

                v_zeri = 1.0/((2*nsfsites*kbol*bind1/(pi*pi*mH*mass_i))**0.5)
                v_zerj = 1.0/((2*nsfsites*kbol*bind2/(pi*pi*mH*mass_j))**0.5)

                ! Assuming: E_b = 0.3*E_d, Hasegawa,Herbst&Leung.
                t_hop_i  = v_zeri*exp(bind1*0.3/Td(zone))
                t_hop_j  = v_zerj*exp(bind2*0.3/Td(zone))
                R_diff_i = 1.0/(Nsites*t_hop_i)
                R_diff_j = 1.0/(Nsites*t_hop_j)

                ndustdens = grnfrac*ndens*ndust  ! Grain abundance * number H2 * epsilon
                if (freezeeffic .ne. 1.0) then
                    ndustdens = ndustdens * (freezeeffic**(3.5/1.5))  ! change in grain number density per spatial volume with growth.
                end if

c               if ((strg2 .eq. 'OD(gr)').and.(strg1 .eq. 'H(gr)')) then
c                   write(*,'(A30)') '-----------------------------'c
c                   write(*,'(F3.0,F7.2)') zone,Td(zone)
c                   write(*,'(A13,A13)') strg1,strg2
c                   write(*,'(E13.3,E13.3)') bind1,bind2
c                   write(*,'(E13.3,1x,E13.3,1x,E13.3,1x,E13.3)') kappa_ij,
c     &							mu_i_j/mh, R_diff_i, R_diff_j
c                   write(*,*) mass_i, mass_j
c                   write(*,*) t_hop_i,t_hop_j
c                   write(*,'(A30)') '-----------------------------'
c               end if

                if (ndens .gt. 1E+1) then

                    calcrateq = 0.0

                    select case(strg1)

                        case('H(gr)') ! Quantum tunneling for H
                            t_quant = v_zeri*exp((2*a_spacing/(hp/(2*pi)))
     &                          *(2*mass_i*mH*kbol*bind1*0.3)**0.5)
                            R_diff_i_q = 1.0/(Nsites*t_quant)
                            if (strg2 .eq. 'H(gr)') then
                                calcrateq=kappa_ij*(2.0*R_diff_i_q)
     &                              /ndustdens
                            else if (strg2 .eq. 'D(gr)') then
                                calcrateq=kappa_ij*(2.0*R_diff_i_q)
     &                              /ndustdens
                            else
                                calcrateq=kappa_ij*(R_diff_i_q)/ndustdens
                            endif
                            calcrateo = kappa_ij*(R_diff_i+R_diff_j)
     &							/ndustdens

                        case('D(gr)') ! Quantum tunneling for D
                            t_quant = v_zeri*exp((2*a_spacing/(hp/(2*pi)))
     &                          *(2*mass_i*mH*kbol*bind1*0.3)**0.5)
                            R_diff_i_q = 1.0/(Nsites*t_quant)
                            if (strg2 .eq. 'H(gr)') then
                                calcrateq=kappa_ij*(2.0*R_diff_i_q)
     &                              /ndustdens
                            else if (strg2 .eq. 'D(gr)') then
                                calcrateq=kappa_ij*(2.0*R_diff_i_q)
     &                              /ndustdens
                            else
                                calcrateq=kappa_ij*(R_diff_i_q)/ndustdens
                            endif
                            calcrateo = kappa_ij*(R_diff_i+R_diff_j)
     &							/ndustdens

                        case default ! Sweeping time for all others.
                            calcrateo = kappa_ij*(R_diff_i+R_diff_j)
     &							/ndustdens

					end select

					calcrate = max(calcrateo,calcrateq)
					if (calcrate .lt. 1e-30) calcrate = 0.0
                else
                    calcrate = 0.0  ! Not enough gas to compute grain chemistry.
                endif
            else
                calcrate = 0.0  ! In the upper atmosphere, no grain chem computed.
            endif


            monoly = compmono()
            if (monoly .le. 1.0) monoly = 1.0

            occupation1 = occupyice(r1)
            occupation2 = occupyice(r2)

            if ((strg2 .eq. 'H(gr)').and.(strg1 .eq. 'H(gr)')) monoly = 1.0
            if ((strg2 .eq. 'H(gr)').and.(strg1 .eq. 'D(gr)')) monoly = 1.0
            if ((strg2 .eq. 'D(gr)').and.(strg1 .eq. 'H(gr)')) monoly = 1.0
            if ((strg2 .eq. 'D(gr)').and.(strg1 .eq. 'D(gr)')) monoly = 1.0

c           if ((strg2 .eq. 'H(gr)').and.(strg1 .eq. 'H(gr)')) print *,strg1,strg2,
c     &         calcrateq,calcrateo,calcrate,monoly,occupation1,occupation2
C           print *, r1,r2,calcrate
C           print *,strg1,strg2,monoly,occupation1,occupation2,
C     &         calcrateq,calcrateo,calcrate

            if (occupation1 .gt. 1.0) occupation1 = 1.0
            if (occupation2 .gt. 1.0) occupation2 = 1.0

            calcrate = calcrate * alpha * 1.0/monoly *
     &          occupation1 * occupation2

            if (calcrate .lt. 1e-30) calcrate = 0.0

C -45) X-ray thermal desorption (not implemented)
C      using 21) Classical evaporation terms / thermal desorption (e = 101), from Bergin's gasgr.f code
		ELSE IF (rtype.EQ.-45) THEN

			IF ((Td(zone) .lt. 140).and.(zone.gt.2)) then
				nu0 = ((2.0 * gamma * kbol * alpha)/(pi*pi*mH*beta))**0.5
c				print *, 'n0', nu0,alpha, kbol,gamma
				temperaturevec = dtaudt(:,1)
				dervec = dtaudt(:,2)
				v1_n = (nu0*dexp(-alpha/temperaturevec(1:63)))*dervec(1:63)
				normrate = (nu0*dexp(-alpha/temperaturevec(63)))
				v2_n = dervec(1:63)
				numer_v1 = 0
				denom_v2 = 0
				del_h = temperaturevec(2:63)-temperaturevec(1:62)

c				Integrate eq 13 of Hase. and Herb. via trapezoidal int. over the cool-
c				down period.  Then weight by time spent in cool down vs. time spent
c				in quiescent x-ray heating-less period.

				do qn = 1,62
					numer_v1 = numer_v1+0.5*(v1_n(qn)+v1_n(qn+1))*del_h(qn)
					denom_v2 = denom_v2+0.5*(v2_n(qn)+v2_n(qn+1))*del_h(qn)
				enddo

c				Photons/second from rate_xray:
				rate_xray = xrayheat_rate()
c				If the time between x-ray photons is longer than the hubble time set rate to zero.
				if (1.0/rate_xray .gt. 4e17) then
					calcrate = 0.0
				else
c					write(*,*) 'vals ind: ',numer_v1,denom_v2
C					Length of time between subsequent x-ray hits:
c					Seconds per photon:
					time_xray_hit = 1.0/rate_xray
					time_xray_cooldown = dtaudt(63,3)
c					write(*,*) 'vals: ',numer_v1/denom_v2, time_xray_cooldown/time_xray_hit
					calcrate = numer_v1/denom_v2*time_xray_cooldown/time_xray_hit
					write(*,*) 'th:',time_xray_hit,'tc:',time_xray_cooldown
					write(*,*) 'x-ray calcrate:',calcrate
					write(*,*) 'thermal calcrate:',normrate
				endif
            ELSE
				calcrate = 0.0
c				calcrate = 3.210335712098724E-002
			ENDIF

			if (calcrate .lt. 1e-30) calcrate = 0.0

C -48) Radionuclide ionization (scaling taken from cosmic ray ionization rates). with time decay
C     Turn off/on in 5flags.inp
C     LIC 8/4/2013.
C         ->  Alpha is the difference in ionization rates between H2 and
C             species X (He is 0.84, for example -- Umebayashi et al. 2013).
C         ->  Beta is the 26 Al abundance.

		ELSE IF (rtype.EQ.-48) THEN

	    	calcrate = 0.0
	    	IF (incl_radionuc .eqv. .TRUE.) then
C				print *,'Time Decay of RN:',0.5**(timeval/year/1e6*1.04)
				calcrate = alpha*RNrate(zone) * 0.5**(timeval/year/1e6*0.84)
	    	END IF


		END IF

        if (calcrate .lt. 1e-30) calcrate = 0.0

		RETURN

	END IF

C..............................................................................
C Time-independent reactions
C..............................................................................
C 0) Grain Reaction (r2 = GRAIN*) of E+GRAIN0.
C      <*> Depends on grain size.
	IF (rtype.EQ.0) THEN
		calcrate = freezeeffic*ndust*alpha*(Tg(zone)/300.0D0)**beta

C 1) Cosmic ray ionization (r2 = CRP):
	ELSE IF (rtype.EQ.1) THEN
		calcrate = alpha * (zetaCR(zone)/1.3D-17) * CRatten(zone)
		if (.NOT. CRionization) calcrate = 0

C 2-12) Two-body reaction:
	ELSE IF (rtype .GE. 2 .AND. rtype .LE. 12) THEN
		calcrate = alpha*(Tg(zone)/300.0D0)**beta*dexp(-gamma/Tg(zone))

C 13) Photoreaction (r2 = PHOTON):
C     UV rates. Considers lam < 1500 A in the integrated UV field.
	ELSE IF (rtype.EQ.13) THEN
		calcrate = photorate_nocalc(r1, alpha)

C 14) X-ray ionization (r2 = XRAY)
	ELSE IF (rtype.EQ.14) THEN
		calcrate = xrayrate(zone)
		write (*,'(A13,A30,I3,A2,E10.3)') r1,' X-ray ionization rate at zone ',
     &		zone,': ',calcrate
		if (calcrate .lt. 1e-30) calcrate = 0.0

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
        if (calcrate .lt. 1e-30) calcrate = 0.0

C 20) Depletion terms (e = 100), basically same as 0 but from EAB's gasgr.f code
C      <*> Depends on grain size.

	ELSE IF (rtype.EQ.20) THEN
		vel = ((8.0 * kbol * Tg(zone))/(pi * beta * mH))**0.5
		stick = gamma
C		calcrate = alpha * vel * stick * freezeeffic * ndust
                calcrate = alpha * vel * stick * freezeeffic * ndust * sigadjust

C 21) Classical evaporation terms / thermal desorption (e = 101), from EAB's gasgr.f code
	ELSE IF (rtype.EQ.21) THEN
		bind1 = 0.0
		CALL ispecies(r1, size(bindspec), bindspec, beind1)
		bind1 = bindenerg(beind1)  ! look up binding energy from input table.
		if (bind1 .EQ. 0 .or. beind1 .eq. 0) then
			write(*,*) "Can't find B.E. in 6grainbe.inp -Therm.",r1
			stop
		end if

		! Increasing H(gr) and D(gr) B.E. to simulate energy between chemi and physisorption.
        ! This can also be treated as time dependent for reac type -21.
		if ((r1 == 'H(gr)        ').and.(Td(zone).GE.physlim)) bind1 = chemtemp
		if ((r1 == 'D(gr)        ').and.(Td(zone).ge.physlim)) bind1 = chemtemp

		nu0 = ((2.0 * gamma * kbol * bind1)/(pi*pi*mH*beta))**0.5
		calcrate = nu0 * dexp(-bind1 / Td(zone))

		IF (.NOT. thermaldesorp) calcrate = 0

C 22) Cosmic-ray desorption  from Hasegawa et al. - just evaporation
C     rate at 70 K times the time spent at 70 K when hit by cosmic
C     ray.  so Tcr = 70 K.  time spent at 70 K = frac (e = 102)
C     beta = 0 is the old e = -2 (modified rates from Bringa & Johnson 04)
C	  From Bergin's gasgr.f code
C     LIC:  Rates computed from the original Morfill et al. 1976 spectrum,
C     which has an equivalent H2 ionization rate of 9.66e-18 s^-1, using
C     the CR ionization cross sections + secondary ionization formalism of
C     Padovani et al 2009.
C
C     beta != 0.0 is the Bringa & Johnson 2004 Work
C     beta == 0.0 appears never to be used.

	ELSE IF (rtype.EQ.22) THEN

		bind1 = 0.0
		CALL ispecies(r1, size(bindspec), bindspec, beind1)
		bind1 = bindenerg(beind1)
		if (bind1 .EQ. 0 .or. beind1 .eq. 0) then
			write(*,*) "Can't find B.E. in 6grainbe.inp -CR.",r1
			stop
		end if
		 ! Increasing H(gr) and D(gr) B.E. to simulate energy between chemi and physisorption.
		 if ((r1 == 'H(gr)        ').and.(Td(zone).GE.physlim)) bind1 = chemtemp
		 if ((r1 == 'D(gr)        ').and.(Td(zone).ge.physlim)) bind1 = chemtemp

		IF (beta.EQ.0) THEN
			calcrate = (alpha + beta * Td(zone) ** (-1 * gamma)) *
     &					(zetaCR(zone)/1.3D-17) * CRatten(zone)
		ELSE
			nu0 = ((2.0 * gamma * kbol * bind1 )/(pi*pi*mH*beta))**0.5
			calcrate = nu0 * frac * dexp(-bind1 / Tcr) *
     &					(zetaCR(zone)/9.7D-18) * CRatten(zone)
		END IF

		IF (.NOT. CRdesorp) calcrate = 0.0

C 23) Depletion terms for ions.  From Bergin's gasgr.f code (mostly)
C	  Cfac = factor which takes into account the charge of the accreting particle.
C	  Eqn. 9 from Willacy et al. 1998 (a = 1000 A = 0.1 um = 1e-5 cm)
C      <*> Depends on grain size.

	ELSE IF (rtype.EQ.23) THEN
		vel = ((8.0 * kbol * Tg(zone))/(pi * beta * mH))**0.5
		stick = gamma

C       Back out the grain size from freezeout efficiency:
C   Changed to base on local dust surface area correction, KRS 5/20/22
C		grnfoc = 1D-5*freezeeffic**(-1.0/1.5)
    grnfoc = 1D-5 * sigadjust**0.5
C		stickingco = 0.3
		stickingco = 1.0

		Cfac = 1 + 16.71D-4 / (grnfoc * Tg(zone))
		calcrate = Cfac * alpha * vel * stick * freezeeffic * stickingco * ndust

C 30) Manually calculated photodissociation rates.  Otherwise rtype = 13
C     numbers > 30 are for multiple paths.  Path # = rtype - 29
	ELSE IF (rtype.GE.30 .AND. rtype.LT.40) THEN
		calcrate = photorate(r1, rtype)

C 41) Adding in grain surface reactions. Grain properties need to be calculated.
C     Give E_des and E_act.
C     alpha = branching ratio
C     beta = nothing
C	  gamma = E_act
C      <*> Depends on grain size.

 	ELSE IF (rtype.EQ.41) THEN
	    calcrate = 0.0
	    if (zone .gt. 3) then
			Nsites = 1E6		! This should really be an input parameter, also used for uvfield.f
            if (freezeeffic .ne. 1.0) then
C                sigmagr = pi*(rgr*freezeeffic**(-1.0/1.5))**2
								sigmagr = sigadjust * pi * rgr
                Nsites = Nsites * sigmagr/(pi*rgr**2)  ! scale number of sites by increase/decrease in surface area per grain (grain growth)
            end if

			nsfsites = 1.5e15   ! cm^-2 (TA)
			a_spacing = 1E-8    ! Barrier between adjacent sites (Ang).

			strg1=trim(r1)
			strg1=adjustl(strg1)
			strg2=trim(r2)
			strg2=adjustl(strg2)

			! Take grain off the end
			mass_i = masscalculator(strg1)
			mass_j = masscalculator(strg2)

			mu_i_j = mass_i*mass_j/(mass_i+mass_j)*mH
			kappa_ij = exp(-2*(a_spacing/(hp/(2*pi)))*
     &			 (2*mu_i_j*kbol*gamma)**0.5)

	        beind1 = -1
			beind2 = -1
			bind1 = 0.0
			bind2 = 0.0

            CALL ispecies(r1, size(bindspec), bindspec, beind1)
            CALL ispecies(r2, size(bindspec), bindspec, beind2)
            bind1 = bindenerg(beind1)
            bind2 = bindenerg(beind2)
!                  print *, beind1,beind2,bind1,bind2

			if (bind1 .EQ. 0 .or. bind2 .EQ. 0) then
				write(*,*) "Can't find B.E. in 6grainbe.inp -gr.",r1,r2
				stop
			end if

			v_zeri = 1.0/((2*nsfsites*kbol*bind1/(pi*pi*mH*mass_i))**0.5)
			v_zerj = 1.0/((2*nsfsites*kbol*bind2/(pi*pi*mH*mass_j))**0.5)

			! Assuming: E_b = 0.3*E_d, HHL
			t_hop_i  = v_zeri*exp(bind1*0.3/Td(zone))
			t_hop_j  = v_zerj*exp(bind2*0.3/Td(zone))
			R_diff_i = 1.0/(Nsites*t_hop_i)
			R_diff_j = 1.0/(Nsites*t_hop_j)

            ndustdens = grnfrac*ndens*ndust  ! Grain abundance * number H2 * epsilon
            if (freezeeffic .ne. 1.0) then
                ndustdens = ndustdens * (freezeeffic**(3.5/1.5))  ! change in grain number density per spatial volume with growth.
            end if

C			write(*,'(A30)') '-----------------------------'
C			write(*,'(A13,A13)') strg1,strg2
C			write(*,'(E13.3,E13.3)') bind1,bind2
C			write(*,'(E13.3,1x,E13.3,1x,E13.3,1x,E13.3)') kappa_ij, mu_i_j/mh, R_diff_i, R_diff_j
C			write(*,*) mass_i, mass_j
C			write(*,*) t_hop_i,t_hop_j
C			write(*,'(A30)') '-----------------------------'

			if (ndens .gt. 1E+1) then
				calcrateq = 0.0
				! Quantum tunneling for H
				select case(strg1)
                    case('H(gr)') ! Quantum tunneling for H
                        t_quant = v_zeri*exp((2*a_spacing/(hp/(2*pi)))
     &                       *(2*mass_i*mH*kbol*bind1*0.3)**0.5)
                        R_diff_i_q = 1.0/(Nsites*t_quant)

                        if (strg2 .eq. 'H(gr)') then
                            calcrateq=kappa_ij*(2.0*R_diff_i_q)
     &                           /ndustdens
                        else if (strg2 .eq. 'D(gr)') then
                            calcrateq=kappa_ij*(2.0*R_diff_i_q)
     &                          /ndustdens
                        else
                            calcrateq=kappa_ij*(R_diff_i_q)/ndustdens
                        endif

                        calcrateo = kappa_ij*(R_diff_i+R_diff_j)
     &						/ndustdens

                    case('D(gr)') ! Quantum tunneling for D
                        t_quant = v_zeri*exp((2*a_spacing/(hp/(2*pi)))
     &                      *(2*mass_i*mH*kbol*bind1*0.3)**0.5)
                        R_diff_i_q = 1.0/(Nsites*t_quant)
                        if (strg2 .eq. 'H(gr)') then
                            calcrateq=kappa_ij*(2.0*R_diff_i_q)
     &                          /ndustdens
                        else if (strg2 .eq. 'D(gr)') then
                            calcrateq=kappa_ij*(2.0*R_diff_i_q)
     &                          /ndustdens
                        else
                                calcrateq=kappa_ij*(R_diff_i_q)/ndustdens
                        endif

                        calcrateo = kappa_ij*(R_diff_i+R_diff_j)
     &						/ndustdens

					case default
						calcrateo = kappa_ij*(R_diff_i+R_diff_j)
     &							/ndustdens
					end select
					calcrate = max(calcrateo,calcrateq)
					if (calcrate .lt. 1e-30) calcrate = 0.0
			else
				calcrate = 0.0  ! Not enough gas to compute grain chemistry.
			endif
		else
		    calcrate = 0.0  ! In the upper atmosphere, no grain chem computed.
		endif

		calcrate = calcrate*alpha
                print *,'41 ',calcrate,r1,r2,alpha
		if (calcrate .lt. 1e-30) calcrate = 0.0

C 42) Grain surface reactions. Temperature dependent barrier (H + CO...)
C     Give E_des and E_act.
C     alpha = E_desi
C     beta = E_desj
C	  gamma = E_act

 	ELSE IF (rtype.EQ.42) THEN
        print *, 'Grain reactions with T-dependent barrier. Coming soon!'
        stop
		calcrate = 0.0

C 43) Ice photodissociations (not used).
C     Approximate cross section 3e-18 (ranges from 1-5e-18).
C
C     alpha = E_desi
C     beta = E_desj
C	  gamma = E_act
C
c     Soon make alpha beta gamma involve branching ratios and cross sections

 	ELSE IF (rtype.EQ.43) THEN
		calcrate = totflux(zone)*3e-18
		if (calcrate .lt. 1e-30) calcrate = 0.0

C 46) Read CR rate ionization directly in from env file.
C    Same as 1) Cosmic ray ionization (r2 = CRP):

 	ELSE IF (rtype.EQ.46) THEN
		calcrate = alpha * (zetaCR(zone)/1.3D-17)
		if (calcrate .lt. 1e-30) calcrate = 0.0
		if (.NOT. CRionization) calcrate = 0.0

C 47) Cosmic-ray desorption  from Hasegawa et al. - just evaporation
C     rate at 70 K times the time spent at 70 K when hit by cosmic
C     ray.  so Tcr = 70 K.  time spent at 70 K = frac (e = 102)
C     beta = 0 is the old e = -2 (modified rates from Bringa & Johnson 04)
C	  From Bergin's gasgr.f code

C     LIC:  Rates computed from the original Morfill et al. 1976 spectrum,
C     which has an equivalent H2 ionization rate of 9.66e-18 s^-1, using
C     the CR ionization cross sections + secondary ionization formalism of
C     Padovani et al 2009.
C
C     beta != 0.0 is the Bringa & Johnson 2004 Work
C     beta == 0.0 appears never to be used.

 	ELSE IF (rtype.EQ.47) THEN

		bind1 = 0.0
		CALL ispecies(r1, size(bindspec), bindspec, beind1)
		bind1 = bindenerg(beind1)
		if (bind1 .EQ. 0 .or. beind1 .eq. 0) then
			write(*,*) "Can't find B.E. in 6grainbe.inp.",r1
			stop
		end if

		! Increasing H(gr) and D(gr) B.E. to simulate energy between chemi and physisorption.
		! Using physisorption here
		if ((r1 == 'H(gr)        ').and.(Td(zone).GE.physlim)) bind1 = chemtemp
		if ((r1 == 'D(gr)        ').and.(Td(zone).ge.physlim)) bind1 = chemtemp

	 	IF (beta.EQ.0) THEN
			calcrate = (alpha + beta * Td(zone) ** (-1 * gamma))*(zetaCR(zone)/1.3D-17)
		ELSE
			nu0 = ((2.0 * gamma * kbol * bind1)/(pi*pi*mH*beta))**0.5
			calcrate = nu0 * frac * dexp(-bind1 / Tcr)*(zetaCR(zone)/9.7D-18)
		END IF

		IF (.NOT. CRdesorp) calcrate = 0.0

		if (calcrate .lt. 1e-30) calcrate = 0.0

C 48) Radionuclide ionization (scaling taken from cosmic ray ionization rates).
C     Turn off/on in 5flags.inp.
C     LIC 8/4/2013.
C         ->  Alpha is the difference in ionization rates between H2 and
C             species X (He is 0.84, for example -- Umebayashi et al. 2013).
C         ->  Beta is the 26 Al abundance.

 	ELSE IF (rtype.EQ.48) THEN
	    calcrate = 0.0
	    IF (incl_radionuc .eqv. .TRUE.) then
			calcrate = alpha*RNrate(zone)
		    if (calcrate .lt. 1e-30) calcrate = 0.0
		END IF

C 49) X-ray ionization of Helium (r2 = XRAY):
	ELSE IF (rtype.EQ.49) THEN
		calcrate = xrayionHe_rate()
		print *, 'He X-ray ionization rate at zone ', zone, ': ', calcrate
		if (calcrate .lt. 1e-29) calcrate = 0.0

C 50) Reactions with ortho/para dependent pathways.  H2D+ and CH2D+ are supported.
C
C     alpha = reaction rate
C     beta = reaction energy barrier
C     Hugo et al. 2009, first column is rate, second column is energy required
C     third column is in the form of (H2D+,H2)
C        1 = (o,o)
C        2 = (o,p)
C        3 = (p,p)
C        4 = (p,o)

	ELSE IF (rtype.EQ.50) THEN
C	    opratH2 = 9.0*dexp(-170.5/Tg(zone))
        calcrate = 0.0

        opratH2 = 5.3534e-3+(3.0346-5.3534e-3)/(1.0+(96.5330/Tg(zone))**3.5096)

		if (Tg(zone) .lt. 80.0) opratH2=9.0*dexp(-170.5/Tg(zone))
		if (opratH2 .lt. 1E-3) opratH2 = 1E-3  ! Based on Flower et al. 2006 Fig1

	    f_orthoH2 = opratH2/(1.0+opratH2)
		f_paraH2 = 1.0-f_orthoH2

        if (r1(1:4) .eq. 'H2D+') THEN

            opratH2Dp = -1.6977e-2 + (3.0375+1.6977e-2)/
     &			(1.0+(47.9640/Tg(zone))**3.0692)
            if (opratH2Dp .lt. 0.0) opratH2 = 0.0  ! From Lee & Bergin 2014, but cutting off below zero.
            f_orthoH2Dp = opratH2Dp/(1.0+opratH2Dp)
            f_paraH2Dp = 1.0-f_orthoH2Dp

            if (gamma.eq.1.0) then
                fro = f_orthoH2Dp*f_orthoH2
            else if (gamma.eq.2.0) then
                fro = f_orthoH2Dp*f_paraH2
            else if (gamma.eq.3.0) then
                fro = f_paraH2Dp*f_paraH2
            else if (gamma.eq.4.0) then
                fro = f_paraH2Dp*f_orthoH2
            else
            write(*,*) 'Cannot identify o/p identities for reactants.',
     &			gamma,r1
                stop
            end if

            calcrate = alpha*dexp(-beta/Tg(zone))*fro

        else if (r1(1:5) .eq. 'CH2D+') THEN

C  Based on Roeuff 2013, only considering the o/p of H2:

            if (gamma.eq.1.0) then
                fro = f_orthoH2
            else if (gamma.eq.2.0) then
                fro = f_paraH2
            else
                write(*,*) 'Cannot identify o/p identities for reactants.',
     &			gamma,r1
                stop
            end if

            calcrate = alpha*dexp(-beta/Tg(zone))*fro
        else
            print *, 'No ortho/para information provided.'
            stop
        end if

C 51) Calculate RN loss on the fly (assumes maximal escape and no neighboring contribution to the rate.)
C     Currently obsolete.
 	ELSE IF (rtype.EQ.51) THEN

	    X26AL = 1.75E-10

		E_Mg = 1.808 ! g/cm^2
		E_Bt =  0.66 ! g/cm^2
		E_Bt = 0.473 ! g/cm^2
		E_An = 0.512 ! g/cm^2

c Beta particle range taken from electron data in Padovani 2009 Figure 8
C Photon cross sections are from Finocchi and Gail 1997

		sig_Mg = 12.6
		sig_Bt = 0.14
		sig_An = 6.8

C First escape upwards:
c		sigH = RNatten(zone,1)
		if (zone .eq. 1) then
			sigH = 0.0
		else if (zone .eq. Nz) THEN
			sigH = rho(zone) * 2.0*(zAU(zone-1)-zAU(zone))*1.496e13
		else
c			colden_tot = colden_tot + rho(i) * 2 * (zcm(i) - zcm(i-1))
			sigH = rho(zone) * ((zAU(zone-1)-zAU(zone))/2.0+
     &			 (zAU(zone)-zAU(zone+1))/2.0)*1.496e13
		endif

C treating escape of beta particle and photon separately:
		Eatten1 = E_Mg*(1.0-exp(-sigH/sig_Mg)) + E_Bt*(1.0-exp(-sigH/sig_Bt)) +
     &			 E_An*(1.0-exp(-sigH/sig_An))*(1.0-exp(-sigH/sig_Bt))
		f26Al1 = 0.5*Eatten1/(E_Mg+E_Bt+E_An)

		Eatten2 = E_Mg*(1.0-exp(-sigH/sig_Mg))
		f26Al2 = Eatten2/E_Mg

		calcrate = alpha * RadNuc * beta/X26AL * (f26Al1*0.82+f26Al2*0.18)
		if (calcrate .lt. 1e-30) calcrate = 0.0

C 99) alpha = k, no other calculations needed
	ELSE IF (rtype.EQ.99) THEN
		calcrate = alpha

C Ignore Time-dependent Reaction for now
 	ELSE IF (rtype .LT. 0) THEN
		calcrate = 0.0

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
c        write(*,*) 'colden1: ',colden1(i),'colden2: ',colden2(i)
		CRatten(i) = 0.5 * (exp(-1 * colden1(i)/sigmaCR) +
     &				 exp(-1 * colden2(i)/sigmaCR))

c		print *, 'CR Atten for zone ',i,' = ', CRatten(i)
	END DO

	RETURN
	END

C..............................................................................
C
C RN Attenuation
C LIC - 4/1/2013
C
C Calculate the vertical surface density for RN loss.
C
C..............................................................................
	SUBROUTINE calcRNatten()

	INCLUDE "environ.h"
	INCLUDE "constants.h"

	INTEGER i
	DOUBLE PRECISION colden_tot, colden1(Nz), colden2(Nz)

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
c        write(*,*) 'colden1: ',colden1(i),'colden2: ',colden2(i)
		RNatten(i,1) = colden1(i)
		RNatten(i,2) = colden2(i)

c		print *, 'CR Atten for zone ',i,' = ', CRatten(i)
	END DO

	RETURN
	END


C..............................................................................
C
C X-Ray Grain Heating
C
C returns photons/second/grain absorbed by grain lattice
C
C
C Assumptions:
C kTxr = 2 keV
C E1 = 0.8 keV, E2 = 20 keV
C Nsec (number of secondary ionizations per unit photoelectron energy) = 30
C E range 1 - 30 keV
C..............................................................................
	DOUBLE PRECISION FUNCTION xrayheat_rate()
	IMPLICIT NONE

	DOUBLE PRECISION xrayinteg
	DOUBLE PRECISION xrayheatinteg
	DOUBLE PRECISION kTmin, kTmax
	double precision energyvec(100),fluxvec(100)
	double precision a
	integer nptsE
	integer i,b

	PARAMETER (kTmin = 1.0, kTmax = 20.0)
c	write(*,*) 'present xrayheat_rate, begin'

c	CALL qromb(xrayheatinteg, kTmin, kTmax, xrayheat_rate)


	nptsE = 100
	a = (kTmax-kTmin)/(nptsE-1)

	energyvec = (/(((i-1)*a+kTmin),i=1,nptsE)/)   ! assign to y in increments of 1.5 starting at 1.5

	do b=1,nptsE

		fluxvec(b) = xrayheatinteg(energyvec(b))

	enddo

c	write(*,*) energyvec, fluxvec
	call avint(fluxvec, energyvec, nptsE, kTmin, kTmax, xrayheat_rate)
c	write(*,*) 'xry heat: ',xrayheat_rate
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
	DOUBLE PRECISION FUNCTION xrayheatinteg(E)

C Common Blocks
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	DOUBLE PRECISION E

c	DOUBLE PRECISION Lv0, Lv, F0, F
c	DOUBLE PRECISION KTXR, E1, E2, Nsec, c
c	PARAMETER (KTXR = 2, E1 = 0.8, E2 = 10.0)
c	PARAMETER (c = 2.99792458D10)

C External Functions
	DOUBLE PRECISION sigma_hot, xray_flux


c	Lv0 = Lxray/KTXR * (dexp(-E1/KTXR) - dexp(-E2/KTXR))
c	Lv = Lv0 * dexp(-E/KTXR)
c	F0 = Lv/(8*pi*Rs(zone)**2*AU**2 * E * keV)
c	xrayinteg = Nsec * Habssigma(E) * F0 * dexp(-sigmaE(E)*Nrz(zone))

c	xrayinteg = Nsec * Habssigma(E) * nxray_photons(zone) * c
c    write(*,*) 'xrayheatinteg',E
	xrayheatinteg = sigma_hot(E) * xray_flux(E)

	if (xrayheatinteg .LE. 1.0e-30) then
		xrayheatinteg = 1.0e-30
	endif

c		write(*,*) 'present 4',E,xrayheatinteg
c	print *, E, xrayinteg, sigmaE(E), Nrz(zone)
c	write(*,*) E, sigma_hot(E), xray_flux(E)
c   Want to integrate sigma_hot(E) * xray_flux(E) between

	RETURN
	END

C..............................................................................
C
C X-Ray heating absorption cross x
C
C xcross: xray cross section from Glassgold, Najita, and Igea in cm^-2.
C given E=ephot (keV);
C
C
C flux_int: contains all the pieces for the integral of the flux over
C X-ray energy.  This does not include the attenuation or any constants.
C However, it does include the absorption probability and the exponential
C term.  See notes.
C ***/
c
C..............................................................................
	DOUBLE PRECISION FUNCTION sigma_hot(ephot)
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	double PRECISION cross,ephot
	double PRECISION lac, agr_lac, agr_lac2
	double PRECISION rho_gr, mugr, pabs
c	REAL*8 e_sil, xs_sil
c	DIMENSION e_sil(2000), xs_sil(2000)
	integer nln, n
c	write(*,*) 'present 6'
C   Read in Cross-X data.  May want to hard-code these arrays in for
C   efficiency sake later.
c	open (unit=59,file='xs_format.data',status='old')
cc	write(*,*) 'present 7'
c	fend = 0
c	nln = 1
c	do while (nln .lt. 2001)
c		read(59,130,iostat=fend) e_sil(nln), xs_sil(nln)
c		nln = nln + 1
cc		write(*,*) e_sil(nln), xs_sil(nln)
c	enddo
c	close(59)
c  130    format(1x,e13.6,2x,e13.6)
	n = 2000

	call linint(ephot,cross,e_sil/1000.0,xs_sil,n)
c	write(*,*) e_sil/1000.0, xs_sil
c	write(*,*) 'E, cross ',ephot,cross
c	cross = 2.22e-22 * (ephot)**(-2.485)
	rho_gr = 3.3
	mugr = 172.0*mh
c	/** 1000 A grain **/ ?? not 50 A?
	agr = 50.0 * 1.0e-8

C   calculating probability that an x-ray photon will interact with a 50 A.
C   subunit!
C   lac is the linear attenuation coeffiecient -- see Dwek and Smith eq 1,2
c		write(*,*) 'present 9'
	lac = rho_gr/mugr * cross
c	write(*,*) 'cross: ',cross
	agr_lac = agr * lac
	agr_lac2 = 2.0 * agr_lac

	if (agr_lac .lt. 0.01) then
		pabs = agr_lac2*agr_lac + exp(-agr_lac2)*(agr_lac2+1.0) - 1.0
		pabs = pabs/(agr_lac*agr_lac2)
	else
		pabs = 1.333*agr_lac
	endif

	sigma_hot = pi * agr**2 * pabs

	RETURN
	END

C..............................................................................
c     Linear Interpolation
c     Given a value of x return a value of y based on interpolation
c     within a table of y values (ytab) corresponding to the x values
c     contained in the array xtab.  The subroutine assumes that the
c     values in xtab increase monotonically
c
c    John Mahaffy 2/12/95
C..............................................................................
      subroutine linint(x,y,xtab,ytab,n)
	  double precision x,y
	  integer n,i,i1
      double precision xtab(n),ytab(n)

      if (x.lt.xtab(1).or.x.gt.xtab(n)) then
         write(6,*) 'x = ', x, '  is out of interpolation range'
         stop
      endif

      do 100 i=2,n
         if (x.le.xtab(i)) go to 200
  100    continue

  200 i1=i-1
      wx=(x-xtab(i1))/(xtab(i1+1)-xtab(i1))
      y=(1-wx)*ytab(i1)+wx*ytab(i1+1)
      return
      end

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
	double precision energyvec(100),fluxvec(100)
	double precision a
	integer nptsE
	integer i,b

c	CALL qromb(xrayinteg, kTmin, kTmax, xrayion_rate)

	nptsE = 100
	a = (kTmax-kTmin)/(nptsE-1)

	energyvec = (/(((i-1)*a+kTmin),i=1,nptsE)/)   ! assign to y in increments of 1.5 starting at 1.5

	do b=1,nptsE
		fluxvec(b) = xrayinteg(energyvec(b))
	enddo

	call avint(fluxvec, energyvec, nptsE, kTmin, kTmax, xrayion_rate)
c        xrayion_rate = xrayion_rate * 10.0
c        write(*,*) 'ENHANCED X-rays! 1/25/14'
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
	DOUBLE PRECISION FUNCTION xrayionHe_rate()
	IMPLICIT NONE

	DOUBLE PRECISION xrayintegHe
	EXTERNAL xrayintegHe
	DOUBLE PRECISION kTmin, kTmax
	PARAMETER (kTmin = 1.0, kTmax = 20.0)
	double precision energyvec(100),fluxvec(100)
	double precision a
	integer nptsE
	integer i,b

c	CALL qromb(xrayinteg, kTmin, kTmax, xrayion_rate)

	nptsE = 100
	a = (kTmax-kTmin)/(nptsE-1)

	energyvec = (/(((i-1)*a+kTmin),i=1,nptsE)/)   ! assign to y in increments of 1.5 starting at 1.5

	do b=1,nptsE
		fluxvec(b) = xrayintegHe(energyvec(b))
	enddo

	call avint(fluxvec, energyvec, nptsE, kTmin, kTmax, xrayionHe_rate)
c        xrayionHe_rate = xrayionHe_rate * 10.0

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

c	DOUBLE PRECISION Lv0, Lv, F0, F
	DOUBLE PRECISION KTXR, E1, E2, Nsec, c
c	PARAMETER (KTXR = 2, E1 = 0.8, E2 = 10.0, Nsec = 30)
c	PARAMETER (c = 2.99792458D10)
	PARAMETER (Nsec = 30)
C External Functions
	DOUBLE PRECISION sigmaE, Habssigma, xray_flux


c	Lv0 = Lxray/KTXR * (dexp(-E1/KTXR) - dexp(-E2/KTXR))
c	Lv = Lv0 * dexp(-E/KTXR)
c	F0 = Lv/(8*pi*Rs(zone)**2*AU**2 * E * keV)
c	xrayinteg = Nsec * Habssigma(E) * F0 * dexp(-sigmaE(E)*Nrz(zone))

c	xrayinteg = Nsec * Habssigma(E) * nxray_photons(zone) * c
c	print *, 'a',xray_flux(E)
	xrayinteg = Nsec * Habssigma(E) * xray_flux(E)
c	print *, 'b', xray_flux(E)
c	print *, Habssigma(E)
	if (xrayinteg .LE. 1.0e-30) then
		xrayinteg = 1.0e-30
	endif

c	print *, E, Nsec, xrayinteg

	RETURN
	END


C..............................................................................
C Integrand for X-ray ionization rate calculation for Helium
C
C Assumptions:
C kTxr = 2 keV
C E1 = 0.8 keV, E2 = 10 keV
C Nsec (number of secondary ionizations per unit photoelectron energy) = 30
C E range 1 - 10 keV
C..............................................................................
	DOUBLE PRECISION FUNCTION xrayintegHe(E)

C Common Blocks
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	DOUBLE PRECISION E

c	DOUBLE PRECISION Lv0, Lv, F0, F
	DOUBLE PRECISION KTXR, E1, E2, Nsec, c
c	PARAMETER (KTXR = 2, E1 = 0.8, E2 = 10.0, Nsec = 30)
c	PARAMETER (c = 2.99792458D10)
	PARAMETER (Nsec = 30)
C External Functions
	DOUBLE PRECISION sigmaE, Hesigma, xray_flux

	xrayintegHe = Nsec * Hesigma(E) * xray_flux(E)
c	print *, 'He: ', xray_flux(E), Hesigma(E)
c	print *, Hesigma(E)
	if (xrayintegHe .LE. 1.0e-30) then
		xrayintegHe = 1.0e-30
	endif

c	print *, E, Nsec, xrayinteg

	RETURN
	END
C..............................................................................
C
C He Sigma
C
C Calculates the He photoionization cross section using the fitting formula of Yan, Sadeghpour and Dalgarno 1998.
C
C..............................................................................
	DOUBLE PRECISION FUNCTION Hesigma(Einp)

	DOUBLE PRECISION Einp
	DOUBLE PRECISION seqterm, ahe(6),xe
	INTEGER i

	xe = Einp*1e3/24.58

	ahe(1) = -4.7416
	ahe(2) =  14.8200
	ahe(3) = -30.8678
	ahe(4) =  37.3584
	ahe(5) = -23.4585
	ahe(6) =  5.9133

	seqterm = 0.0
	do i =1,6
	   seqterm = seqterm + ahe(i)/xe**(real(i)/2.0)
	enddo

	Hesigma = 733.0/(Einp**(7.0/2.0))*(1.0+seqterm)*1D-24
c	print *,'Helium: ',Einp,Hesigma

	RETURN
	END


C..............................................................................
C
C This function reads in the uv field from the photons file.  Mostly taken from
C Ted's readuv.f file.
C
C..............................................................................
C
C Input parameter(s):
C
C inpfile			== photons.dat filename
C
C..............................................................................
	subroutine readxray(inpfile)

C Include the environmental variables and xsect variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

	CHARACTER*80 inpfile

	double precision findrad, radius, UVtemp
	double precision radius_low, radius_high, radius_match
	character*900 buff
	character*20 words(90)
	integer nwords, read_flg, j, i, inw, nj, sflg
	integer header, nlen, fend

c External Functions
	double precision atod

C Format(s)
 20		format(A900)

c Use the first radius
	radius = Rs(zone)

	maxzone = 0
	read_flg = 0
	sflg = 0
	header = 0
	i = 1

	open(unit = 99, file = trim(inpfile))
c	print *, 'READXRAY: opening photon file:'
	read(99,20,end=100) buff
	write(*,'(a70)') buff

C JF 12/8/10 - Read in until found radius is larger than given radius.  Then
C pick the radius that is closest to the given radius
	fend = 0
	radius_low = 0
	radius_high = 0
	do while (fend .EQ. 0)
		read(99, '(A900)', iostat=fend) buff
		if (index(buff, 'Radius(AU)') .ne. 0) then
			call gwordslong(buff, nwords, words)
			findrad = atod(words(2))
			if (findrad .GE. radius) then
				radius_high = findrad
				exit
			else
				radius_low = findrad
			endif
		else
			cycle
		endif
	enddo
	close(99)

C Find radius closest to given radius
	if ((radius - radius_low) .LT. (radius_high - radius)) then
		radius_match = radius_low
	else
		radius_match = radius_high
	endif

c Match last radius if within 1%
	if ( (radius_high .EQ. 0) .AND.
     +		((radius - radius_low)/radius .LE. 1e-2) ) then
		radius_match = radius_low
c Match first radius if within 1%
	else if ( (radius_low .EQ. 0) .AND.
     +		((radius_high - radius)/radius .LE. 1e-2) ) then
		radius_match = radius_high
	endif

C Now read in the UV field for the found radius
	open(unit = 99, file = trim(inpfile))
	read(99,20,end=100) buff

	fend = 0
	do while (fend .EQ. 0)
	   read(99, '(A900)', iostat=fend) buff
	   call gwordslong(buff, nwords, words)

c  loop exit construct
		if ((words(1) .eq. 'Radius(AU)').and.(sflg .eq. 1)) then
c		    print *, 'readxray: done reading in x-ray photons'
		    xnlam = i
		    exit
		end if

		if (words(1) .eq. 'Radius(AU)') then
			read_flg = read_flg + 1
			findrad = atod(words(2))
			if (findrad .EQ. radius_match) then
c				print *, 'readxray: isolated radius in file'
				sflg = 1
			end if
		end if

		if (sflg .eq. 1) header = header + 1

		if ((sflg .eq. 1).and.(header .gt. 3)) then
			xraylevels(i) = atod(words(1))
			nj = nwords
			Itemp = 0.0

		    do inw = 2, nj
				j = inw - 1
				xrayfield(j,i) = atod(words(inw)) ! WACKY!
c				print *, j,i, xrayfield(j,i)

c set the last zone uv field to be equal to the zone above it to avoid
c extrapolation errors.
c JF, 12/8/10 - Corrected so that it only is used if the last zone has a higher
c flux than the one immediately above it.

c				if (inw .eq. nj .and. xrayfield(j-1,i) .LT. xrayfield(j,i)) then
c					if (ndust .eq. 1) then
c						xrayfield(j,i) = xrayfield(j-1,i)
c					else
c						xrayfield(j,i) = 0.0
c					endif
c				endif

c Set erroneous values to 0
				if (xrayfield(j,i) .lt. 0.0) xrayfield(j,i) = 0.0

		    end do
		    i = i + 1
		end if
	  end do
100	continue
	if (sflg .eq. 0) then
		print *, 'ERROR in reading x-ray photons.dat file'
		print *, 'unable to find radius in x-ray photons.dat file.  Exiting'
		stop
	end if
c	print *, 'MULTIPLYING X-RAYS by 2.0 -- NOTE!'
	xnlam = xnlam - 1
	nj = nj - 1
c	write(*,*) 'xnalm',xnlam,nlam
	close(99)
c	stop

	return
	end


C..............................................................................
C
C This function reads in the uv field from the photons file.  Mostly taken from
C Ted's readuv.f file.
C
C..............................................................................
C
C Input parameter(s):
C
C inpfile			== photons.dat filename
C
C..............................................................................
	subroutine readsimonxray(inpfile)

C Include the environmental variables and xsect variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

	CHARACTER*80 inpfile
	double precision radius
	integer j
	integer header, nlen, fend, nln
	double precision arrr,arrz,zeta_cr_simon

c External Functions
	double precision atod

C Format(s)
 20		format(A900)

c Use the first radius
	radius = Rs(zone)
c	inpfile = 'hd100546_xrgrid.txt'
	j = 1
c	print *, inpfile
	open (unit=67,file=trim(inpfile),status='old')
	do while (fend .EQ. 0)
		read(67,130, iostat=fend) arrr, arrz, zeta_cr_simon
		if (abs(arrr - radius)/radius .le. 5e-3) then
			xrayratesimon(j) = zeta_cr_simon
			j = j + 1
c			print *, 'saved: '
		endif
		if ((abs(arrr - radius)/radius .ge. 5e-3).and.(arrr .gt. radius)) exit
c		print *, abs(arrr - radius)/radius, arrr, arrz, zeta_cr_simon
	enddo
	close(67)

 130    format(e13.7,1x,e13.7,1x,e13.7)

	print *, 'READXRAY: done reading xrayfile from Simon B.'
c	print *, xrayratesimon
	return
	end


C..............................................................................
C
C This function reads in the radionuclide calculations.
C
C..............................................................................
C
C Input parameter(s):
C
C inpfile			== photons.dat filename
C
C..............................................................................
	subroutine readRN(inpfile)

C Include the environmental variables and xsect variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

	CHARACTER*80 inpfile

	double precision findrad, radius
	double precision radius_low, radius_high, radius_match
	character*900 buff
	character*20 words(90)
	integer nwords, read_flg, j, i, inw, nj, sflg
	integer header, nlen, fend
	integer qct
	double precision z_tmp(150), zeta_tmp(150)
	double precision zhere, zetav

c External Functions
	double precision atod

C Format(s)
 20		format(A900)

c Use the first radius
	radius = Rs(zone)

	maxzone = 0
	read_flg = 0
	sflg = 0
	header = 0
	i = 1

	open(unit = 99, file = trim(inpfile))
c	print *, 'READRN: opening radionuclide file:'
	read(99,20,end=100) buff
c	write(*,'(a70)') buff

C JF 12/8/10 - Read in until found radius is larger than given radius.  Then
C pick the radius that is closest to the given radius
	fend = 0
	radius_low = 0
	radius_high = 0
	do while (fend .EQ. 0)
		read(99, '(A900)', iostat=fend) buff
		if (index(buff, 'Radius(AU)') .ne. 0) then
			call gwordslong(buff, nwords, words)
c			print *,'words:',words,words(2)
			findrad = atod(words(2))
			if (findrad .GE. radius) then
				radius_high = findrad
				exit
			else
				radius_low = findrad
			endif
		else
			cycle
		endif
	enddo
	close(99)

C Find radius closest to given radius
	if ((radius - radius_low) .LT. (radius_high - radius)) then
		radius_match = radius_low
	else
		radius_match = radius_high
	endif

c Match last radius if within 1%
	if ( (radius_high .EQ. 0) .AND.
     +		((radius - radius_low)/radius .LE. 1e-2) ) then
		radius_match = radius_low
c Match first radius if within 1%
	else if ( (radius_low .EQ. 0) .AND.
     +		((radius_high - radius)/radius .LE. 1e-2) ) then
		radius_match = radius_high
	endif

C Now read in the UV field for the found radius
	open(unit = 99, file = trim(inpfile))
	read(99,20,end=100) buff

	fend = 0
	do while (fend .EQ. 0)
	   read(99, '(A900)', iostat=fend) buff
	   call gwordslong(buff, nwords, words)

c  loop exit construct
		if ((words(1) .eq. 'Radius(AU)').and.(sflg .eq. 1)) then
c		    print *, 'ReadRN: done reading in radionuclide file.'
		    exit
		end if

		if (words(1) .eq. 'Radius(AU)') then
			read_flg = read_flg + 1
			findrad = atod(words(2))
			if (findrad .EQ. radius_match) then
c				print *, 'readrn: isolated radius in file',radius_match
				sflg = 1
			end if
		end if

		if (sflg .eq. 1) header = header + 1

  		if ((sflg .eq. 1).and.(header .gt. 1)) then

			do qct=2,nwords
				z_tmp(qct-1) = atod(words(qct))
			enddo

			read(99, '(A900)', iostat=fend) buff
	        call gwordslong(buff, nwords, words)

			do qct=2,nwords
				zeta_tmp(qct-1) = atod(words(qct))
			enddo

			do j=1,Nz
			    zhere = zAU(j)
				do qct=1,nwords-1
				  if ((zhere .ge. z_tmp(qct)*0.99).and.(zhere .LT. z_tmp(qct+1))) then
					zetav=(z_tmp(qct+1)-zhere)/(z_tmp(qct+1)-z_tmp(qct))*zeta_tmp(qct)
     &					+ (zhere-z_tmp(qct))/(z_tmp(qct+1)-z_tmp(qct))*zeta_tmp(qct+1)
                    RNrate(j)=zetav
				  elseif (zhere.gt.z_tmp(nwords-1)) then
				    RNrate(j)=zeta_tmp(nwords-1)
				  endif
				enddo
                    !        RNrate(j) = RNrate(j)*10.0
			enddo
		endif
	  end do
100	continue
	if (sflg .eq. 0) then
		print *, 'ERROR in reading radionuclide file'
		print *, 'unable to find radius in radionuclide file.  Exiting'
		stop
	end if

	close(99)

	return
	end


C..............................................................................
C
C Return X-ray flux at given E value
C
C..............................................................................
	DOUBLE PRECISION FUNCTION xray_flux(Einp)

	INCLUDE "environ.h"
	INCLUDE "rates.h"

	DOUBLE PRECISION Einp
	INTEGER i, j

C Find the right energy level
	DO i=1,xnlam
		IF (xraylevels(i) .EQ. Einp) THEN

			xray_flux = xrayfield(zone, i)
c			write(*,*) 'xraylevels(i) ',Einp,xray_flux
			RETURN
		ELSE IF (Einp .LT. xraylevels(i)) THEN
			j = i-1
c			write(*,*) 'OR: xraylevels(i)',xraylevels(i),Einp,j
			EXIT
		ENDIF
	ENDDO
c write(*,*) '1**', Einp,xray_flux
c   Interpolate between Energies to get flux:
	xray_flux = ((Einp-xraylevels(j))*xrayfield(zone,j+1) +
     &			 (xraylevels(j+1)-Einp)*xrayfield(zone,j)) /
     &			 (xraylevels(j+1)-xraylevels(j))
c	write(*,*) 'xraylevels(i) ', Einp,xray_flux
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
	DOUBLE PRECISION x

C	Evals(1) = 1.0
C	Evals(2) = 1.5
C	Evals(3) = 2.0
C	Evals(4) = 3.0
C	Evals(5) = 4.0
C	Evals(6) = 5.0
C	Evals(7) = 6.0
C	Evals(8) = 8.0
C	Evals(9) = 10.0
C	Evals(10) = 15.0
C	Evals(11) = 20.0
C
C	sigmavals(1) = 1.14D+01
C	sigmavals(2) = 2.93D+00
C	sigmavals(3) = 1.11D+00
C	sigmavals(4) = 2.81D-01
C	sigmavals(5) = 1.05D-01
C	sigmavals(6) = 4.91D-02
C	sigmavals(7) = 2.63D-02
C	sigmavals(8) = 9.82D-03
C	sigmavals(9) = 4.56D-03
C	sigmavals(10) = 1.13D-03
C	sigmavals(11) = 4.18D-04
C
C
C	DO i=1,nvals
C		IF (Evals(i) .EQ. Einp) THEN
C			Habssigma = 1D-24 * sigmavals(i)
C			RETURN
C		ELSE IF (Einp < Evals(i)) THEN
C			j = i-1
C			EXIT
C		ENDIF
C	ENDDO
C
C	Habssigma = ((Einp-Evals(j))*sigmavals(j+1) +
C     &			 (Evals(j+1)-Einp)*sigmavals(j)) /
C     &			 (Evals(j+1)-Evals(j))
C
C	Habssigma = Habssigma * 1d-24

C Now using cross section from Yan, Sadeghpour and Dalgarno et al. 1998
	x = Einp*1000.0/15.4
	Habssigma = (45.57*(1 - 2.003*x**-0.5 - 4.806/x + 50.577*x**-1.5 -
     &				171.044*x**-2 + 231.608*x**-2.5 -
     &				81.885*x**-3)/Einp**3.5)*1e-24 ! Barns -> cm^2

c	print *,'Hydrogen: ',Einp,Habssigma

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
      PARAMETER (EPS=1.e-3, JMAX=25, JMAXP=JMAX+1, K=5, KM=K-1)
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
	CHARACTER*10 molstring, dstr, dn,f, fis
	DOUBLE PRECISION atomarr(8), masses(8)
	INTEGER ct, x,xold, oldatomindex, numberofatoms, oldnumberofatoms
	DOUBLE PRECISION C_n, O_n, N_n, He_n, H_n, Mg_n, Na_n, Ca_n, S_n, D_n, Y_n, Z_n
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
	D_n = 2.0
	Y_n = 13.0
	Z_n = 18.0
C changed Z to same mass as O KRS 2/5/15
	masses = [H_n, C_n, O_n, N_n, S_n, D_n, Y_n, Z_n]

	wasdigit = .false.

c	          [H,C,O,N,S,D,Y,Z]
	atomarr = [0,0,0,0,0,0,0,0]
	ct = 1
	dstr = 'h' ! Initialize while loop.

	do while ((dstr .ne. '').and.(dstr .ne. '('))
	     dstr = molstring(ct:ct)
         select case(dstr)
         	case('0':'9')
			if (wasdigit) then
			    xold = x
			    read( dstr, '(i10)' )  x
                formt='('//trim('i1,i1')//')'
                write(str,formt) xold,x
                str=adjustl(str)
				read(str, '(i20)') numberofatoms
			    atomarr(oldatomindex) = atomarr(oldatomindex)-
     &				 xold+numberofatoms
		        oldnumberofatoms = numberofatoms
				wasdigit = .true.
			else
			    read( dstr, '(i10)' )  x
			    numberofatoms = x
			    atomarr(oldatomindex) = atomarr(oldatomindex)
     &				 + numberofatoms-1
		        oldnumberofatoms = numberofatoms
				wasdigit = .true.
			endif

         	case default
		    wasdigit = .false.
	          select case (dstr)       ! number of type integer
	          case ('H')            ! all values below 0
				   dn = molstring(ct+1:ct+1)
				   if (dn .eq. 'e') then
						!write(*,*) 'He!'
						ct = ct + 1
			       else
	               atomarr(1) = atomarr(1) + 1
				   oldatomindex = 1
				   endif
	          case ('C')
	               atomarr(2) = atomarr(2) + 1
				   oldatomindex = 2
	          case ('O')
	               atomarr(3) = atomarr(3) + 1
				   oldatomindex = 3
			  case ('N')
	               atomarr(4) = atomarr(4) + 1
				   oldatomindex = 4
			  case ('S')
				   dn = molstring(ct+1:ct+1)
				   if (dn .eq. 'i') then
						!write(*,*) 'Si!'
						ct = ct + 1
			       else
	               atomarr(5) = atomarr(5) + 1
				   oldatomindex = 5
				   endif
	          case ('D')
				   dn = molstring(ct+1:ct+1)
	               atomarr(6) = atomarr(6) + 1
				   oldatomindex = 6
	          case ('Y')
	               atomarr(7) = atomarr(7) + 1
				   oldatomindex = 7
	          case ('Z')
	               atomarr(8) = atomarr(8) + 1
				   oldatomindex = 8
			  !CASE ('+')
			  ! write(*,*) d
			  !CASE ('E')
		      CASE DEFAULT

	          END SELECT
		 end select
		 if (sum(atomarr) .eq. 0.0) then
			write(*,*) 'Could not compute mass of ice species: ',molstring,'stopping.'
			stop
		 endif

c		 write(*,*) d
		 ct = ct + 1

	enddo

	masscalculator = sum(atomarr*masses)

	return
	end

	subroutine avint(ftab, xtab, ntab, a, b, result)
      !
      !***********************************************************************
      !
      !! AVINT estimates the integral of unevenly spaced data.
      !
      !
      !  Discussion:
      !
      !    The method uses overlapping parabolas and smoothing.
      !
      !  Reference:
      !
      !    Philip Davis and Philip Rabinowitz,
      !    Methods of Numerical Integration,
      !    Blaisdell Publishing, 1967.
      !
      !    P E Hennion,
      !    Algorithm 77,
      !    Interpolation, Differentiation and Integration,
      !    Communications of the Association for Computing Machinery,
      !    Volume 5, page 96, 1962.
      !
      !  Modified:
      !
      !    30 October 2000
      !
      !  Parameters:
      !
      !    Input, real FTAB(NTAB), the function values,
      !    FTAB(I) = F(XTAB(I)).
      !
      !    Input, real XTAB(NTAB), the abscissas at which the
      !    function values are given.  The XTAB's must be distinct
      !    and in ascending order.
      !
      !    Input, integer NTAB, the number of entries in FTAB and
      !    XTAB.  NTAB must be at least 3.
      !
      !    Input, real A, the lower limit of integration.  A should
      !    be, but need not be, near one endpoint of the interval
      !    (X(1), X(NTAB)).
      !
      !    Input, real B, the upper limit of integration.  B should
      !    be, but need not be, near one endpoint of the interval
      !    (X(1), X(NTAB)).
      !
      !    Output, real RESULT, the approximate value of the integral.
      !
        implicit none
      !
        integer ntab
      !
        double precision a
        double precision atemp
        double precision b
        double precision btemp
        double precision ca
        double precision cb
        double precision cc
        double precision ctemp
        double precision ftab(ntab)
        integer i
        integer ihi
        integer ilo
        integer ind
        double precision result
        double precision sum1
        double precision syl
        double precision term1
        double precision term2
        double precision term3
        double precision x1
        double precision x2
        double precision x3
        double precision xtab(ntab)

        if ( ntab < 3 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'AVINT - Fatal error!'
          write ( *, '(a,i6)' ) '  NTAB is less than 3.  NTAB = ', ntab
          stop
        end if

        do i = 2, ntab

          if ( xtab(i) <= xtab(i-1) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'AVINT - Fatal error!'
            write ( *, '(a)' ) '  XTAB(I) is not greater than XTAB(I-1).'
            write ( *, '(a,i6)' ) '  Here, I = ', I
            write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
            write ( *, '(a,g14.6)' ) '  XTAB(I) =   ', xtab(i)
            stop
          end if

        end do

        result = 0.0E+00

        if ( a == b ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'AVINT - Warning!'
          write ( *, '(a)' ) '  A = B, integral=0.'
          return
        end if
      !
      !  If A > B, temporarily switch A and B, and store sign.
      !
        if ( a > b ) then
          syl = b
          b = a
          a = syl
          ind = -1
        else
          syl = a
          ind = 1
        end if
      !
      !  Bracket A and B between XTAB(ILO) and XTAB(IHI).
      !
        ilo = 1
        ihi = ntab

        do i = 1, ntab
          if ( xtab(i) >= a ) then
            exit
          end if
          ilo = ilo + 1
        end do

        ilo = max ( 2, ilo )
        ilo = min ( ilo, ntab-1 )

        do i = 1, ntab
          if ( b >= xtab(i) ) then
            exit
          end if
          ihi = ihi - 1
        end do

        ihi = min ( ihi, ntab-1 )
        ihi = max ( ilo, ihi-1 )
      !
      !  Carry out approximate integration from XTAB(ILO) to XTAB(IHI).
      !
        sum1 = 0.0E+00

        do i = ilo, ihi

          x1 = xtab(i-1)
          x2 = xtab(i)
          x3 = xtab(i+1)

          term1 = ftab(i-1) / ((x1-x2)*(x1-x3))
          term2 = ftab(i) / ((x2-x1)*(x2-x3))
          term3 = ftab(i+1) / ((x3-x1)*(x3-x2))

          atemp = term1 + term2 + term3
          btemp = -(x2+x3)*term1-(x1+x3)*term2-(x1+x2)*term3
          ctemp = x2*x3*term1+x1*x3*term2+x1*x2*term3

          if ( i <= ilo ) then
            ca = atemp
            cb = btemp
            cc = ctemp
          else
            ca = 0.5E+00 * ( atemp + ca )
            cb = 0.5E+00 * ( btemp + cb )
            cc = 0.5E+00 * ( ctemp + cc )
          end if

          sum1 = sum1
     &	        + ca * ( x2**3 - syl**3 ) / 3.0E+00
     &	        + cb * 0.5E+00 * ( x2**2 - syl**2 )
     &	        + cc * ( x2 - syl )

          ca = atemp
          cb = btemp
          cc = ctemp

          syl = x2

        end do

        result = sum1
     &	        + ca * ( b**3 - syl**3 ) / 3.0E+00
     &	        + cb * 0.5E+00 * ( b**2 - syl**2 )
     &	        + cc * ( b - syl )
      !
      !  Restore original values of A and B, reverse sign of integral
      !  because of earlier switch.
      !
        if ( ind /= 1 ) then
          ind = 1
          syl = b
          b = a
          a = syl
          result = -result
        end if

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
