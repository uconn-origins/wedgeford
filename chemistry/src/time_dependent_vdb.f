C..............................................................................
C
C This function checks to see if the on grain abundance is larger than
C the grain abundance and, if so, reduces the rate of the reaction
C
C..............................................................................
	DOUBLE PRECISION FUNCTION grain_abun_check(ralpha, rbeta, species)
	implicit none
	
C Initialization of common blocks:	
	INCLUDE "Fcn.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	DOUBLE PRECISION ralpha, rbeta
	DOUBLE PRECISION nr1, ngr
	DOUBLE PRECISION calcrate
	CHARACTER*13 species
	INTEGER r1index, grindex

	calcrate = ralpha*(Tk(zone)/300.0D0)**rbeta

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
c		write(*,*) species
c		write(*,'(1pE10.3, 1X, 1pE10.3)') nr1, ngr
c		write(*, '(1pE10.3, 1X, 1pE10.3)') calcrate, calcrate * ngr / nr1
		calcrate = calcrate * ngr / nr1
	ENDIF
	
	grain_abun_check = calcrate
	
	RETURN
	END

C..............................................................................
C
C This function calculates the H2 photodissociation rate and takes into account
C the H2 self-shielding.  
C Based heavily on Ted Bergin's self_shieldH2.f code.
C
C "Taken from Draine and Bertoldi (1996), ApJ, 468, 269.
C  Also in Hollenbach and Tielens (1999), Rev. Mod. Phys., 71, 174"
C
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
	DOUBLE PRECISION H2abun, f, maxlam, etau

	DATA colden/0.0/, colden_h2/0.0/, colden_h2_tot/ntime * 0.0/


c	avshock = AvIS(zone)
c	sigmad_1000 = 2.0e-21
c	colden = avshock*2.0e21
c	taud_1000 = sigmad_1000 * colden

	k_Av0 = 4.0e-11 * g0_100

	IF (zone .EQ. 1) zcm(zone-1) = 0.0
	H2abun = rho(zone) / mH
	colden_h2 = H2abun * (zcm(zone) - zcm(zone-1))
	colden_h2_tot(timestep) = colden_h2_tot(timestep) + colden_h2

	IF (colden_h2_tot(timestep) .LE. 1.0e14) THEN
		f = 1.0
	ELSE
		f = (colden_h2_tot(timestep)/1.0e14)**(-0.75)
	END IF

c  changed to I/I0 at 1000 Ang
	maxlam = 1000.0
	call find_etau(maxlam, etau, zone)
	self_shield_H2 = k_Av0 * f * etau

c write values to selfshield.test
	IF (shieldtest) THEN
		OPEN(unit=02, file=fnameSS, status='old',
     &			ACCESS='APPEND')
		write(2,110), Rs(zone), etau, self_shield_H2, 
     &			colden_h2_tot(timestep), colden, f, avshock, timestep
		CLOSE(02)
	END IF
 110	format(1x, 'H2 ', 1(1x,e10.3), 2x, f8.3, 5(1x,e8.3), i4)

	RETURN
	END


C..............................................................................
C
C This function calculates the CO photodissociation rate and takes into account
C the CO self-shielding.  Isotopes of C (Y = 13C) and O (Z = 18O) are included.
C Based heavily on Ted Bergin's self_shield_v2.f code.
C
C "will compute the self-shielding coefficients for the CO molecule 
C using the rates given in van Dishoeck and Black 1988"
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

	double precision colden_co, colden_13co, colden_c18o
	double precision colden_h2_12, colden_h2_13, colden_h2_18
	double precision colden_co_tot(ntime), colden_13co_tot(ntime)
	double precision colden_c18o_tot(ntime)
	double precision colden_h2_tot_12(ntime), colden_h2_tot_13(ntime)
	double precision colden_h2_tot_18(ntime)
	DATA colden_co_tot/ntime * 0.0/, colden_13co_tot/ntime * 0.0/
	DATA colden_c18o_tot/ntime * 0.0/ 
	DATA colden_h2_tot_12/ntime * 0.0/, colden_h2_tot_13/ntime * 0.0/
	DATA colden_h2_tot_18/ntime * 0.0/ 

	double precision thetatot
	double precision nCO, H2abun
	double precision k_Av0, etau, maxlam
	double precision lnh2, lnco, t, u
	double precision y1,y2,y3,y4

	integer i, j, i_h2, i_co
	integer flag
	integer COindex

	IF (shieldtest) THEN
		OPEN(unit=02, file=fnameSS, status='old', ACCESS ='APPEND')
		IF (timestep.EQ.1) THEN
			WRITE(02,*)
			WRITE(02,*) 'Zone: ', zone
		END IF
	END IF

	do i = 1,ih2
	  do j = 1, ico
		theta(i,j) = 0.0
		dtheta(i,j) = 0.0
	  end do
	end do

	if (iso .eq. 12)  call vdb_coshield(theta, xnco, xnh2)
	if (iso .eq. 13)  call vdb_13coshield(theta, xnco, xnh2)
	if (iso .eq. 18)  call vdb_c18oshield(theta, xnco, xnh2)
	
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
	do i = 1, ico
		theta_highh2(i) = theta(ih2,i)
	end do
	call spline(xnco,theta_highh2,ico,1.0e30,1.0e30,dtheta_highh2)
	do i = 1, ih2
		theta_highco(i) = theta(i,ico)
	end do
	call spline(xnh2,theta_highco,ih2,1.0e30,1.0e30,dtheta_highco)

C Get the current CO abundance
	CALL ispecies('CO           ', ns, s, COindex)
	IF (timestep .EQ. 1) THEN
		nCO = abundances(COindex, timestep)
	ELSE
		nCO = abundances(COindex, timestep-1)
	END IF

c/***
c   compute the sheilding due to the dust continuum
c***/
	if (iso .eq. 12) k_Av0 = 2.039e-10 * g0_CO
	if (iso .eq. 13) k_Av0 = 2.034e-10 * g0_CO
	if (iso .eq. 18) k_Av0 = 2.035e-10 * g0_CO

c/***
c   here we have to save the column density of H2 and CO at each time step
c   and add it to the column density in the previous zone.
c   this requires the time steps to be the same in each zone!!!! which
c   should be the case - but put in a check just in case.  the check is 
c   located in the main program.
c***/
	if (zone.eq.1) zcm(zone-1) = 0.0
	if (iso .eq. 12) then
		colden_co = nCO * (zcm(zone) - zcm(zone-1))
		if (nCO .LE. 0.0) colden_co = 1.0e4 
		colden_co_tot(timestep) = colden_co_tot(timestep) + colden_co
		lnco = dlog10(colden_co_tot(timestep))
		H2abun = rho(zone) / mH
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
		H2abun = rho(zone) / mH
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
		H2abun = rho(zone) / mH
		colden_h2_18 = H2abun * (zcm(zone) - zcm(zone-1))
		colden_h2_tot_18(timestep) = colden_h2_tot_18(timestep) 
     &		+ colden_h2_18
		lnh2 = dlog10(colden_h2_tot_18(timestep))
	end if

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
		call splint(xnh2,theta_highco,dtheta_highco,
     +          ih2, lnh2, thetatot)
		self_shield_CO = thetatot * k_Av0
		flag = 1
	else if (lnh2 .ge. 23.0) then
		call splint(xnco,theta_highh2,dtheta_highh2,
     +          ico, lnco, thetatot)
		self_shield_CO = thetatot * k_Av0
		flag = 1
	end if


	if (flag .eq. 0) then
		do i = 1,ih2
			if (lnh2 .GT. xnh2(i)) i_h2 = i
		end do
		do i = 1,ico
			if (lnco .GT. xnco(i)) i_co = i
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

		self_shield_CO = k_Av0 * thetatot
	end if
	flag = 0

c ETAU
	maxlam = 960.0
	call find_etau(maxlam, etau, zone)
	
	self_shield_CO = self_shield_CO * etau

	IF (shieldtest) THEN
		if (iso .eq. 12) then
			write(2, 112) iso,zcm(zone),AvIS(zone), self_shield_CO, 
     &			colden_co_tot(timestep), colden_h2_tot_12(timestep),
     &			thetatot, etau, nCO, timestep
		end if
		if (iso .eq. 13) then
			write(2, 112) iso,zcm(zone), AvIS(zone), self_shield_CO, 
     &			colden_13co_tot(timestep), colden_h2_tot_13(timestep),
     &			thetatot, etau, nCO, timestep
		end if
		if (iso .eq. 18) then
			write(2, 112) iso,zcm(zone), AvIS(zone), self_shield_CO, 
     &			colden_c18o_tot(timestep), colden_h2_tot_18(timestep),
     &			thetatot, etau, nCO, timestep
		end if

		CLOSE(02)
	END IF

 112     format(1x,f3.0,2x,8(1x,e8.3),2x, i4)

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
c	if (x .EQ. 1500) then
c		print *, 'IN ETAU: ', x, jzone, photuv, photuv_max
c	endif
	etau = photuv/photuv_max

	return
	end

C..............................................................................
C VDB_COSHIELD - returns xnco, xnh2 and theta values for CO
C
C From van Dishoeck?
C..............................................................................
	subroutine vdb_coshield(theta, xnco, xnh2)
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
	theta(1,2) = 9.68e-1
	theta(1,3) = 7.764e-1
	theta(1,4) = 3.63e-1
	theta(1,5) = 7.013e-1
	theta(1,6) = 1.295e-2 
	theta(1,7) = 1.738e-3
	theta(1,8) = 9.985e-5 

	theta(2,1) = 8.215e-1  
	theta(2,2) = 7.916e-1  
	theta(2,3) = 6.160e-1 
	theta(2,4) = 2.749e-1  
	theta(2,5) = 5.351e-2  
	theta(2,6) = 1.065e-2 
	theta(2,7) = 1.519e-3  
	theta(2,8) = 8.818e-5 

	theta(3,1) = 7.160e-1  
	theta(3,2) = 6.900e-1  
	theta(3,3) = 5.360e-1 
	theta(3,4) = 2.359e-1  
	theta(3,5) = 4.416e-2  
	theta(3,6) = 8.769e-3 
	theta(3,7) = 1.254e-3  
	theta(3,8) = 7.558e-5 

	theta(4,1) = 3.500e-1  
	theta(4,2) = 3.415e-1  
	theta(4,3) = 2.863e-1 
	theta(4,4) = 1.360e-1  
	theta(4,5) = 2.500e-2  
	theta(4,6) = 4.983e-3 
	theta(4,7) = 7.151e-4  
	theta(4,8) = 3.796e-5 

	theta(5,1) = 4.973e-2  
	theta(5,2) = 4.877e-2  
	theta(5,3) = 4.296e-2 
	theta(5,4) = 2.110e-2  
	theta(5,5) = 4.958e-3  
	theta(5,6) = 9.245e-4 
	theta(5,7) = 1.745e-4  
	theta(5,8) = 8.377e-6 

	theta(6,1) = 1.310e-4  
	theta(6,2) = 1.293e-4  
	theta(6,3) = 1.160e-4 
	theta(6,4) = 6.346e-5  
	theta(6,5) = 1.822e-5  
	theta(6,6) = 6.842e-6 
	theta(6,7) = 3.622e-6  
	theta(6,8) = 3.572e-7 

	return
	end

C..............................................................................
C VDB_13COSHIELD - returns xnco, xnh2 and theta values for 13CO
C
C From van Dishoeck?
C..............................................................................
	subroutine vdb_13coshield(theta, xnco, xnh2)
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
	theta(1,2) = 0.9887 
	theta(1,3) = 0.9159  
	theta(1,4) = 0.6845   
	theta(1,5) = 0.2610  
	theta(1,6) = 0.06032  
	theta(1,7) = 7.788e-3
	theta(1,8) = 3.402e-4 

	theta(2,1) = 0.8181    
	theta(2,2) = 0.8083    
	theta(2,3) = 0.7463   
	theta(2,4) = 0.5324    
	theta(2,5) = 0.2185    
	theta(2,6) = 0.04961  
	theta(2,7) = 6.431e-3  
	theta(2,8) = 2.859e-4 

	theta(3,1) = 0.7011    
	theta(3,2) = 0.6922    
	theta(3,3) = 0.6386   
	theta(3,4) = 0.4540    
	theta(3,5) = 0.1835    
	theta(3,6) = 0.04160  
	theta(3,7) = 5.556e-3  
	theta(3,8) = 2.404e-4

	theta(4,1) = 0.3599    
	theta(4,2) = 0.3573    
	theta(4,3) = 0.3392   
	theta(4,4) = 0.2585    
	theta(4,5) = 0.1202    
	theta(4,6) = 0.02767  
	theta(4,7) = 3.389e-3  
	theta(4,8) = 1.346e-4 

	theta(5,1) = 6.037e-2  
	theta(5,2) = 5.993e-2  
	theta(5,3) = 5.929e-2 
	theta(5,4) = 5.423e-2  
	theta(5,5) = 3.320e-2  
	theta(5,6) = 6.691e-3
	theta(5,7) = 7.129e-4  
	theta(5,8) = 1.858e-5 

	theta(6,1) = 8.019e-4  
	theta(6,2) = 8.014e-4  
	theta(6,3) = 7.979e-4 
	theta(6,4) = 7.640e-4  
	theta(6,5) = 5.197e-4  
	theta(6,6) = 1.115e-4 
	theta(6,7) = 1.500e-5  
	theta(6,8) = 6.254e-7 

	return
	end
	
C..............................................................................
C VDB_C18OSHIELD - returns xnco, xnh2 and theta values for C18O
C
C From van Dishoeck?
C..............................................................................
	subroutine vdb_c18oshield(theta, xnco, xnh2)
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
	theta(1,2) = 0.9897 
	theta(1,3) = 0.9243  
	theta(1,4) = 0.6673   
	theta(1,5) = 0.2921  
	theta(1,6) = 0.09464  
	theta(1,7) = 0.001451
	theta(1,8) = 7.45e-4  

	theta(2,1) = 0.8088    
	theta(2,2) = 0.8       
	theta(2,3) = 0.7450   
	theta(2,4) = 0.5405    
	theta(2,5) = 0.2383    
	theta(2,6) = 0.07686  
	theta(2,7) = 0.01194   
	theta(2,8) = 6.310e-4 

	theta(3,1) = 0.7032    
	theta(3,2) = 0.6953    
	theta(3,3) = 0.6477    
	theta(3,4) = 0.4708    
	theta(3,5) = 0.2091    
	theta(3,6) = 0.06811  
	theta(3,7) = 0.01042   
	theta(3,8) = 5.071e-4

	theta(4,1) = 0.3611    
	theta(4,2) = 0.3587    
	theta(4,3) = 0.3424   
	theta(4,4) = 0.2655    
	theta(4,5) = 0.1371    
	theta(4,6) = 0.04805  
	theta(4,7) = 0.006614  
	theta(4,8) = 2.436e-4 

	theta(5,1) = 6.093e-2  
	theta(5,2) = 6.059e-2  
	theta(5,3) = 6.005e-2 
	theta(5,4) = 5.592e-2  
	theta(5,5) = 4.069e-2  
	theta(5,6) = 1.480e-2
	theta(5,7) = 1.640e-3  
	theta(5,8) = 3.276e-5 

	theta(6,1) = 9.061e-4  
	theta(6,2) = 9.061e-4  
	theta(6,3) = 9.042e-4 
	theta(6,4) = 8.855e-4  
	theta(6,5) = 7.410e-4  
	theta(6,6) = 2.968e-4 
	theta(6,7) = 3.616e-5  
	theta(6,8) = 8.619e-7 

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

	DOUBLE PRECISION a_gr, Mlayer, ngr
	INTEGER i, grindex

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

	a_gr = 1.0e-5
	Mlayer = 4 * pi * a_gr**2 / ((3e-8)**2)
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

