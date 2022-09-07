c--------------------------------------------------------------------

C..............................................................................
C
C dtaudtemp
C lic 1/22
C calculate the cooldown rate due to x-ray thermal desorption
C..............................................................................

	SUBROUTINE dtaudtemp(itime,zonedefT,tmpder1,nomo_tmpder1,Numr1)
	IMPLICIT NONE

C Common Blocks
	INCLUDE "Fcn.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	INTEGER n_mols
	PARAMETER (n_mols = 4)
c	DOUBLE PRECISION nh2(n_mols), theta(n_mols)
	INTEGER itime, imol
	CHARACTER*13 species
	CHARACTER*13 specarr(n_mols)
c	CHARACTER*13 dataarr(n_mols,3)
	REAL dataarr(n_mols,3),dataarr_rc(n_mols,3)
	INTEGER r1index, grindex,p
	DOUBLE PRECISION nr1, ngr, Numr1(n_mols)
	DOUBLE PRECISION bind_factor,zonedefT
	PARAMETER (bind_factor = 1.47)
	DOUBLE PRECISION polint,numtot
	double precision tmpder1(67,3),nomo_tmpder1(67,3)
C	Calculate Total Number of X per GRain:
C		CO
C		NH3
C		O
C		C
C		H
c	write(*,*) 'really, why did i ever get here?'
	specarr(1) ='CO(gr)       '
	specarr(2) ='NH3(gr)      '
	specarr(3) ='O(gr)        '
	specarr(4) ='C(gr)        '
c	specarr(5) ='H(gr)        '
c	print *,itime,zonedefT
C Molecular data:
	dataarr(1,1:2) = [28.0,1210.0]
	dataarr(2,1:2) = [17.0,1110.0]
	dataarr(3,1:2) = [16.0,800.0]
	dataarr(4,1:2) = [12.0,800.0]
c	dataarr(5,1:2) = [1.0,350.0]

	numtot = 0.0
C Calculate grain species:
	do imol = 1, n_mols
		species = specarr(imol)
		CALL ispecies(species, ns, s, r1index)
		CALL ispecies('GRAIN        ', ns, s, grindex)
		IF (itime .le. 2) THEN
			nr1 = abundances(r1index, timestep)
			ngr = abundances(grindex, timestep)
		ELSE
			nr1 = abundances(r1index, timestep-1)
			ngr = abundances(grindex, timestep-1)
		END IF
		Numr1(imol) = nr1/ngr
		numtot = numtot + Numr1(imol)
		dataarr(imol,3) = Numr1(imol)
c		write(*,*) species,itime,Numr1(imol)
	enddo
	if (itime.eq. 1) write(*,*) 'Total number surf spec: ',numtot

c	old_Numr1
c	old_tmpder1
c	nomo_Numr1
c	nomo_tmpder1

	if (itime.eq.1) then
c		Initialize Array:
	    do p = 1,67
			nomo_tmpder1(p,1:3) = [0.,0.,0.]
		enddo

		dataarr_rc = dataarr

c		Set N_spec = Zero for all Species
		do p = 1,n_mols
			dataarr_rc(imol,3) = 0.0
		enddo

		call diffusion(dataarr,bind_factor,specarr,zonedefT,tmpder1)
		nomo_tmpder1 = tmpder1
	endif
c	write(*,*) 'before diff'
	call diffusion(dataarr,bind_factor,specarr,zonedefT,tmpder1)
c	write(*,*) 'after diff'
c	IF (itime .EQ. 40) stop
	RETURN
	END


C..............................................................................
C
C JN's xray heat diffusion code modified to calculate dt/dT
C lic 1/23
C
C..............................................................................

	SUBROUTINE diffusion(dataarr,bind_factor,specnamearr,
     +      	zonedefT,tmpder1)
	IMPLICIT NONE

c     	calculate cooling of spot heated grain via thermal diffusion
c     	using an implicit method.  include evaporative cooling by CO.

c     	notes:  to avoid solution of a non-linear equation (the heat
c     	capacity is temperature-dependent), i am letting the heat
c     	capacity lag in the time domain.  see if this is stable!
c	6/30/00: it works w/o evap cooling!
c	7/7/00:  it works w/ evap and rad cooling!
c	1/21/12: limit number of molecules that are available to evaporate

	logical done, writeout
	real pi,kappa,kB,kA,JpeV,radiusgran,dx,rcs,fkapitz,fspec
	real rlameff,Qe,sigSB,Asite, v0,Tbind,dTemp,Tmax,Aeff,grvol
	real Edvol,Edep,Rad,Tspot,time,dt,Tdmin,Eevap, Erad,dEevap
	real dErad,dEheat,rCv,Hcap,therm_cond,sp_heat,alpha,beta
	real game,game_mol,gamr,Eheat,Tdiff,Tdx,Tdxi,Ttarg,Tdiffi
	double precision zonedefT, sigadjust
	integer ngr, nhot, nmol,num_mol,i,imol,iwcount,ct, imoln
	real mirror, dt0, EdeV

	sigadjust = 2.5e-13
	C corresponds to old radiusgran=50e-10
	if (incl_locdust) then
		sigadjust = locdust(zone)
	end if


	parameter (mirror=0)	 ! endpoint: 1 for insulated, 0 for sink.

	parameter (dt0=5.e-11)   ! initial time step size (s)
c	parameter (dt0=5.e-11)   ! initial time step size (s)

	parameter (pi=3.14159)
	parameter (EdeV=300.)    ! Energy deposited in hot grain (eV)
	parameter (JpeV=1.6e-19) ! Joules per eV
c	parameter ()     ! equilibrium grain temperature (K)
c	parameter (radius=1.5e-8) ! grain radius (m)
	parameter (radiusgran=sqrt((1e-14*dustfrac)) ! grain radius (m)
c dustfrac = <r^2>/(0.1 um)^2
	parameter (dx=2.*radiusgran) ! grid size = 2*grain radius (m)
	parameter (rcs=0.02)	 ! ratio r_contact/radius   Can only be 0.02 or 0.1?
	parameter (fkapitz=1.0)  ! reduction factor for conductivity
	parameter (fspec=1.0)  ! reduction factor for specific heat
	parameter (rlameff=1.0)  ! eff. lambda for radiative cooling (lam/100um)
	parameter (Qe=0.01*radiusgran/1.e-7/rlameff) ! radiative efficiency
	parameter (sigSB=5.67e-8)! SB constant (J/m2/s/K^4)

	parameter (Asite=7.e-20) ! area per binding site for CO (m^2)
	parameter (v0=2.e12)     ! nu0 for CO (s^-1)
	parameter (Tbind=960.)	 ! binding energy for CO (K)
	parameter (kB=1.38e-23)  ! Boltzmann constant (W/K)

	parameter (ngr=101, nhot=51)
	parameter (nmol=5)

	real Tmp(ngr),Tmpn(ngr),Acf(ngr),Bcf(ngr),Ccf(ngr),Dcf(ngr)
	real v0_in(nmol), Tbind_in(nmol),tempHotOld
	real amw(nmol), bind_factor
        real Revap(nmol),Nevap(nmol),Evap(nmol),Nevapcent(nmol)
	character*13 name(nmol)
	REAL dataarr(nmol,3),dtmpev

c  LIC+++
	real Tbind_in_tmp(nmol),Nevap_store(nmol)
	real Revap_tmp(nmol),Nevap_tmp(nmol),Evap_tmp(nmol)
	character*13 name_tmp(nmol)
	real maxspec(nmol),maxspec_tmp(nmol)
	CHARACTER*13 specarr(nmol),specnamearr(nmol)
	integer num_mol_orig,tempercounter
	double precision timetemp(300000,2)
	real derivTemp
	integer tc,numinter,z
	real Tend
	double precision timeevalu,tmpevalu,timestep
	double precision,dimension(:,:), allocatable :: derstore
c	double precision,dimension(:,:), allocatable :: tmpder1
	double precision tmpder1(67,3)
c  ++++++
	parameter (dTemp=10.,Tmax=150.)


	Tend=zonedefT
	do i = 1, num_mol
		Tbind_in(i) = dataarr(i,2)
		amw(i) = dataarr(i,1)
		maxspec(i) = dataarr(i,3)
		name(i) = specnamearr(i)

		Tbind_in(i) = Tbind_in(i) * bind_factor
		v0_in(i) = 2.0 * 1.5e15 * 1.38e-16 * Tbind_in(i)
        v0_in(i) = v0_in(i) / (pi*pi*amw(i)*1.6e-24)
		v0_in(i) = v0_in(i)**0.5

c       for CO need to go to measured value of vibration freq.
        if (name(i) .eq. 'CO(gr)       ') then
            v0_in(i) = 2.0e12 * (Tbind_in(i)/960.0)**0.5
		end if
        print *, name(i), 'nu0= ', v0_in(i), ' Eb = ', Tbind_in(i)
	end do
		num_mol_orig = num_mol

c	Using Joel's correction for effective area
	if (rcs.eq.0.02) then
		Aeff = 0.0587*radiusgran*radiusgran
	elseif (rcs.eq.0.1) then
		Aeff = 0.265*radiusgran*radiusgran
	else
		print *, 'rcs out of bounds'
		stop
	endif

c	rcontact=rcs*radiusgran
c	Aeff=pi*rcontact*rcontact     ! eff. cross sec area for heat transport (m^2)
	Aeff=Aeff*fkapitz

	grvol = 4.*pi/3.*(radiusgran**3)  ! grain volume
	Edep = EdeV*JpeV	      ! energy deposited in hot grain (J)
	Edvol = Edep/grvol

c this many binding sites w/ this vibration freq
c Revap * Eperevap, evergy removed by evaporations as a function of time
	do imol = 1, num_mol
		Revap(imol) = 4.*pi*(radiusgran**2)/Asite*v0_in(imol)
        Evap(imol) = 4.*pi*(radiusgran**2)/
     +      Asite*v0_in(imol)*Tbind_in(imol)*kB
	end do
	Rad = 4.*pi*(radiusgran**2)*Qe*sigSB

c energy lost to radiation?

	do imoln = 1, num_mol
		write(*,*) Revap(imoln), Evap(imoln), Rad
	end do


c     	Initialize grain temperature array
	do i=1,ngr
		Tmp(i)=Tend
	enddo
c	Tmp(nhot)=Tspot
	Tmp(nhot)=Tspot(Edvol,Tend,fspec)
c       write(*,*) Tend, Tmp(nhot)

c     	Start heat diffusion loop:
	time=0.
	dt=dt0
	Tdmin=0.0001 	!(for heat sink; 50A, 300eV)
c	Tdmin=0.000005 	!(for heat sink; 50A, 300eV)
	do imol= 1, num_mol
           Nevap(imol)=0.
		   Nevapcent(imol) = 0.
	end do
	Eevap=0.
	Erad=0.
	done=.false.
	iwcount=0

	tempercounter = 1
	tempHotOld = Tmp(nhot)
c	write(*,*) 'prehra'
	do while (.not.done)

	    dEevap=0.
	    dErad=0.
	    dEheat=0.

c      	    Setup heat diffusion coefficient arrays
c        write(*,*) 'wee1a'
	    do i=1,ngr

			rCv = sp_heat(Tmp(i),fspec)
			Hcap=grvol*rCv
			kappa = therm_cond(Tmp(i))
			kA=kappa*Aeff

			alpha = dt/dx*kA/Hcap
	        beta  = dt/Hcap
			game = 0.0
c			;write(*,*) 'wee1b',i
			do imol = 1, num_mol
				game_mol  = Evap(imol)*exp(-Tbind_in(imol)/Tmp(i))
				game  = game + game_mol
			enddo

	        gamr  = Rad*Tmp(i)**4

			if (i.eq.1) then
	    	    Acf(1) = 0.			! not used in Tridag
	    	    Bcf(1) = (1.+2.*alpha)
				if (mirror.eq.0) then 	! fixed end temperature
	    	        Ccf(1) = -alpha
	    	        Dcf(1) = Tmp(1) + alpha*Tend - beta*(game+gamr)
				else 			! zero derivative (insulated)
					Ccf(1) = -2.*alpha
					Dcf(1) = Tmp(1) - beta*(game+gamr)
				endif
			elseif (i.eq.ngr) then
	    	    Bcf(ngr) = (1.+2.*alpha)
	    	    Ccf(ngr) = 0.		! not used in Tridag
				if (mirror.eq.0) then 	! fixed end temperature
	    	        Acf(ngr) = -alpha
	    	        Dcf(ngr) = Tmp(ngr) + alpha*Tend - beta*(game+gamr)
				else 			! zero derivative (insulated)
	    	        Acf(ngr) = -2.*alpha
	    	        Dcf(ngr) = Tmp(ngr) - beta*(game+gamr)
				endif
			else
				Acf(i) = -alpha
				Bcf(i) = (1.+2.*alpha)
				Ccf(i) = -alpha
				Dcf(i) = Tmp(i) - beta*(game+gamr)
			endif
c           write(*,*) 'wee1c',i
	        do imol = 1, num_mol
	          Nevap(imol)=Nevap(imol)+
     +              Revap(imol)*exp(-Tbind_in(imol)/Tmp(i))*dt
c				print *, 'Nevap = ', Nevap(imol)
			end do

	        dEevap= dEevap+ game
	        dErad = dErad + gamr
	        Eevap = Eevap + game*dt
	        Erad  = Erad  + gamr*dt
    		dEheat= dEheat + Eheat(Tmp(i),Tend,fspec)*grvol
c			if (i.eq.nhot) write(*,*) 'Evap: ', game*dt,'Erad: ', gamr*dt
	    enddo

c	    Now solve the set of linear equations
c	    cwrite(*,*) 'rahh!'
	    call tridag(Acf,Bcf,Ccf,Dcf,Tmpn,ngr)

	    time=time+dt

c	    Test if done?

	    Tdiff=0.

	    do i=1,ngr
		     Tdiffi = abs(Tmp(i)-Tmpn(i))
		     Tdiff = Tdiff + Tdiffi**2
		     Tmp(i)=Tmpn(i)
	    enddo

	    Tdiff = sqrt(Tdiff)

	    Tdx=0.

	    do i=1,ngr-1
		     Tdxi = abs(Tmp(i+1)-Tmpn(i))
		     Tdx = Tdx + Tdxi**2
	    enddo


		derivTemp = dt/(Tmp(nhot)-tempHotOld)
		if (derivTemp .gt. 1e30) derivTemp=1e30
		if (derivTemp .lt. 1e-30) derivTemp=1e-30

	    Tdx = sqrt(Tdx)
        timetemp(tempercounter, 1:2) = [time, Tmp(nhot)]
c		write (*,*) timetemp(tempercounter, 1:2)

		if (Tmp(nhot).lt.120.) dt=dt0*10.
		if (Tmp(nhot).lt.90.) dt=dt0*100.
		if (Tmp(nhot).lt.80.) dt=dt0*1000.
	    if (Tmp(nhot).lt.60.) dt=dt0*10000.
	    if (Tmp(nhot).lt.40.) dt=dt0*100000.
	    if (Tmp(nhot).lt.20.) dt=dt0*1000000.
	    if (Tdiff.lt.Tdmin) done=.true.
c	    if (Tdx.lt.0.001) done=.true.

c		If within 2% of the end, you are done!
        if (Tmp(nhot) .le. Tend*1.02) done=.true.

c	    Ttarg=Tmax-(iwcount*dTemp)

c		write(*,*) Tmp(nhot), Tend, Tmax, dt,Tmp(nhot)-tempHotOld

c		if (Tmp(nhot).le.Ttarg) then
c		writeout=.true.



c		iwcount=iwcount+1
c	    else
c		writeout=.false.
c	    endif
c        write(*,*) 'wee1'
c	    topsecret=14.588
c	    if (Tmp(nhot).le.topsecret) then
c		done=.true.
c	    endif


c 800	format (f9.4,1x)
        ct = 1
		do imol = 1, num_mol
			if (Nevap(imol).lt.maxspec(imol)) then
c			Nevap(imol)=Nevap(imol)+
c     +			Revap(imol)*exp(-Tbind_in(imol)/Tmp(i))*dt
c			print *, 'Still evaporate = ',Nevap(imol),maxspec(imol),name(imol)
				Nevap_tmp(ct) = Nevap(imol)
				Revap_tmp(ct) = Revap(imol)
				Tbind_in_tmp(ct) = Tbind_in(imol)
				name_tmp(ct) = name(imol)
				Evap_tmp(ct) = Evap(imol)
				maxspec_tmp(ct) = maxspec(imol)
				ct = ct + 1
			else
				Nevap_store(imol) = Nevap(imol)
			endif
		end do
		num_mol = ct - 1
		Nevap = Nevap_tmp
		Revap = Revap_tmp
		Tbind_in = Tbind_in_tmp
		name = name_tmp
		Evap = Evap_tmp
		maxspec = maxspec_tmp

		tempercounter = tempercounter + 1

		tempHotOld = Tmp(nhot)
c		write(*,*) 'wee2'
c		print *, tempercounter, Tmp(nhot)
        if (tempercounter .eq. 300000) then
			print *, 'X-ray heat diffusion did not converge, stopping.'
			stop

		endif

	enddo
	print *, 'cooldown took X steps: X=',tempercounter-1
c	diffusion = timetemp
c	write(*,*) 'posthra'
	tc = tempercounter - 1
	numinter = 69
	allocate(derstore(numinter,3))
c	allocate(tmpder1(numinter,2))

	timestep = (log10(maxval(timetemp(1:tc,1))*0.99)
     +      -log10(timetemp(1,1)*1.01))/(numinter-1)*0.95
	write(*,*) 'posthra'
c	write(*,*) timestep


	do z = 1,numinter
	    timeevalu = log10(timetemp(1,1)*1.05) + timestep*(z-1)

		call nearestneigh(log10(timetemp(1:tc,1)),timetemp(1:tc,2),tc,
     +      timeevalu,tmpevalu)
		derstore(z,1:2) = [10**(timeevalu), tmpevalu]
c		write(*,*) 10**(timeevalu), tmpevalu
	enddo

	tmpder1(:,2) = (derstore(5:numinter,1)-
     +      derstore(4:numinter-1,1)) /
     +      (derstore(5:numinter,2)-
     +      derstore(4:numinter-1,2))

	tmpder1(:,1) = (derstore(5:numinter,2)+
     +      derstore(4:numinter-1,2))/2.0

	tmpder1(:,3) = (derstore(5:numinter,1)+
     +      derstore(4:numinter-1,1))/2.0

c	write(*,*) '-----------------------'
c	write(*,*) tmpder1(1:64,1)
c	write(*,*) '-----------------------'
c	write(*,*) tmpder1(1:64,2)
c	write(*,*) '-----------------------'
	return

c	do imol = 1, num_mol_orig
c	Nevap_store(imol) = Nevap(imol)
c	enddo
c	print *,Nevap_store

	end


c--------------------------------------------------------------------
	real function therm_cond(temp)

	parameter (kord=6)	! order of fit to kappa
 	real ak(0:kord),logT
 	data ak / -3.72, -0.192, 8.03, -15.31, 11.81, -4.04, 0.512 /
c 	data ak / -3.73, 2.53, -3.34, 1.97, -0.369 /

c	logT = alog10(temp)
c	fx = ak(6)*logT**6 + ak(5)*logT**5 + ak(4)*logT**4 +
c     *		ak(3)*logT**3 + ak(2)*logT**2 + ak(1)*logT + ak(0)
c	therm_cond = 10.**fx 	! in W/cm/K
c	therm_cond = therm_cond*100.	! in W/m/K

	therm_cond = 0.3

	return
	end

c--------------------------------------------------------------------

	real function sp_heat(temp,fspec)

	parameter (Tbhc=50.)	 ! T bdry for change in Leger heat capacity law
	parameter (Cvc1=140.)    ! heat capacity = Cvc1 T^p1 (J/K/m^3); T<50K
	parameter (p1=2.)
	parameter (Cvc2=2200.)   ! heat capacity = Cvc2 T^p2 (J/K/m^3); T=50-150K
	parameter (p2=1.3)

	if (temp.le.Tbhc) then
	    sp_heat=Cvc1*fspec*temp**p1
	else
	    sp_heat=Cvc2*fspec*temp**p2
	endif

	return
	end


c--------------------------------------------------------------------

	real function Tspot(Edvol,Tend,fspec)

c	NOTE: This function should be consistent with the specific heat routine.

	parameter (bndry=5.787e6)  ! boundary for heat cap choice (J/m3)
	parameter (Tbhc=50.)	   ! T bdry for change in Leger heat capacity law
	parameter (Cvc1=140.)      ! heat capacity = Cvc1 T^p1 (J/K/m^3); T<50K
	parameter (p1=2.)
	parameter (Cvc2=2200.)     ! heat capacity = Cvc2 T^p2 (J/K/m^3); T=50-150K
	parameter (p2=1.3)

	p1p1 = p1+1
	p2p1 = p2+1
	if (Edvol.gt.bndry) then
	    aa = p2p1/(Cvc2*fspec)
	    bb = (Cvc1*fspec)/p1p1
	    cc = Edvol - bb*(Tbhc**p1p1-Tend**p1p1)
	    Tspot=(aa*cc + Tbhc**p2p1)**(1./p2p1)
	else
	    bb = p1p1/(Cvc1*fspec)
	    Tspot = (bb*Edvol + Tend**p1p1)**(1./p1p1)
	endif

	return
	end

c--------------------------------------------------------------------

	real function Eheat(Tempi,Tend,fspec)

c	NOTE: This function should be consistent with the specific heat routine.
c	calculates energy/vol stored in grain as heat

	parameter (Tbhc=50.)	   ! T bdry for change in Leger heat capacity law
	parameter (Cvc1=140.)      ! heat capacity = Cvc1 T^p1 (J/K/m^3); T<50K
	parameter (p1=2.)
	parameter (Cvc2=2200.)     ! heat capacity = Cvc2 T^p2 (J/K/m^3); T=50-150K
	parameter (p2=1.3)

	p1p1 = p1+1
	p2p1 = p2+1
	if (Tempi.gt.Tbhc) then
	    e1 = (Cvc1*fspec)/p1p1*(Tbhc**p1p1 - Tend**p1p1)
	    e2 = (Cvc2*fspec)/p2p1*(Tempi**p2p1 - Tbhc**p2p1)
	    Eheat = e1+e2
	else
	    e1 = (Cvc1*fspec)/p1p1*(Tempi**p1p1 - Tend**p1p1)
	    Eheat = e1
	endif

	return
	end

c--------------------------------------------------------------------
	logical function output(iwcount, Temp)

	parameter (dTemp=10.,Tmax=150.)

	Ttarg=Tmax-(iwcount*dTemp)
	if (Temp.le.Ttarg) then
		output=.true.
		iwcount=iwcount+1

	else
		output=.false.
	endif

	return
	end

c--------------------------------------------------------------------
      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      REAL a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
      INTEGER j
      REAL bet,gam(NMAX)
      if(b(1).eq.0.)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END

c--------------------------------------------------------------------
	  SUBROUTINE nearestneigh(xa,ya,n,x,yeval)
      INTEGER n,NMAX
      double precision x,y,yeval,xa(n),ya(n)
c     Just find the nearest value in the array as a crude interpolation
c     ponzi scheme.
      INTEGER i,m,ns
      double precision den,dif,dift
      ns=1
	  ns2=1
	  yeval = 0
	  y1 = 0
	  y2 = 0
      dif1=abs(x-xa(1))
	  dif2=dif1
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif1) then
		  ns1 = i
		  dif1 = dift
		elseif ((dift.lt.dif2).and.(dift.gt.dif1)) then
		  ns2 = i
		  dif2 = dift
        endif
11    continue
      y1 = ya(ns1)
	  y2 = ya(ns2)
	  yeval = dif1/(dif1+dif2)*y2+dif2/(dif1+dif2)*y1
	  if (x .le. minval(xa)) yeval = ya(1)
	  if (x .ge. maxval(xa)) yeval = ya(n)

      return
      END
