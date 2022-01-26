
C------------------------------------------------------------------------------
C
C "Driver" subroutine to the stiff SODE solver "DVODEPK".
C
C------------------------------------------------------------------------------
C Input parameter(s):
C 
C neq      == number of equations,
C
C tfirst   == initial time step,
C
C nt       == amount of time steps to be taken,
C
C yy(neq)  == initial abundances at t=0
C
C tlast    == last time moment,
C
C rtol1    == relative tolerance parameter,
C
C atol1    == absolute tolerance parameter,
C
C------------------------------------------------------------------------------
C
C Output parameter(s):
C
C abundances(ns,nt) == species abundances
C
C------------------------------------------------------------------------------
C Used subroutine(s) and function(s) (alphabetically): 
C 
C Function(s): calcrate, Fcn, Jacobian
C Subroutine(s): dvodpk
C 
C------------------------------------------------------------------------------
      subroutine run_chemistry(neq,tfirst,nt,yy,tlast,rtol1,atol1)
      implicit none

C Initialization of common blocks:
      INCLUDE "Fcn.h"
      INCLUDE "dvode.h"
	  INCLUDE "environ.h"
	  INCLUDE "BindingE.h"

C Global variable(s):
      integer neq, nt
      real*8 yy, tfirst, tlast, rtol1(nspec), atol1(nspec)
      dimension yy(nspec)

C Local variable(s):
	  DOUBLE PRECISION calcrate

      INTEGER n_mols,dummyy
	  PARAMETER (n_mols = 4)
	  real grnfrac
      integer mf, itol, itask, istate, iopt, lrw, liw, IPAR, i, k, 
     1   neqz, timedep, j
      real*8 RPAR, t, tout, tstep, atol(nspec), rtol(nspec)
      dimension RPAR(nspec), IPAR(nspec*nss)
	  double precision dtaudt(67,3),nomo_dtaudt(67,3),Ngrsp(n_mols)
	  double precision n_dtaudt(67,3),old_Ngrsp(n_mols),old_dtaudt(67,3)
	  double precision diffabun
	  INTEGER imol
	  CHARACTER*13 species
	  CHARACTER*13 specarr(n_mols)
	  DOUBLE PRECISION numtot, nr1, ngr, Numr1(n_mols)
	  INTEGER r1index, grindex
	  
C External routines:
      external Fcn, Jacobian, psol

C Initial parameters:
        mf       = 21  !21/29 SPIGM/direct method
        rtol     = rtol1
        atol     = atol1
        itol     = 4					 !whether atol,rtol are scalers or arrays
        itask    = 1
        istate   = 1
        iopt     = 1
        neqz     = neq
        lrw      = neqz*nss+721+32*neqz  !Declared size of RWORK (61+17*n+LWP)
        liw      = 30+5*neqz             !Declared size of IWORK (30+LIWP)
		IWORK(1) = neqz*nss              !LWP - size of real array for precond. 
		IWORK(2) = neqz*5                !LIWP - size of integer array for precond.
		IWORK(3) = 1                     !JPRE=0,1,2,3
		IWORK(4) = 1                     !JACFLG=0,1
		IWORK(6) = 20000                 !MXSTEP
        t = 0.0d0
		tout = tfirst
        tstep= 1.0D+01**(dlog10(tlast/tfirst) / (nt-1.0d0))
		grnfrac = 6e-12
C Call DVODEPK:
      do i = 1, nt

C Set the time-dependent reaction rates

		timestep = i
		timedep = 1
	    DO j = 1, nre
		  IF (rtype(j) .LT. 0) THEN

				ak(j) = calcrate(alpha(j), beta(j), gamma(j), rtype(j), 
     &				r1(j),r2(j),grnfrac, timedep, dtaudt,tout)
c			if (ak(j) .LE. 1.0d-100) then
c				write(98,*) j, rtype(j), r1(j), ak(j)
c			endif
		  END IF
	    END DO
  	  	     
        times(i) = tout 
        write(*,10) 'zone = ', zone, 't1 = ',tout / 3.155D+07
   10   FORMAT(1X,A7,I2,2X,A5,1pE10.3)

 20		call dvodpk(Fcn,neqz,yy,t,tout,itol,rtol,atol,itask,
     &     istate,iopt,rwork,lrw,iwork,liw,Jacobian,psol,mf,rpar,ipar)

        tout = tout * tstep

C An error occured during the integration:
C If something was wrong, print warning:
C istate   == 2  if "dvode" was successful, negative otherwise,
C            -1 means excess work done on this call. (Perhaps wrong mf.)
C            -2 means excess accuracy requested. (Tolerances too small.)
C            -3 means illegal input detected. (See printed message.)
C            -4 means repeated error test failures. (Check all input.)
C            -5 means repeated convergence failures. (Perhaps bad
C               Jacobian supplied or wrong choice of mf or tolerances.)
C            -6 means error weight became zero during problem. (Solution
C               component i vanished, and atol or atol(i) = 0.)
C
c      if (istate.lt.0) then
c         write (*,*) 'istate = ',istate
c      end if
      if (istate.eq.-1) then
         write(*,*) 'excessive amount of work was done, continue...'
c         write(*,*) 'number of steps will be increased by 1.25'
c         iwork(6) = iwork(6) * 1.25
c         if (iwork(6) .gt. 1e6) then
c			write(*,*) 'Too many steps (> 1e6)'
c			stop
c		 endif
         istate = 2
	   tout = tout / tstep
         goto 20
      else if (istate.eq.-2) then
         write(*,*) 'requested accuracy too much,'
         write(*,*) 'relative tolerance will increase by 2'
         do j = 1,nspec
	     rtol(j) = rtol(j) * 2.0D0
	     if (rtol(j).gt.1.0D-5) stop
	 enddo
c         atol = 2.23D-16
         istate = 3
	   tout = tout / tstep
         goto 20
      else if (istate.eq.-3) then
         write(*,*) 'illegal input or infinite loop of calls, stop'
         stop
      else if (istate.eq.-4) then
         write(*,*) 'repeated error test failures (singularity?)'
c		 Write out species with largest error associated with them
		 write(*,*) IWORK(16), s(IWORK(16))
         istate = 2
	   tout = tout / tstep
         goto 20
      else if (istate.eq.-5) then
c         write(*,*) 'repeated convergence test failures,'
c         write(*,*) 'perhaps Jacobian is not accurate enough'
         istate = 2
	   tout = tout / tstep
         goto 20 
      else if (istate.eq.-6) then
c         write(*,*) 'a solution component "i" is vanished, but'
c         write(*,*) 'pure absolute error control ATOL = 0 was requested'
         istate = 2
	   tout = tout / tstep
         goto 20
      end if


C TEST - Rate Tests
C j = species index, k = reacs index
		IF (ratetest) THEN
			do j=1,Nratetest
				open(unit=30+j, file=fnameRT(j), status='old', 
     +				access='append')
				write(30+j,100,advance='no') zone,tout/tstep/3.155D+07
				do k=1,NRTreacs(j)
					write(30+j, 101, advance='no') 
     +					RTreacs(j,k), RTrates(j, k)
				enddo
				write(30+j, *)

				close(30+j)
			enddo
		END IF
 100	FORMAT(1X,I6,3X,1pe8.2)
 101	FORMAT(2X,I6,2X,1pe9.2)

C Save computed abundances at i-th time step:
        do k = 1, neq
c  add a break on the populations to not allow them to go too low
			IF (yy(k) .LT. MINABUN) THEN
				yy(k) = MINABUN
			END IF
			abundances(k,i) = yy(k)
        end do

c end of dvode loop
      end do

C Exit:
      return
      end

