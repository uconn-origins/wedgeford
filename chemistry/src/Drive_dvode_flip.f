
C------------------------------------------------------------------------------
C
C "Driver" subroutine to the stiff SODE solver "DVODEPK".
C
C runs the chemistry for a single time step at once
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
C
C
C 
C
C------------------------------------------------------------------------------
C Used subroutine(s) and function(s) (alphabetically):
C
C Function(s): calcrate, Fcn, Jacobian
C Subroutine(s): dvodpk
C
C------------------------------------------------------------------------------
      subroutine run_chemistry_dt(neq,yy,time_j,time_jp1,rtol0,atol0,abunzone)
      implicit none

C Initialization of common blocks:
      INCLUDE "Fcn.h"
      INCLUDE "dvode.h"
      INCLUDE "environ.h"
      INCLUDE "BindingE.h"

C Global variable(s):
      integer neq
      real*8 yy, time_j, time_jp1, rtol0(nspec), atol0(nspec),abunzone
      dimension yy(nspec), abunzone(nspec,nzone)
      
C Local variable(s):
      DOUBLE PRECISION calcrate

      INTEGER n_mols,dummyy
      PARAMETER (n_mols = 4)
      real grnfrac
      integer mf, itol, itask, istate, iopt, lrw, liw, IPAR, i, k,
     1   neqz, timedep, j, ISAV
      real*8 RPAR, atol(nspec), rtol(nspec), t0, t1, RSAV
      dimension RPAR(nspec), IPAR(nspec*nss), ISAV(100), RSAV(100)
      double precision dtaudt(67,3),nomo_dtaudt(67,3),Ngrsp(n_mols)
      double precision n_dtaudt(67,3),old_Ngrsp(n_mols),old_dtaudt(67,3)
      double precision diffabun
      INTEGER imol
      CHARACTER*13 species
      CHARACTER*13 specarr(n_mols)
      DOUBLE PRECISION numtot, nr1, Numr1(n_mols)    
      INTEGER r1index, grindex

C External routines:
      EXTERNAL Fcn, Jacobian, psol

C Initial parameters:
        mf       = 21  !21/29 SPIGM/direct method
        rtol     = rtol0
        atol     = atol0
        itol     = 4                          !whether atol,rtol are scalers or arrays
        itask    = 1
       
        iopt     = 1
        neqz     = neq
        lrw      = neqz*nss+721+32*neqz  !Declared size of RWORK (61+17*n+LWP)
        liw      = 30+5*neqz             !Declared size of IWORK (30+LIWP)
          
        grnfrac = 6e-12
        t0      = time_j
        t1      = time_jp1
        
        IF (t0 .LT. 1) THEN
          istate   = 1
        ELSE
          RSAV = rsavzone(zone,:)
          ISAV = isavzone(zone,:)
          CALL dvksrc(RSAV,ISAV,2)
          istate   = 2
        END IF
        
        IWORK(1) = neqz*nss              !LWP - size of real array for precond.
        IWORK(2) = neqz*5                !LIWP - size of integer array for precond.
        IWORK(3) = 1                     !JPRE=0,1,2,3
        IWORK(4) = 1                     !JACFLG=0,1
        IWORK(6) = 20000                 !MXSTEP
        IWORK(7) = 1                     !maximum number of error messages when iterative step is too small
        

        DO j = 1, nre
          timedep = 1
          IF (rtype(j) .LT. 0) THEN
            ak(j) = calcrate(alpha(j), beta(j), gamma(j), rtype(j),r1(j),r2(j),grnfrac, timedep, dtaudt,t1)
          END IF
        END DO
        DO k = 1, neq
           yy(k) = abunzone(k,zone)
c           atol(k) = MAX(atol(k)*yy(k), MINABUN)
        END DO
        
        WRITE(*,10) 'zone = ', zone, 't0 = ',t0 / 3.155D+07
10      FORMAT(1X,A7,I2,2X,A5,1pE10.3)

20      CALL dvodpk(Fcn,neqz,yy,t0,t1,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,Jacobian,psol,mf,rpar,ipar)


        
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

      ELSE IF (istate.eq.-5) THEN
         print *, 'repeated convergence test failures, problem species:'
         print *, IWORK(16), s(IWORK(16)) !write out species that are gumming up the works
         istate = 2
         j = IWORK(16)
         atol(j) = atol(j) * 2.0D0
         IF (atol(j) .gt. yy(j) * 2.0D0) THEN
             write(*,*) 'ERROR in abundance on the order of abundance :/'
         END IF
c         goto 20
      ELSE IF (istate.eq.-6) THEN
         print *, 'a solution component "i" is vanished, but'
         print *, 'pure absolute error control ATOL = 0 was requested'
         do j = 1,nspec
             atol(j) = MINABUN*1e-2 ! reset absolute tolerances
         end do
         istate = 2
c         goto 20
      END IF !end of error check
      
C     save common block parameters for the zone
      CALL dvksrc(RSAV,ISAV,1)
      IF (istate .eq. -3) THEN !if this zone doesn't work, remove it from running
        print *, 'illegal input or infinite loop of calls, removing zone from calculation'
        runzone(zone) = .FALSE.
        GO TO 40
      ELSE IF (istate .eq. -2) THEN !if the tolerances are too low, use the scaling factor to adjust globally
        DO j = 1, ns
         atol1(j) = RWORK(14) * atol1(j)
         rtol1(j) = RWORK(14) * rtol1(j)
        END DO
      ELSE IF (istate .eq. -4) THEN 
        
        
        rsavzone(zone,:) = RSAV
        isavzone(zone,:) = ISAV


      

C TEST - Rate Tests
C j = species index, k = reacs index
      IF (ratetest) THEN
        DO j=1,Nratetest
          OPEN(unit=30+j, file=fnameRT(j), status='old',access='append')
           WRITE(30+j,100,advance='no') zone, t1/3.155D+07
            DO k=1,NRTreacs(j)
              WRITE(30+j, 101, advance='no') RTreacs(j,k), RTrates(j, k)
            END DO
           WRITE(30+j, *)
          CLOSE(30+j)
        END DO
      END IF !end of rate test
100   FORMAT(1X,I6,3X,1pe8.2)
101   FORMAT(2X,I6,2X,1pe9.2)

C Save computed abundances in the zone:
      DO j = 1, neq
c  add a break on the populations to not allow them to go too low
      IF (yy(j) .LT. MINABUN) THEN
        yy(j) = MINABUN
      END IF ! end of abundances check
      END DO 
40    abunzone(:,zone) = yy
      RETURN
      END
