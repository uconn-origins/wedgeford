C..............................................................................
C
C This function creates the Jacobian used in solving the ODEs, called
C by the ODE solver
C
C..............................................................................
      subroutine Jacobian(F,n,t,y,ysv,rewt,fty,V,hrl1,JacS,IWP,IER,
     1     RPAR, IPAR)
      implicit none


c Global variable:
      integer n, IPAR, IWP, IER
      real*8 y, ysv, t, RPAR, rewt, fty, v, hrl1, JacS
      dimension y(*), IPAR(*), RPAR(*), ysv(*), rewt(*), fty(*), 
     1   V(*), JacS(*), IWP(*)

C Common blocks:
      INCLUDE "Fcn.h"
	
C Local variables:
	integer j,k,jr1, jr2, jp1, jp2, jp3, jp4, jp5, 
     1    ir, iter, kk, jpvt
	real*8 term1, term2, U, Jac
	DIMENSION Jac(nspec,nspec),jpvt(8*nspec),ir(nspec*nss/2)

C External function:
      EXTERNAL F

C Construct JAC "on the fly":
      kk = 0

C Initialization:
	DO j = 1, ns
        DO k = 1, ns
          Jac(j,k) = 0.0d0
	  ENDDO
	ENDDO

C Loop by reactions:
      DO j = 1, nre
	  jr1 = ir1(j)
	  jr2 = ir2(j)
	  jp1 = ip1(j)
	  jp2 = ip2(j)
	  jp3 = ip3(j)
	  jp4 = ip4(j)
	  jp5 = ip5(j)	  

C No second reactant:
        IF (jr2.EQ.0) THEN
	    term1 = ak(j)
	    Jac(jr1,jr1) = Jac(jr1,jr1) - term1
	    Jac(jp1,jr1) = Jac(jp1,jr1) + term1
	    if (jp2.ne.0) Jac(jp2,jr1) = Jac(jp2,jr1) + term1
	    if (jp3.ne.0) Jac(jp3,jr1) = Jac(jp3,jr1) + term1
	    if (jp4.ne.0) Jac(jp4,jr1) = Jac(jp4,jr1) + term1
	    if (jp5.ne.0) Jac(jp5,jr1) = Jac(jp5,jr1) + term1	    
	  ELSE
C Reactions with second reactant:
          term1 = ak(j)*y(jr1)
          term2 = ak(j)*y(jr2) 
          Jac(jr1,jr1) = Jac(jr1,jr1) - term2
          Jac(jr1,jr2) = Jac(jr1,jr2) - term1
	    Jac(jr2,jr1) = Jac(jr2,jr1) - term2
          Jac(jr2,jr2) = Jac(jr2,jr2) - term1
	    Jac(jp1,jr1) = Jac(jp1,jr1) + term2
          Jac(jp1,jr2) = Jac(jp1,jr2) + term1
	    if (jp2.ne.0) Jac(jp2,jr1) = Jac(jp2,jr1) + term2
          if (jp2.ne.0) Jac(jp2,jr2) = Jac(jp2,jr2) + term1
	    if (jp3.ne.0) Jac(jp3,jr1) = Jac(jp3,jr1) + term2
          if (jp3.ne.0) Jac(jp3,jr2) = Jac(jp3,jr2) + term1
	    if (jp4.ne.0) Jac(jp4,jr1) = Jac(jp4,jr1) + term2
          if (jp4.ne.0) Jac(jp4,jr2) = Jac(jp4,jr2) + term1
	    if (jp5.ne.0) Jac(jp5,jr1) = Jac(jp5,jr1) + term2
          if (jp5.ne.0) Jac(jp5,jr2) = Jac(jp5,jr2) + term1          
	  END IF

C End loop by reactions:
      END DO

C Transform to the form I - Jac*hrl1 and initialize non-zero elements:
	DO j = 1, ns
	  DO k = 1, ns	
	    Jac(j,k)= -hrl1*Jac(j,k)
          if (k.eq.j) Jac(j,j) = Jac(j,j) + 1.0d0	    	   	     
	    if (Jac(j,k).ne.0.0d0) then
	      kk = kk + 1
	      JacS(kk)= Jac(j,k) 
            ir(kk)  = j
            IPAR(kk)  = k
	    endif
	  ENDDO
	ENDDO

	nJ = kk
	iter = iter + 1

C Factorize the matrix:
	  U = 1.0d0
       
	  call MA28AD(n,nJ,JacS,nJ*4,ir,nJ*2,IPAR,U,IWP,jpvt,RPAR,IER)

      return
      end 
C
C Solution of the system Ax=b:
C
      SUBROUTINE PSOL(N, T, Y, FTY, WK, HRL1, WP, IWP, V, LR, IER,
     1                  RPAR, IPAR)

C Global parameters:
      INTEGER IWP, IER, IPAR, N, LR
      DOUBLE PRECISION T, Y, FTY, WK, HRL1, WP, V, RPAR
      DIMENSION Y(*), FTY(*), WK(*), WP(*), IWP(*), V(*), 
     1	RPAR(*), IPAR(*)

C Common blocks:
      INCLUDE "Fcn.h"

C Direct solution of the system Ax=b:
	CALL MA28CD(n,WP,nJ*4,IPAR,IWP,V,RPAR,1)

      RETURN
      END 