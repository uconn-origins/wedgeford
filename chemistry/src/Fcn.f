C..............................................................................
C
C This function creates the ODEs to be solved, called by the ode solver function
C
C..............................................................................
C
C Input parameter(s):
C
C n			== maximum number of species
C
C t			== time value
C 
C y(n)		== species abundance
C 
C ydot(n)	== species ODE 
C 
C RPAR		== ODE solver return variable
C
C IPAR		== ODE solver return variable
C
C..............................................................................

      SUBROUTINE Fcn(n,t,y,ydot,RPAR,IPAR)
      IMPLICIT NONE
	                                
C Global variable(s):
      INTEGER n, IPAR
      REAL*8 t, y, ydot, RPAR                                       
      DIMENSION y(n),ydot(n), RPAR(*), IPAR(*) 

C Common blocks:
      INCLUDE "Fcn.h"
	                                           
C Local variable(s):
      INTEGER i, jr1, jr2, jp1, jp2, jp3, jp4, jp5
      REAL*8 term
      CHARACTER*13 rr1, rr2

C TEST Variables - Max Rates
	  INTEGER itest(Nratetest), j
	  REAL*8 totdes, totfor, limit
	  totdes = 0
	  totfor = 0
	  limit = 1e-20
	  DO j=1,Nratetest
		NRTreacs(j) = 0
		CALL ispecies(RTspec(j), ns, s, itest(j))
	  END DO
	  
C Initialization of the array:
	  ydot = 0.0D0

C Construct ODE system "on the fly":
C dy/dt:
      DO i = 1, nre
        
        rr1 = r1(i)
        rr2 = r2(i)          
	    jr1 = ir1(i)
	    jr2 = ir2(i)
	    jp1 = ip1(i)
	    jp2 = ip2(i)
	    jp3 = ip3(i)
	    jp4 = ip4(i)
	    jp5 = ip5(i)

C Destruction of species: 
      IF (jr2.EQ.0) THEN     
	    term = ak(i)*y(jr1) 
        ydot(jr1)=ydot(jr1)-term      
      ELSE
	    term = ak(i)*y(jr1)*y(jr2)
        ydot(jr1)=ydot(jr1)-term
        ydot(jr2)=ydot(jr2)-term        
      END IF

C TEST - Max Rates
	  IF (ratetest) THEN
		  DO j=1,Nratetest
			IF (term .LT. limit) EXIT
			IF (jr1.EQ.itest(j) .OR. jr2.EQ.itest(j)) THEN
				NRTreacs(j) = NRTreacs(j) + 1
				RTreacs(j,NRTreacs(j)) = i
				RTrates(j,NRTreacs(j)) = -term
			ENDIF
		  END DO
	  END IF

C Formation of species:      
	  ydot(jp1)=ydot(jp1)+term
      IF (jp2.NE.0) ydot(jp2)=ydot(jp2)+term
      IF (jp3.NE.0) ydot(jp3)=ydot(jp3)+term
      IF (jp4.NE.0) ydot(jp4)=ydot(jp4)+term
      IF (jp5.NE.0) ydot(jp5)=ydot(jp5)+term

C TEST - Max Rates
	  IF (ratetest) THEN
		  DO j=1,Nratetest 
			IF (term .LT. limit) EXIT
			IF (jp1.EQ.itest(j) .OR. jp2.EQ.itest(j) .OR.
     +			jp3.EQ.itest(j) .OR. jp4.EQ.itest(j)) THEN

				NRTreacs(j) = NRTreacs(j) + 1
				RTreacs(j,NRTreacs(j)) = i
				RTrates(j,NRTreacs(j)) = term
			ENDIF
		  END DO
	  END IF

      END DO

C Exit:
      RETURN
      END
