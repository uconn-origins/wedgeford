C..............................................................................
C
C This subroutine read a file with chemical reactions.
C
C..............................................................................
C
C Version 1.2 (11/03/2002)
C
C..............................................................................
C
C Input parameter(s):
C
C input      == a name of a file with chemical reactions,
C
C..............................................................................
C
C Common block(s):
C
C..............................................................................
C
C BL1:
C
C ns           == amount of species,
C
C s(ns)        == a species set,
C
C nre          == amount of reactions in the input file,
C
C index(nre)    == indices of the corresponding chemical reactions,
C
C r1(nre)       == names of the first reactants,
C
C ir1(nre)      == array of indices of 'r1' species in extended species set 's',
C
C r2(nre)       == names of the second reactants,
C
C ir2(nre)      == array of indices of 'r2' species in extended species set 's',
C
C p1(nre)       == names of the first products, 
C
C ip1(nre)      == array of indices of 'p1' species in extended species set 's',
C
C p2(nre)       == names of the second products, 
C
C ip2(nre)      == array of indices of 'p2' species in extended species set 's',
C
C p3(nre)       == names of the third products, 
C
C ip3(nre)      == array of indices of 'p3' species in extended species set 's',
C
C p4(nre)       == names of the fourth products, 
C
C ip4(nre)      == array of indices of 'p4' species in extended species set 's',
C
C p5(nr)        == names of the 5-th products, 
C
C ip5(nr)       == array of indices of 'p5' species in extended species set 's',
C
C alpha(nre)    == first components of the rate coeffs.,
C
C beta(nre)     == second components of the rate coeffs.,
C
C gamma(nre)    == third components of the rate coeffs.
C
C..............................................................................
C
C Used subroutines(s) (alphabetically): ispecies, len_tri2
C
C..............................................................................
      SUBROUTINE readr(input)
      IMPLICIT NONE

C Global variable(s):
      CHARACTER*80 input

C Local variable(s):
      INTEGER i, nlen

C Initialization of common blocks:
      INCLUDE "Fcn.h"

C Format(s):
 100  format (1X,I4,1X,2(A10),10X,5(A10),E8.2,1X,F5.2,1X,F8.1,1X,I2)

C Open input file:
      CALL len_tri2(input,80,nlen)
      OPEN (unit=10,file=input(1:nlen),status='old',
     &      access='sequential')

C..............................................................................
C Read chemical network data from input:
C..............................................................................
      READ (10, *) nre
      
      DO i = 1, nre
      
      READ (10,100) index(i),r1(i),r2(i),p1(i),p2(i),p3(i),
     &        p4(i), p5(i), alpha(i), beta(i), gamma(i), rtype(i)

C Search positions of species 'r1, r2, p1, p2, p3, p4' in 's':
         CALL ispecies(r1(i),ns,s,ir1(i))
         CALL ispecies(r2(i),ns,s,ir2(i))
         CALL ispecies(p1(i),ns,s,ip1(i))
         CALL ispecies(p2(i),ns,s,ip2(i))
         CALL ispecies(p3(i),ns,s,ip3(i))
         CALL ispecies(p4(i),ns,s,ip4(i))   
         CALL ispecies(p5(i),ns,s,ip5(i)) 
         
      END DO
      
C Close all files:
      CLOSE (10)
      
C Exit:
      RETURN
      END
C..............................................................................
C
C This subroutine reads the names of chemical species from "specs".
C
C..............................................................................
C
C Version 1.2 (07/12/2004)
C
C..............................................................................
C
C Input parameter(s): 
C
C input   == a name of a file contains the names of the species,
C
C..............................................................................
C
C Output parameter(s):
C
C y(1:ny) == the names of species,
C
C..............................................................................
C
C Global parameter(s):
C
C nspec       == maximal amount of chemical species to be considered,
C
C..............................................................................
      subroutine reads(input)
      implicit none

C Global variable(s):
      character*80 input  ! the name of the input file,

C Local variable(s):
      integer i  ! counters for loops
      integer nlen !, first  ! a length of a string, the beginning of a string,

C Initialization of common blocks:
      INCLUDE "Fcn.h"
C
C Open input file(s):
C
      call len_tri2(input,80,nlen)
      open (unit=07,file=input(1:nlen),status='old',access='sequential')
      rewind 07
C
C Read 'y' from the file in the next loop:
C
      read (07,*) ns
      
      do i=1,ns
      
        read (07,'(a10)') s(i)
        
      enddo  
      
C Close all file(s):
      close (07)      

C Exit:
      return
      end