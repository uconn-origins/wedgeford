
C..............................................................................
C
C This subroutine read a file with chemical reactions.
C
C..............................................................................
C
C Input parameter(s):
C
C input      == a name of a file with chemical reactions,
C
C..............................................................................
C
C Used subroutines(s) (alphabetically): ispecies
C
C..............................................................................
      SUBROUTINE readr(input)
      IMPLICIT NONE

C Global variable(s):
      CHARACTER*80 input

C Local variable(s):
      INTEGER i, nlen, nlines, fend
	  CHARACTER*150 line

C Initialization of common blocks:
      INCLUDE "Fcn.h"

C Format(s):
 100  format (I5,1X,2(A13),4(A13),2X,E8.2,1X,E9.2,1X,E9.2,1X,I3)
C %4d %-12s %-12s %-12s %-12s %-12s %-12s   %8.2e %9.2e %9.2e %3d
 
C Open input file:
      OPEN (unit=10,file=trim(input),status='old')

C..............................................................................
C Read chemical network data from input:
C..............................................................................
	  fend = 0
	  nlines = -1
	  nre = 0
	  i = 0

      DO WHILE (fend .EQ. 0)
		READ(10, '(A150)', IOSTAT=fend) line
		nlines = nlines + 1
C Skip commented out lines (indicated by #)
		IF (line(1:1) .NE. '#') THEN
c		    print *, line
			nre = nre + 1
			i = i + 1
C Read in the reaction
			READ (line,100) reacindex(i),r1(i),r2(i),p1(i),p2(i),p3(i),
     &				p4(i), alpha(i), beta(i), gamma(i), rtype(i)
c	        print *,reacindex(i),r1(i),r2(i),alpha(i),beta(i),gamma(i),rtype(i)
		ENDIF 

C Search positions of species 'r1, r2, p1, p2, p3, p4' in 's':
         CALL ispecies(r1(i),ns,s,ir1(i))
         CALL ispecies(r2(i),ns,s,ir2(i))
         CALL ispecies(p1(i),ns,s,ip1(i))
         CALL ispecies(p2(i),ns,s,ip2(i))
         CALL ispecies(p3(i),ns,s,ip3(i))
         CALL ispecies(p4(i),ns,s,ip4(i))   
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
C Input parameter(s): 
C
C input   == a name of a file contains the names of the species,
C
C..............................................................................
C
C Output parameter(s):
C
C y(ny)   == the names of species,
C
C..............................................................................
C
C Global parameter(s):
C
C nspec   == maximal amount of chemical species to be considered,
C
C..............................................................................
      subroutine reads(input)
      implicit none

C Global variable(s):
      character*80 input  ! the name of the input file,

C Local variable(s):
      integer i, j  ! counters for loops
      integer nlen !, first  ! a length of a string, the beginning of a string,
	  integer fend
	  CHARACTER*150 line

C Initialization of common blocks:
      INCLUDE "Fcn.h"
C
C Open input file(s):
C
      open (unit=07,file=trim(input),status='old')

C
C Read 'y' from the file in the next loop:
C
	  fend = 0
	  ns = 0
	  i = 0
	  j = 0


      DO WHILE (fend .EQ. 0)
		READ(07, '(A150)', IOSTAT=fend) line
C Skip commented out lines (indicated by #)
		IF (line(1:1) .NE. '#') THEN
			ns = ns + 1
			i = i + 1

C Read in the species
			READ (line,'(A13)') s(i)

C If it is a grain species, store species number
			IF ('(gr)' .EQ. s(i)(len_trim(s(i))-3:)) THEN
				grspecs(j) = i
				j = j + 1
			ENDIF

		ENDIF 
	  END DO

C Marker for end of grain species in array. 
	  grspecs(j) = -1

C Close all file(s):
      close (07)      

C Exit:
      return
      end