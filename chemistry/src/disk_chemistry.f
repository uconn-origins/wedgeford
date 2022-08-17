C..............................................................................
C
C Program to model the chemical evolution under specified initial and physical
c conditions
C
C..............................................................................
C
C copyright (c) 2002-3, Dmitry Semenov, AIU Jena
C copyright (c) 2004-7, Dmitry Semenov, MPIA Heidelberg
C copyright (c) 2008-11, Jeffrey Fogel, University of Michigan
C
C..............................................................................
C  Modifications added by Jeffrey Fogel, University of Michigan
C  Grain reactions, X-ray photoionization, Self-shielding, uv field,
C  Lyman alpha radiation, photodesorption, time dependent reactions
C..............................................................................
C
C Input parameter(s):  file: 0io.inp
C
C specs       == the name of a file contains names of chemical elements
C                involved in a set of chemical reactions,
C
C rates       == the name of a file contains a system of chemical
C                reactions for species from 'specs',
C
C uvfile       == the name of a file contains the UV field, as calculated
C				  from Nuria Calvet's codes.
C
C out         == a name of output file contains calculated abundances,
C
C
C input parameter(s):  file: 2times.inp
C
C tlast       == ending time of the evolution of the system [years],
C
C tfirst      == first time step [years],
C
C nstep       == number of time steps.
C
C
C input parameter(s):  file: 3abunds.inp
C
C nfrc        == number of initially present species
C
C spec0(nfrc) == Species name
C
C frc(nfrc)   == and their abundances (relative to H).
C
C
C input parameter(s):  file: 4toleran.inp
C
C relerr      == relative error of the calculations,
C
C abserr      == absolute error of the calculations
C
C
C input parameter(s):  file: 5flags.inp
C
C
C g0_100	  == G0 value, normalized at 100 AU
C
C g0_CO		  == G0 value for CO (calculated ignoring Lyman alpha radiation)
C
C Lxray		  == X-ray luminosity [ergs/s]
C
C RadNuc	  == Radionuclide Ionization Rate [1/s]
C
C fg          == dust to gas mass ratio (normal is 100)
C
C freezeeffic == parameter to approximate freeze out effects in grain growth. Normal (0.1 mu grains) is 1.
C
C CRdesorp	  == turn on/off cosmic ray desorption
C
C CRionization == turn on/off cosmic ray ionization
C
C photodesorp == turn on/off UV photodesorption
C
C LyAphotodesorp	== turn on/off lyman alpha photodesorption
C
C thermaldesorp		== turn on/off thermal desorption
C
C include_LyA == include lyman alpha radiation or remove it from UV field
C
C xraydust		== turn on/off dust-dependent x-ray opacities
C
C incl_radionuc == turn on/off radionuclide ionization.  Must include external file w/ RN calculations.
C
C ratetest	  == print out rates of given species or not
C
C RTspec(maxratetest)	== species to print rates of
C
C shieldtest  == print out self-shielding test files or not
C
C
C input parameter(s):  file: 1environ.inp
C
C Nr          == number of radial grid points [SET TO 1],
C
C Nz          == number of vertical grid points,
C
C Rs(Nz)      == radius (AU),
C
C rho(Nz)     == density [g/cm^3],
C
C ngr(Nz)     == dust number density [1/cm^3],
C
C Tg(Nz)      == gas temperature [K],
C
C Td(Nz)      == dust temperature [K],
C
C zAU(Nz)	  == Height above the midplane [AU]
C
C zcm(Nz)	  == Height from disk surface, surface = 1.0e10 cm [cm]
C
C Nrz(Nz)     == Radial column density [1/cm^2]
C
C ZetaCR(Nz)  == cosmic ray ionization rate.
C
C
C..............................................................................
C
C Output file(s):
C
C [out].out (set in '0io.inp')
C
C..............................................................................
C
C Global parameter(s):
C
C nreac       == maximal amount of chemical reactions to be read,
C
C nspec       == maximal amount of chemical species to be considered,
C
C ntime       == maximal amount of taken time steps,
C
C..............................................................................
C
C Common block(s):
C
C Look in *.h files
C
C..............................................................................
C
C Important variable(s):
C
C s(1:ns)     == an array with the names of species,
C
C y(1:ns)     == an array with calculated number densities of species,
C
C gdens       == number density of hydrogen nuclei,
C
C..............................................................................
C
C Used subroutines(s) (alphabetically): calcCRatten, calcrate, ini_abunds,
C		ini_abunds2D, readr, reads, readuv, run_chemistry
C
C..............................................................................
      PROGRAM disk_chemistry
      IMPLICIT NONE

C Initialization of common blocks:
      INCLUDE "Fcn.h"
      INCLUDE "constants.h"
      INCLUDE "environ.h"
      INCLUDE "rates.h"
      INCLUDE "BindingE.h"

C Global variable(s):
      INTEGER nfrc, nstep, istart, iend, fend
      INTEGER lastzer
      REAL*8 tlast, tfirst, relerr(nspec), abserr(nspec), frc,
     1  yy, gdens
      CHARACTER*13 spec0, rr2
      CHARACTER*80 specs, rates, out, ver, cmgtm, uvfile, xrayfile
      CHARACTER*80 radionucfile
      CHARACTER*80 isrffile
      CHARACTER*80 xrayfilesb
      CHARACTER*80 dustn
      CHARACTER*80 abun2dfile
      CHARACTER*80 outabunfile
      CHARACTER*28 file2d,newfile2d
      DIMENSION frc(nspec2), spec0(nspec2), yy(nspec)
      real grnfrac

C External Functions
      DOUBLE PRECISION xrayion_rate, calcrate, dtaudtemp, polint

C Local variable(s):
      DOUBLE PRECISION errval, lastyr
      INTEGER i, j, l, nlen, nait, ipon, iana,
     1  iff, lx, is, li, lt, timedep, pos
      CHARACTER*5 ait, siend, eqlsign
      CHARACTER*160 line,testfi
      CHARACTER*13 errspec, lts
      CHARACTER*20 flagname
      CHARACTER*140 flagval
      CHARACTER*100 iofile, envfile
      CHARACTER*80 preadd, radstring, postadd, text1, text2, endenv
      INTEGER intflagval, dot, errspecindex
c	  real dtaudt(9000,3)
      integer nln,indv
      INTEGER jg ! iterates final zone species for lower error, temporary
      LOGICAL sameness


C Current version of the code:
      PARAMETER (ver = ' version 2.7.3 (02/17/2015)')

C Initialization of variables for 'timer':
      cmgtm = 'disk_chemistry,' //ver
      print *, cmgtm

C..............................................................................
C Read input data:
C 1) input & output file names:
C..............................................................................
        print *, 'start'
	  call getarg(1,iofile)
      open (unit=01,file=iofile,status='old',access='sequential')
        read (01,*)
        read (01,'(A40)') specs
        read (01,'(A40)') rates
        read (01,'(A80)') uvfile
            uvfile = uvfile( :index(uvfile,'#')-1)
C            print *, uvfile
        read (01,'(A80)') xrayfile
            xrayfile = xrayfile( :index(xrayfile,'#')-1)
C            print *, xrayfile
        read (01,'(A80)') isrffile
            isrffile = isrffile( :index(isrffile,'#')-1)
C            print *, isrffile
C        read (01,'(A80)') xrayfilesb
C            xrayfilesb = xrayfilesb( :index(xrayfilesb,'#')-1)
        read (01,'(A80)') radionucfile
            radionucfile = radionucfile( :index(radionucfile,'#')-1)
            print *, radionucfile
        read (01,'(A80)') abun2dfile
            abun2dfile = abun2dfile( :index(abun2dfile,'#')-1)
C            print *, abun2dfile
      close (01)
      print *, "abun2dfile!!!!!!!! ",abun2dfile
C Read species set from 'specs':
        CALL reads (specs)
        print *, "read specs"
        print *, rates
        CALL readr (rates)
      print *,"read specs and rates"
C..............................................................................
C 2) Time evolution:
C..............................................................................
      print *, '2times'
      open (unit=01,file='2times.inp',status='old',access='sequential')
        read (01,*)
        read (01,*) tlast
          tlast = tlast * year
        read (01,*) tfirst
          tfirst = tfirst * year
        read (01,*) nstep

      close (01)
      print *, "read 2times"

C..............................................................................
C 3) Initial abundances:
C..............................................................................
      if (tfirst/year .LT. 9.000D-00) then
      	open (unit=01,file='3abunds.inp',status='old')
		fend = 0
		nfrc = 0
		do while (fend .EQ. 0)
			read(01,'(A160)',iostat=fend) line
			if (line(1:1) .NE. '#') then
				nfrc = nfrc + 1
				read(line,*) spec0(nfrc), frc(nfrc)
			endif
		enddo

      	close(01)
      	endif
        print *, "read 3abunds"

C..............................................................................
C 4) Tolerance:
C..............................................................................
c       write(*,*) 'LIC here 0a'
	   open (unit=01,file='4toleran.inp',status='old')
        read (01,*)
        read (01,*) errval
		do i=1,nspec
			relerr(i) = errval
		enddo

        read (01,*) errval
		do i=1,nspec
			abserr(i) = errval
		enddo

c       set species-specific tolerances
		fend = 0
		do while (fend .EQ. 0)
			read(01,'(A160)',iostat=fend) line
			if (line(1:1) .EQ. '#') cycle

			read(line,*) flagname, errspec, errval
			call ispecies(errspec, ns, s, errspecindex)
			if (trim(flagname) .EQ. 'rel') then
				relerr(errspecindex) = errval
			else if (trim(flagname) .EQ. 'abs') then
				abserr(errspecindex) = errval
			else
				write(*,*) 'Illegal value in 4toleran.inp'
				stop
			endif
		enddo

      close (01)
      print *, "read 4tol"

      !  Initialize freezeeffic (read in below)
	  freezeeffic = 1.0

C..............................................................................
C 5) Flags:
C..............................................................................
		open (unit=01,file='5flags.inp',status='old')
		fend = 0
		do while (fend .EQ. 0)
			read(01,'(A160)',iostat=fend) line
			if (line(1:1) .NE. '#') then
				read(line,*) flagname, eqlsign, flagval
c 100			format(A13,3X,A140)

C Parameters
			if (trim(flagname).EQ.'fg') then
				read(flagval, '(E8.2)') fg
				do i=1,nfrc
					if (trim(spec0(i)).NE.'H2'.AND.trim(spec0(i)).NE.'He') then
						frc(i) = frc(i) * 100.0 / fg
					endif
				enddo

			else if (trim(flagname).EQ.'freezeeffic') then
				read(flagval, '(E8.2)') freezeeffic

			else if (trim(flagname).EQ.'epsilon') then
				read(flagval, '(E8.2)') ndust

C Flags
			else if (trim(flagname).EQ.'CRdesorp') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					CRdesorp = .FALSE.
				else
					CRdesorp = .TRUE.
				endif
			else if (trim(flagname).EQ.'CRionization') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					CRionization = .FALSE.
				else
					CRionization = .TRUE.
				endif
			else if (trim(flagname).EQ.'photodesorp') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					photodesorp = .FALSE.
				else
					photodesorp = .TRUE.
				endif
			else if (trim(flagname).EQ.'LyAphotodesorp') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					LyAphotodesorp = .FALSE.
				else
					LyAphotodesorp = .TRUE.
				endif
			else if (trim(flagname).EQ.'thermaldesorp') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					thermaldesorp = .FALSE.
				else
					thermaldesorp = .TRUE.
				endif
			else if (trim(flagname).EQ.'include_lya') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					include_LyA = .FALSE.
				else
					include_LyA = .TRUE.
				endif
			else if (trim(flagname).EQ.'xraydust') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					xraydust = .FALSE.
				else
					xraydust = .TRUE.
				endif
			else if (trim(flagname).EQ.'incl_radionuc') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					incl_radionuc = .FALSE.
				else
					incl_radionuc = .TRUE.
				endif
				if (radionucfile .eq. 'None') incl_radionuc = .FALSE.
				if (radionucfile .eq. 'NONE') incl_radionuc = .FALSE.
				if (radionucfile .eq. 'none') incl_radionuc = .FALSE.

			else if (trim(flagname).EQ.'incl_isrf') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					incl_isrf = .FALSE.
				else
					incl_isrf = .TRUE.
				endif
				if (isrffile .eq. 'None') incl_isrf = .FALSE.
				if (isrffile .eq. 'NONE') incl_isrf = .FALSE.
				if (isrffile .eq. 'none') incl_isrf = .FALSE.
			else if (trim(flagname).EQ.'incl_2dabun') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					incl_2dabun = .FALSE.
					print *, "1D abundances!!!!!!!!!!!!!!!!!!!!!"
				else
					incl_2dabun = .TRUE.
					print *, "2D abundances!!!!!!!!!!!!!!!!!!!!!"
				endif
				if (abun2dfile .eq. 'None') incl_2dabun = .FALSE.
				if (abun2dfile .eq. 'NONE') incl_2dabun = .FALSE.
				if (abun2dfile .eq. 'none') incl_2dabun = .FALSE.
			else if (trim(flagname).EQ.'write_2dabun') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					write_2dabun = .FALSE.
				else
					write_2dabun = .TRUE.
				endif
			else if (trim(flagname).eq.'spatial_dust') then
                read(flagval, '(I3)') intflagval
                if (intflagval .EQ. 0) then
                    incl_locdust = .FALSE.
                else
                    incl_locdust = .TRUE.
                endif

C Testing values
			else if (trim(flagname).EQ.'ratetest') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					ratetest = .FALSE.
				else
					ratetest = .TRUE.
				endif
			else if (trim(flagname).EQ.'testspec' .AND. ratetest) then
				flagval = line(index(line,'=')+1:)
				do i=1,maxratetest
					pos = index(flagval, ',')
					IF (pos .EQ. 0) THEN
						RTspec(i) = adjustl(flagval)
						RTspec(i+1) = 'LASTSPEC'
						Nratetest = i
						EXIT
					ENDIF
					RTspec(i) = adjustl(flagval(:pos-1))
					flagval = flagval(pos+1:)
				enddo

			else if (trim(flagname).EQ.'shieldtest') then
				read(flagval, '(I3)') intflagval
				if (intflagval .EQ. 0) then
					shieldtest = .FALSE.
				else
					shieldtest = .TRUE.
				endif
			endif

			endif
		enddo
	  close(01)
C          print *, "read 5flags"

      if ((ndust .ne. 1.0).and.(freezeeffic .ne. 1.0)) then
        print *, 'Warning, both freezeeffic and epsilon set.'
        print *, 'Are you sure you want to do this?'
      end if


C..............................................................................
C 6) Binding energies for grain surface reactions
C..............................................................................
      open (unit=01,file='6grainbe.inp',status='old')
		fend = 0
		nfrc = 0
		do while (fend .EQ. 0)
			read(01,'(A160)',iostat=fend) line
			if (line(1:1) .NE. '#') then
				nfrc = nfrc + 1
				if (index(line,'#').ne.0) then
					line = line( :index(line,'#')-1)
				endif
				read(line,*) bindspec(nfrc), bindenerg(nfrc)
C				write(*,*) bindspec(nfrc), bindenerg(nfrc)
			endif
		enddo
      close(01)

C..............................................................................
C 1) Environmental conditions:
C..............................................................................
	  call getarg(2,envfile)
      open (unit=01,file=envfile,status='old', access='sequential')
	  read (01,*)
	  read (01,*) Nr
	  read (01,*) Nz

	  if (iend.gt.(Nr*Nz)) then
		write(*,*) 'iend = ',iend,' > Nr*Nz = ',Nr*Nz,', stop'
		stop
	  endif
C Read in each zone's parameters
      DO i = 1, Nr*Nz, 1
        if (.not. incl_locdust) then
           print *, 'assume rho_d = rho_g/100'
           read (01,*) Rs(i), rho(i), Tg(i), Td(i), zAU(i), zcm(i),
     &				Nrz(i), zetaCR(i)
        else
           read (01,*) Rs(i), rho(i), ngr(i), Tg(i), Td(i), zAU(i),
     &                           zcm(i), Nrz(i), zetaCR(i),locdust(i)
              print *, 'Local dust fraction adjusted by: ',i,locdust(i)
        end if
C If Tgas <= 0, assume it's not set and set Tg = Td.
		if (Tg(i) .LE. 0) Tg(i) = Td(i)
c                Tg(i) = Tg(i)+20.0
c                Td(i) = Td(i)+20.0
c KRS
	  END DO

c Create dust string:
      lastzer = -1
	  WRITE(dustn, fmt = '(F8.6)') ndust
      tmpdust = dustn
	  DO i=1,len_trim(dustn)
		IF (dustn(i:i) .EQ. '.') THEN
			tmpdust(i:i) = 'p'
        else if ((dustn(i:i) .EQ. '0').and.(lastzer.eq.-1)) THEN
            lastzer = i
		ELSE IF ((dustn(i:i) .NE. '0').and.(lastzer.ne.-1)) THEN
			lastzer = -1
		END IF
	  ENDDO
      if (lastzer.ne.-1) dust = 'e' // trim(tmpdust(1:lastzer-1))
      if (lastzer.eq.-1) dust = 'e' // trim(tmpdust)
      if ( dust(len_trim(dust):len_trim(dust)) .eq. 'p' ) then
        dust=dust(1:len_trim(dust)-1)
      END IF

c Set output file prefix and dust value
	  endenv = envfile(index(envfile, '1environ')+13:len(envfile))
	  dot = index(endenv, '.')
c	  dust = endenv(1:dot-1)
	  out = 'r' // trim(endenv(dot+1:len(endenv))) // '_' // trim(dust)

	  close(01)


C..............................................................................
C Run static chemistry over a 1D grid:
C..............................................................................
	  firstzone = .TRUE.
c	  open(unit=98, file='testrates.dat', status='new')
      DO i = 1, Nr*Nz, 1
		zone = i
		print *, 'zone:',zone
		print *, 'radius:',Rs(i)

		if (firstzone) then
C Calculate the CR attenuation factor
			CALL calcCRatten()
C Calculate the surface density for the RN escape factor
			CALL calcRNatten()
C Read in UV field from uvfield (needs Rs(1)):
			CALL readuv(uvfile)


		IF (incl_isrf) THEN
            print *, 'Including ISRF, filename ',trim(isrffile)
			CALL readisrf(isrffile)
		ENDIF

C file name setup for various test files
c			preadd = '/Nirgal1/fogel/astrochem/semenov/'
			preadd = ''
			radstring = out(2:index(out, '_')-1)

C Create selfshield.test if needed
			IF (shieldtest) THEN
				fnameSS = trim(preadd) // 'selfshield_' // trim(dust) // '_'
				fnameSS = trim(fnameSS) // trim(radstring) // '.test'

				OPEN(unit=02, file=fnameSS, status='replace')
				WRITE(02,*)
				CLOSE(02)
			END IF

C Create <species>_<dust>_<radius>.rout files for rate test
			IF (ratetest) THEN
				postadd = '_' // trim(dust) // '_' // trim(radstring)
				postadd = trim(postadd) // '.rout'
				text2 = '# shell number  time [yr]  reaction number 1'
				text2 = trim(text2) // ' reaction rate 1 [cm-3/s]...'

				DO j=1,Nratetest
					text1 = '# Main ' // trim(RTspec(j)) // ' formation/'
					text1 = trim(text1) // 'destruction routes'
					text1 = trim(text1) //  ' computed by astrochem'
					IF (trim(RTspec(j)) .EQ. 'LASTSPEC') THEN
						EXIT
					ELSE
						fnameRT(j) = trim(preadd) // trim(RTspec(j))
						fnameRT(j) = trim(fnameRT(j)) // postadd
						OPEN(unit=02, file=fnameRT(j), status='replace')
						WRITE(02, '(a80)') text1
						WRITE(02, '(a80)') text2
						WRITE(02,'(a1)') '#'
						CLOSE(02)
					ENDIF
				ENDDO
			ENDIF

C Read in x-ray scattering file
		CALL readxray(xrayfile)

C Read in X-ray Pabs Cross-X data.

			open (unit=59,file='xs_format.data',status='old')
			nln = 1
			do while (nln .lt. 2001)
				read(59,130) e_sil(nln), xs_sil(nln)
				nln = nln + 1
			enddo
			close(59)

  130    format(1x,e12.6,2x,e12.6)

C End of first zone setup stuff
		ENDIF

        write(*,*) 'Grid point #', i

		IF (incl_radionuc) THEN
!		    radionucfile = 'SLR_combined.dat'
			CALL readRN(radionucfile)
			print *, 'Read in radionuclide file: ',radionucfile
		ENDIF

C Convert zone number to string for input/output
        ipon = i
        IF (ipon.LE.9) THEN
          WRITE (ait,'(I1)') ipon
          nait = 1
        ELSE IF (ipon.LE.99) THEN
          WRITE (ait,'(I2)') ipon
          nait = 2
        ELSE  IF (ipon.LE.999) THEN
          WRITE (ait,'(I3)') ipon
          nait = 3
        ELSE  IF (ipon.LE.9999) THEN
          WRITE (ait,'(I4)') ipon
          nait = 4
        ELSE
          WRITE (ait,'(I5)') ipon
          nait = 5
        END IF

C Since Bethell's UV field is not acurate outside of a 400 x 400 AU box, set
C abundances to 0 if we are outside of that box and skip running the
C ODE solver there
		IF (zAU(zone) .GT. ZMAX .OR. Rs(zone) .GT. RMAX .or. Tg(zone).ge.1000.0) THEN
			DO j=1,nspec
				yy(j) = 0.0
			ENDDO
C Ignore cells with Temperature = 0 or Density = 0
        ELSE IF (Tg(zone).EQ.0 .OR. Td(zone).EQ.0 .OR. rho(zone).LE. 1E-24) THEN
			DO j=1,nspec
				yy(j) = 0.0
			ENDDO
C Run zone
		ELSE
C Total gas particle density:
			gdens = rho(zone) / aMp / amu

C Set initial abundances
C If using 2D initial abundances read them in now
            if (incl_2dabun .eqv. .TRUE.) then
            	file2d = trim(out)//'_'//ait(1:nait)//'_'
     &               //trim(adjustl(abun2dfile))//'.inp'
			newfile2d = trim(adjustl(file2d))
			print *,file2d
			print *,newfile2d
			print *,abun2dfile
                open(unit=71,file=file2d,
     &				 status='unknown')
	  			nfrc=0
	  			fend = 0
	  			do while (fend .eq. 0)
	  				read(71,'(A160)',iostat=fend) line
	  				if (line(1:1) .NE. '#' .and. fend .eq. 0) then
	  					read(line,*) spec0(nfrc), frc(nfrc)
	  					print *,spec0(nfrc),frc(nfrc),'!!!!!!!!!!!!!!!!!!!!!!!!!!'
	  					nfrc = nfrc+1
	  				endif
	  			enddo
	  			close(71)

                print *,'nfrc'
                CALL ini_abunds2D(ns,s,nfrc,spec0,frc,gdens,yy)
            else
                CALL ini_abunds(ns,s,nfrc,spec0,frc,gdens,yy)
            endif

C Calculate the x-ray ionization rates
            xrayrate(zone) = xrayion_rate()

            grnfrac = 6e-12

C Calculate the rate coefficients:
			timedep = 0
			DO j = 1, nre
			  rr2 = r2(j)
			  ak(j) = calcrate(alpha(j), beta(j), gamma(j), rtype(j), r1(j),
c     &		r2(j),frc(12),timedep,dtaudtemp)
     &		r2(j),grnfrac,timedep,dtaudtemp,tfirst)
			END DO

C Compute chemistry:
			CALL run_chemistry(ns,tfirst,nstep,yy,tlast,relerr,abserr)

C IF Statement for max z check
		ENDIF

C Write results to the output file:

C Write abundances at the last time step to an additional file
                if (write_2dabun .eqv. .TRUE.) then
                   l = nstep
C                   write(*,*) "writing abundance file",tlast/year
                   lastyr = tlast/year
                   WRITE (lts,'(1PE10.1)') lastyr
                   write(*,*) lts
                   OPEN (unit=21,file=trim(out)//'_'//ait(1:nait)//'_'
     &               //trim(adjustl(lts))//'.inp',status='unknown',
     &               access='append')
                   REWIND 21
                   do j = 1, ns
                      write(21,66) s(j),abundances(j,l)/gdens
                   enddo
                   CLOSE(21)
                   endif
 66             format(1a13,1pe11.3)


C Open output file 'out':
        OPEN (unit=20,file=trim(out)//'_'//ait(1:nait)//'.out',
     &        status='unknown',access='append')
        REWIND 20

      li = ns
      lt = li+1
      iana =1
      write(20,20)
      write(20,75)
 75   format(21x,39('*'))
      write(20,84) lt-1
 84   format(20X,' **       DISK CHEMISTRY             **',/,
     * 20X,' **          ',1I3,' SPECIES SET          **',/,
     * 20X,' **            F77 VERSION            **',/,
     * 20X,' **             02/17/2015            **')
      write(20,75)
      write(20,20)
      write(20,22) (s(j),j=1,ns)
 22   format(8(1x,a13))
      write(20,23)
 23   format(/)
      write(20,36) ns
 36   format(1x,' NUMBER OF VALID SPECIES USED = ',1i4)
      write(20,45) nre
 45   format(1x,' NUMBER OF VALID REACTIONS USED = ',1i6)
      write(20,23)
 20   format(//)

C Write input parameters:
      write(20,64)
 64   format(3x,' INITIAL VALUES: '/)
      WRITE(20,65) ipon, Rs(zone), zAU(zone), gdens, Tg(zone), Td(zone),
     1     rho(zone), ZetaCR(zone), albedo_UV,
     2     times(1)/year, tlast/year
 65   FORMAT(
     * 3X,' Grid point  = ',I5,/,
     * 3X,' Radius      = ',1PE11.3,' AU',/,
     * 3X,' Height (z)  = ',1PE11.3,' AU',/,
     * 3X,' n(H+2H2)    = ',1PE11.3,' cm^(-3)',/,
     * 3X,' Tgas        = ',0PF8.1,'K',/,
     * 3X,' Tdust       = ',0PF8.1,'K',/,
     * 3X,' rho_g       = ',1PE11.3,' g/cm^3',/,
     * 3X,' ZetaCRP     = ',1PE11.3,' s^-1',/,
     * 3X,' albedo(UV)  = ',1PE11.3,/,
     * 3X,' Start time  = ',1PE10.1,' years',/,
     * 3X,' Finish time = ',1PE10.1,' years')
      WRITE(20,20)

C Write calculated abundances:
      write(20,63)
63    format(3X,' CALCULATED ABUNDANCES: '/)
       is = 1
       iff = 6
32    write(20,41) (s(j),j=is,iff)
      write(20,76)
      do l = 1, nstep
         write(20,30) times(l)/year,(abundances(j,l)/gdens,j=is,iff)
      end do
30    format(1x,7(1pe11.3))
      write(20,41)(s(j),j=is,iff)
41    format(6x,'time',6x,6(1a13,1x))
      write(20,20)
      is=is+6
      iff=iff+6


C Interupt?
C      print *,s(iff),s(is),s(ns),li,is,ns
      if (iff.gt.ns) iff=ns
      if (is.gt.li) go to 38
      go to 32

C Last output statement:
38    lx=li
      lt=lt-1
76    format(1x,80('-'))

      CLOSE (20)

	    IF (firstzone) firstzone = .FALSE.

      END DO
c		close(98)

C Close all file(s):
      close (01)

C Exit:
      END

C..............................................................................
C
C This subroutine initializes  abundances of chemical species for t=0.
C
C..............................................................................
C
C Input parameter(s):
C
C ny			 == number of species
C
C y(ny)          == the names of non-conserved species,
C
C n				 == number of species with initial abundances
C
C yr(n)          == the names of chemical species from '3abunds.inp',
C
C fraction(1:n)  == fractions of the corresponding species 'yr',
C
C density        == total amount of hydrogen nuclei
C
C..............................................................................
C
C Output parameter(s):
C
C y0(ny)         == computed initial abundances
C
C..............................................................................
      SUBROUTINE ini_abunds(ny,y,n,yr,fraction,density,y0)
      IMPLICIT NONE

	  INCLUDE "rates.h"
	  INCLUDE "environ.h"

C Global variable(s):
      INTEGER ny, n, namelen
      CHARACTER*13 y, yr
      REAL*8 fraction, density, y0, eabun
      DIMENSION y(ny), y0(ny), yr(n), fraction(n)

C Local variable(s):
      INTEGER i, j, k
      DIMENSION k(n)

C Initialization:
        do i = 1, ny
          y0(i) = 0.0D0
        end do
        do i = 1, n
          k(i) = 0
        end do

C Search of species 'yr' among species 'y':
      i = 0

C A loop by 'yr':
10    continue
        i = i + 1
        if (i.gt.n)  goto 40  ! checked all 'yr', exit,

C A loop by 'y':
        j = 0
20      continue
        j = j + 1
        if (j.gt.ny) goto 10

C Found a necessary species:
        if (yr(i).eq.y(j)) then
          k(i) = j
          goto 10  ! switch to the next species 'yr',
        end if

        goto 20  ! switch to the next species 'y',

40    continue

C Calculate the initial E abundance from the initial ion abundances
	eabun = 0
	do i = 1,n
		namelen = len_trim(yr(i))
		if (yr(i)(namelen:namelen) .EQ. '+') then
			eabun = eabun + fraction(i)
		end if
	end do

	do i = 1,ny
		if (y(i) .EQ. 'E') then
			y0(i) = eabun * density
			print *, 'E abundance calculated: ', y0(i)
		end if
	end do

C Calculate the initial abundance of species on grains and the initial
C grain abundance
	n_ice_init = 0
	do i = 1,n
		IF ('(gr)' .EQ. yr(i)(len_trim(yr(i))-3:)) THEN
			n_ice_init = n_ice_init + fraction(i)*density
		ENDIF
C scale grain abundance by dust settling parameter (except for midplane zone)
C Deprecated, KRS 5/20/22, now read in from 1environ
C		IF ('GRAIN' .EQ. trim(yr(i))) THEN
C			IF (zone .NE. Nz) THEN
C				ngr_init = fraction(i)*density*ndust
C			ELSE
C				ngr_init = fraction(i)*density
C			ENDIF
C		ENDIF
	end do

C Calculation of the initial abundances:
      do i = 1, n
		if (k(i).NE.0) then
			print *,'Initially present is ',y(k(i))
			if (y(k(i)) .EQ. 'GRAIN' .OR. y(k(i)) .EQ. 'GRAIN0') then
				IF (zone .NE. Nz) THEN
					y0(k(i)) = fraction(i)*density*ndust
				ELSE
					y0(k(i)) = fraction(i)*density
				ENDIF
			else
				y0(k(i)) = fraction(i)*density
        C KRS need to change this
			endif
		  endif
      end do

C Exit:
      return
      end

C..............................................................................
C
C This subroutine initializes  abundances of chemical species for t!=0.
C
C..............................................................................
C
C Input parameter(s):
C
C ny			 == number of species
C
C y(ny)          == the names of non-conserved species,
C
C n				 == number of species with initial abundances
C
C yr(n)          == the names of chemical species from 'r*.inp',
C
C fraction(1:n)  == fractions of the corresponding species 'yr',
C
C density        == total amount of hydrogen nuclei
C
C..............................................................................
C
C Output parameter(s):
C
C y0(ny)         == computed initial abundances
C
C..............................................................................
      SUBROUTINE ini_abunds2D(ny,y,n,yr,fraction,density,y0)
      IMPLICIT NONE

	  INCLUDE "rates.h"
	  INCLUDE "environ.h"

C Global variable(s):
      INTEGER ny, n, namelen
      CHARACTER*13 y, yr
      REAL*8 fraction, density, y0, eabun
      DIMENSION y(ny), y0(ny), yr(n), fraction(n)

C Local variable(s):
      INTEGER i, j, k
      DIMENSION k(n)
      print *,"in ini_abunds2D"

C Initialization:
        do i = 1, ny
          y0(i) = 0.0D0
        end do
        do i = 1, n
          k(i) = 0
        end do
        print *,"line 1053"

C Search of species 'yr' among species 'y':
      i = 0

C A loop by 'yr':
10    continue
        i = i + 1
        if (i.gt.n)  goto 40  ! checked all 'yr', exit,

C A loop by 'y':
        j = 0
20      continue
        j = j + 1
        if (j.gt.ny) goto 10

C Found a necessary species:
        if (yr(i).eq.y(j)) then
          k(i) = j
          goto 10  ! switch to the next species 'yr',
        end if

        goto 20  ! switch to the next species 'y',

40    continue


C Calculate the initial abundance of species on grains and the initial
C grain abundance
	n_ice_init = 0
	print *,"line 1083"
	do i = 1,n
		IF ('(gr)' .EQ. yr(i)(len_trim(yr(i))-3:)) THEN
			n_ice_init = n_ice_init + fraction(i)*density
		ENDIF
C scale grain abundance by dust settling parameter (except for midplane zone)
		IF ('GRAIN' .EQ. trim(yr(i))) THEN
			IF (zone .NE. Nz) THEN
				ngr_init = fraction(i)*density*ndust
			ELSE
				ngr_init = fraction(i)*density
			ENDIF
		ENDIF
	end do

C Calculation of the initial abundances:
	  print *,"line 1099"
	  print *,n
      do i = 1, n
      	print *,i,k(1)
		if (k(i).NE.0) then
		    print *,'Initially present is ',y(k(i))
			if (y(k(i)) .EQ. 'GRAIN' .OR. y(k(i)) .EQ. 'GRAIN0') then
				IF (zone .NE. Nz) THEN
					y0(k(i)) = fraction(i)*density*ndust
				ELSE
					y0(k(i)) = fraction(i)*density
				ENDIF
			else
				y0(k(i)) = fraction(i)*density
			endif
			print *,y0(k(i))
		  endif
      end do

C Exit:
      return
      end

C..............................................................................
C
C This subroutine returns index of a given species in a given array of species
C names.
C
C..............................................................................
C
C Input parameter(s):
C
C species      == name of a given species,
C
C na		   == number of species
C
C array(na)    == a given array of species names,
C
C..............................................................................
C
C Output parameter(s):
C
C index      == index of 'species' among 'array(1:na)'
C
C..............................................................................
      SUBROUTINE ispecies(species,na,array,index)
      IMPLICIT NONE

C Global variable(s):
      INTEGER na, index
      CHARACTER*13 species, array
      DIMENSION array(na)

C Local variable(s):
      INTEGER i

C Start search of 'species' in 'array':
      i = 1
 10   CONTINUE
      IF (i.gt.na) GOTO 20      ! No 'species' in 'array', stop,
      IF (species.eq.array(i)) GOTO 30 ! Found 'species', exit,
         i = i + 1
         GOTO 10
C 'species' was not found in 'array', stop:
 20   CONTINUE
cd      write (56,*) species
      index = 0
      RETURN
C 'Species' was found in 'array', initialize 'index':
 30   CONTINUE
      index = i
C Exit:
      RETURN
      END
