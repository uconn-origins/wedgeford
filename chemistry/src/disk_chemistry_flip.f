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
C                  from Nuria Calvet's codes.
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
C g0_100      == G0 value, normalized at 100 AU
C
C g0_CO          == G0 value for CO (calculated ignoring Lyman alpha radiation)
C
C Lxray          == X-ray luminosity [ergs/s]
C
C RadNuc      == Radionuclide Ionization Rate [1/s]
C
C fg          == dust to gas mass ratio (normal is 100)
C
C freezeeffic == parameter to approximate freeze out effects in grain growth. Normal (0.1 mu grains) is 1.
C
C CRdesorp      == turn on/off cosmic ray desorption
C
C CRionization == turn on/off cosmic ray ionization
C
C photodesorp == turn on/off UV photodesorption
C
C LyAphotodesorp    == turn on/off lyman alpha photodesorption
C
C thermaldesorp        == turn on/off thermal desorption
C
C include_LyA == include lyman alpha radiation or remove it from UV field
C
C xraydust        == turn on/off dust-dependent x-ray opacities
C
C incl_radionuc == turn on/off radionuclide ionization.  Must include external file w/ RN calculations.
C
C ratetest      == print out rates of given species or not
C
C RTspec(maxratetest)    == species to print rates of
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
C zAU(Nz)      == Height above the midplane [AU]
C
C zcm(Nz)      == Height from disk surface, surface = 1.0e10 cm [cm]
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
C        ini_abunds2D, readr, reads, readuv, run_chemistry
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
     1  yy, gdens, t, tout, tstep, abunzone, akzone
      CHARACTER*13 spec0, rr2
      CHARACTER*80 specs, rates, out, ver, cmgtm, uvfile, xrayfile
      CHARACTER*80 radionucfile
      CHARACTER*80 isrffile
      CHARACTER*80 xrayfilesb
      CHARACTER*80 dustn
      CHARACTER*80 abun2dfile
      CHARACTER*80 outabunfile
      CHARACTER*28 file2d,newfile2d
      DIMENSION frc(nspec2), spec0(nspec2), yy(nspec), abunzone(nspec,nzone), akzone(nreac,nzone)
      real grnfrac, start_time, finish_time

C External Functions
      DOUBLE PRECISION xrayion_rate, calcrate, dtaudtemp, polint

C Local variable(s):
      DOUBLE PRECISION errval, lastyr
      INTEGER i, j, k, l, nlen, nait, ipon, iana,
     1  iff, lx, is, li, lt, timedep, pos
      CHARACTER*5 ait, siend, eqlsign
      CHARACTER*160 line,testfi
      CHARACTER*13 errspec, lts
      CHARACTER*20 flagname
      CHARACTER*140 flagval
      CHARACTER*100 iofile, envfile
      CHARACTER*80 preadd, radstring, postadd, text1, text2, endenv
      INTEGER intflagval, dot, errspecindex
c      real dtaudt(9000,3)
      integer nln,indv
      INTEGER jg ! iterates final zone species for lower error, temporary
      LOGICAL sameness, old_format


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
      CALL getarg(1,iofile)
      OPEN (unit=01,file=iofile,status='old',access='sequential')
      READ (01,*)
      READ (01,'(A50)') specs
        print *, specs
      READ (01,'(A50)') rates
        print *, rates
      READ (01,'(A90)') uvfile
        uvfile = uvfile( :index(uvfile,'#')-1)
        print *, uvfile
      READ (01,'(A90)') xrayfile
        xrayfile = xrayfile( :index(xrayfile,'#')-1)
        print *, xrayfile
      READ (01,'(A90)') isrffile
        isrffile = isrffile( :index(isrffile,'#')-1)
C             print *, isrffile
      READ (01,'(A90)') radionucfile
        radionucfile = radionucfile( :index(radionucfile,'#')-1)
        print *, radionucfile
      READ (01,'(A90)') abun2dfile
        abun2dfile = abun2dfile( :index(abun2dfile,'#')-1)
C             print *, abun2dfile
      CLOSE (01)
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
      OPEN (unit=01,file='2times.inp',status='old',access='sequential')
       READ (01,*)
       READ (01,*) tlast
       tlast = tlast * year
       READ (01,*) tfirst
       tfirst = tfirst * year
       READ (01,*) nstep
      CLOSE (01)
       print *, "read 2times"

C..............................................................................
C 3) Initial abundances:
C..............................................................................
      IF (tfirst/year .LT. 9.000D-00) THEN
       OPEN (unit=01,file='3abunds.inp',status='old')
        fend = 0
        nfrc = 0
        DO WHILE (fend .EQ. 0)
          READ(01,'(A160)',iostat=fend) line
          IF (line(1:1) .NE. '#') THEN
             nfrc = nfrc + 1
             READ(line,*) spec0(nfrc), frc(nfrc)
         END IF
         END DO

       CLOSE(01)
      END IF
      print *, "read 3abunds"

C..............................................................................
C 4) Tolerance:
C..............................................................................
      OPEN (unit=01,file='4toleran.inp',status='old')
      READ (01,*)
      READ (01,*) errval
      DO i=1,nspec
        relerr(i) = errval
      END DO
      READ (01,*) errval
      DO i=1,nspec
        abserr(i) = errval
      END DO
c        set species-specific tolerances
      fend = 0
      DO WHILE (fend .EQ. 0)
        READ(01,'(A160)',iostat=fend) line
        IF (line(1:1) .EQ. '#') CYCLE        
          READ(line,*) flagname, errspec, errval
          CALL ispecies(errspec, ns, s, errspecindex)
          IF (trim(flagname) .EQ. 'rel') THEN
            relerr(errspecindex) = errval   
          ELSE IF (trim(flagname) .EQ. 'abs') THEN
            abserr(errspecindex) = errval
          ELSE
            print *, 'Illegal value in 4toleran.inp'
            STOP
        END IF
      END DO
      CLOSE (01)
      print *, "read 4tol"

      !  Initialize freezeeffic (read in below)
      freezeeffic = 1.0

C..............................................................................
C 5) Flags:
C..............................................................................
      OPEN (unit=01,file='5flags.inp',status='old')
      fend = 0
      DO WHILE (fend .EQ. 0)
        READ(01,'(A160)',iostat=fend) line
        IF (line(1:1) .NE. '#') then
          READ(line,*) flagname, eqlsign, flagval
c 100            format(A13,3X,A140)

C Parameters
        IF (trim(flagname).EQ.'fg') then
          READ(flagval, '(E8.2)') fg
          DO i=1,nfrc
            IF (trim(spec0(i)).NE.'H2'.AND.trim(spec0(i)).NE.'He') THEN
              frc(i) = frc(i) * 100.0 / fg
            END IF
          END DO
 
        ELSE IF (trim(flagname).EQ.'freezeeffic') THEN
          READ(flagval, '(E8.2)') freezeeffic
        ELSE IF (trim(flagname).EQ.'epsilon') THEN
          READ(flagval, '(E8.2)') ndust

C Flags
        ELSE IF (trim(flagname).EQ.'CRdesorp') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) then
              CRdesorp = .FALSE.
            ELSE
              CRdesorp = .TRUE.
            END IF
        ELSE IF (trim(flagname).EQ.'CRionization') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) then
              CRionization = .FALSE.
            ELSE
              CRionization = .TRUE.
            END IF
        ELSE IF (trim(flagname).EQ.'photodesorp') THEN
          READ (flagval, '(I3)') intflagval
            IF(intflagval .EQ. 0) then
              photodesorp = .FALSE.
            ELSE
              photodesorp = .TRUE.
            END IF
        ELSE IF (trim(flagname).EQ.'LyAphotodesorp') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
              LyAphotodesorp = .FALSE.
            ELSE
              LyAphotodesorp = .TRUE.
            END IF
        ELSE IF (trim(flagname).EQ.'thermaldesorp') THEN
            READ(flagval, '(I3)') intflagval
              IF (intflagval .EQ. 0) THEN
                thermaldesorp = .FALSE.
              ELSE
                thermaldesorp = .TRUE.
              END IF
        ELSE IF (trim(flagname).EQ.'include_lya') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
              include_LyA = .FALSE.
            ELSE
              include_LyA = .TRUE.
            END IF
        ELSE IF (trim(flagname).EQ.'xraydust') THEN
          READ(flagval, '(I3)') intflagval
          IF (intflagval .EQ. 0) THEN
              xraydust = .FALSE.
          ELSE
              xraydust = .TRUE.
          END IF
        ELSE IF (trim(flagname).EQ.'incl_radionuc') then
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) then
              incl_radionuc = .FALSE.
            ELSE
              incl_radionuc = .TRUE.
            END IF
            IF (radionucfile .eq. 'None') incl_radionuc = .FALSE.
            IF (radionucfile .eq. 'NONE') incl_radionuc = .FALSE.
            IF (radionucfile .eq. 'none') incl_radionuc = .FALSE.

        ELSE IF (trim(flagname).EQ.'incl_isrf') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
              incl_isrf = .FALSE.
            ELSE
              incl_isrf = .TRUE.
            END IF
            IF (isrffile .eq. 'None') incl_isrf = .FALSE.
            IF (isrffile .eq. 'NONE') incl_isrf = .FALSE.
            IF (isrffile .eq. 'none') incl_isrf = .FALSE.
        ELSE IF (trim(flagname).EQ.'incl_2dabun') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
              incl_2dabun = .FALSE.
              print *, "1D abundances!"
            ELSE
              incl_2dabun = .TRUE.
              print *, "2D abundances!"
            END IF
            IF (abun2dfile .eq. 'None') incl_2dabun = .FALSE.
            IF (abun2dfile .eq. 'NONE') incl_2dabun = .FALSE.
            IF (abun2dfile .eq. 'none') incl_2dabun = .FALSE.
        ELSE IF (trim(flagname).EQ.'write_2dabun') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
              write_2dabun = .FALSE.
            ELSE
             write_2dabun = .TRUE.
            END IF
        ELSE IF (trim(flagname).eq.'spatial_dust') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
             incl_locdust = .FALSE.
            ELSE
              incl_locdust = .TRUE.
            END IF
        ELSE IF (trim(flagname).eq.'old_format') THEN
           READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
              old_format = .FALSE.
            ELSE
               old_format = .TRUE.
            END IF

C Testing values
        ELSE IF (trim(flagname).EQ.'ratetest') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
              ratetest = .FALSE.
            ELSE
              ratetest = .TRUE.
            END IF
        ELSE IF (trim(flagname).EQ.'testspec' .AND. ratetest) THEN
          flagval = line(index(line,'=')+1:)
          DO i=1,maxratetest
            pos = index(flagval, ',')
            IF (pos .EQ. 0) THEN
              RTspec(i) = adjustl(flagval)
              RTspec(i+1) = 'LASTSPEC'
              Nratetest = i
              EXIT
            END IF
             RTspec(i) = adjustl(flagval(:pos-1))
             flagval = flagval(pos+1:)
          END DO !end rate test

        ELSE IF (trim(flagname).EQ.'shieldtest') THEN
          READ(flagval, '(I3)') intflagval
            IF (intflagval .EQ. 0) THEN
              shieldtest = .FALSE.
            ELSE
              shieldtest = .TRUE.
            END IF
        END IF !flag tests

        END IF !end of check for appropriate lines
      END DO !of reading lines
      CLOSE(01)
C         print *, "read 5flags"

      IF ((ndust .ne. 1.0).and.(freezeeffic .ne. 1.0)) THEN
        print *, 'Warning, both freezeeffic and epsilon set.'
        print *, 'Are you sure you want to do this?'
      END IF


C..............................................................................
C 6) Binding energies for grain surface reactions
C..............................................................................
      OPEN (unit=01,file='6grainbe.inp',status='old')
        fend = 0
        nfrc = 0
        DO WHILE (fend .EQ. 0)
          READ(01,'(A160)',iostat=fend) line
          IF (line(1:1) .NE. '#') THEN
            nfrc = nfrc + 1
            IF (index(line,'#').ne.0) THEN
              line = line( :index(line,'#')-1)
            END IF
            READ(line,*) bindspec(nfrc), bindenerg(nfrc)
          END IF
        END DO
      CLOSE(01)

C..............................................................................
C 1) Environmental conditions:
C..............................................................................
      CALL getarg(2,envfile)
      OPEN (unit=01,file=envfile,status='old', access='sequential')
      READ (01,*)
      READ (01,*) Nr
      READ (01,*) Nz
 
      IF (iend.gt.(Nr*Nz)) THEN
        WRITE(*,*) 'iend = ',iend,' > Nr*Nz = ',Nr*Nz,', stop'
        STOP
      END IF
C Read in each zone's parameters
      IF (old_format) then
        DO i = 1, Nr*Nz, 1
          IF (.not. incl_locdust) THEN
            print *, 'assume rho_d = rho_g/100'
            READ (01,*) Rs(i),rho(i),Tg(i),Td(i),zAU(i),zcm(i),Nrz(i),zetaCR(i)
            ngr(i) = rho(i)*1e-2*0.75/(1.4*pi*1e-21)
          ELSE
            READ (01,*) Rs(i),rho(i),Tg(i),Td(i),zAU(i),zcm(i),Nrz(i),zetaCR(i),locdust(i)
            print *, 'Local dust fraction adjusted by: ',i,locdust(i)
            ngr(i)=rho(i)*1e19*0.75/(1.4*pi*locdust(i)*sqrt(locdust(i)))
          END IF
          IF (Tg(i) .LE. 0) Tg(i) = Td(i) !iset Tg=Td for neg vals
        END DO !old format read loop
c internal density in ngr calc from Krijt & Ciesla 16
      ELSE
        DO i = 1, Nr*Nz, 1
          IF (.not. incl_locdust) THEN
            READ (01,*) Rs(i),rho(i),ngr(i),Tg(i),Td(i),zAU(i),zcm(i),zetaCR(i)
          ELSE
            READ (01,*) Rs(i),rho(i),ngr(i),Tg(i),Td(i),zAU(i),zcm(i),zetaCR(i),locdust(i)
            print *, 'Local dust fraction adjusted by: ',i,locdust(i)
          END IF
C If Tgas <= 0, assume it's not set and set Tg = Td.
          IF (Tg(i) .LE. 0) Tg(i) = Td(i)
        END DO ! new format read loop
      END IF ! all read loops

c Create dust string:
      lastzer = -1
      WRITE(dustn, fmt = '(F8.6)') ndust
      tmpdust = dustn
      DO i=1,len_trim(dustn)
        IF (dustn(i:i) .EQ. '.') THEN
          tmpdust(i:i) = 'p'
        ELSE IF ((dustn(i:i) .EQ. '0').and.(lastzer.eq.-1)) THEN
          lastzer = i
        ELSE IF ((dustn(i:i) .NE. '0').and.(lastzer.ne.-1)) THEN
          lastzer = -1
        END IF
      END DO
      IF (lastzer.ne.-1) dust = 'e' // trim(tmpdust(1:lastzer-1))
      IF (lastzer.eq.-1) dust = 'e' // trim(tmpdust)
      IF ( dust(len_trim(dust):len_trim(dust)) .eq. 'p' ) THEN
        dust=dust(1:len_trim(dust)-1)
      END IF

c Set output file prefix and dust value
      endenv = envfile(index(envfile, '1environ')+13:len(envfile))
      dot = index(endenv, '.')
c      dust = endenv(1:dot-1)
      out = 'r' // trim(endenv(dot+1:len(endenv))) // '_' // trim(dust)

      CLOSE(01)


C..............................................................................
C Setup the initial conditions to run static chemistry over a 1D grid:
C..............................................................................
      firstzone = .TRUE. !computes the initial conditions for the zone
      DO i = 1, Nr*Nz
        zone = i
        print *, 'zone:',zone
        print *, 'radius:',Rs(i)

        IF (firstzone) THEN
C Calculate the CR attenuation factor
          CALL calcCRatten()
C Calculate the surface density for the RN escape factor
          CALL calcRNatten()
C Read in UV field from uvfield (needs Rs(1)):
          CALL readuv(uvfile)

          IF (incl_radionuc) THEN
            CALL readRN(radionucfile)
            print *, 'Read in radionuclide file: ',radionucfile
          END IF !radionucleide

          IF (incl_isrf) THEN
            print *, 'Including ISRF, filename ',trim(isrffile)
            CALL readisrf(isrffile)
          END IF !isrf read

C file name setup for various test files
          preadd = ''
          radstring = out(2:index(out, '_')-1)

C Create selfshield.test if needed
          IF (shieldtest) THEN
            fnameSS = trim(preadd) // 'selfshield_' // trim(dust) // '_'
            fnameSS = trim(fnameSS) // trim(radstring) // '.test'

            OPEN(unit=02, file=fnameSS, status='replace')
             WRITE(02,*)
            CLOSE(02)
          END IF !shieldtest

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
              END IF
            END DO
          END IF !rate test file creation

C Read in x-ray scattering file
          CALL readxray(xrayfile)

C Read in X-ray Pabs Cross-X data.
          OPEN (unit=59,file='xs_format.data',status='old')
            nln = 1
          DO WHILE (nln .lt. 2001)
            READ(59,130) e_sil(nln), xs_sil(nln)
            nln = nln + 1
          END DO
          CLOSE(59)

130      format(1x,e12.6,2x,e12.6)

        END IF !setup for the very first zone only
        
C Ignore cells where grid won't run
        IF (Tg(zone).EQ.0 .OR. Td(zone).EQ.0 .OR. rho(zone).LE. 1E-25 .OR. Tg(zone) .ge. 1000.0) THEN
          runzone(zone) = .FALSE. !this zone will not run
          DO j=1,nspec
            abunzone(j,zone) = 0.0
          END DO
C Run zone
        ELSE
          runzone(zone) = .TRUE. !this zone will be run
          gdens = rho(zone) / aMp / amu
C Set initial abundances
          CALL ini_abunds(ns,s,nfrc,spec0,frc,gdens,yy)
          abunzone(:,zone) = yy
          xrayrate(zone) = xrayion_rate()
          DO j = 1, nre
            timedep = 0
            grnfrac = 6e-12
            rr2 = r2(j)
            akzone(j,zone) = calcrate(alpha(j), beta(j), gamma(j), rtype(j), r1(j), r2(j),grnfrac,timedep,dtaudtemp,tfirst)
          END DO
        END IF !determination of if the zone will be run
        IF (firstzone) firstzone = .FALSE.
      END DO !initial zone stuff
      
C..............................................................................
C  Main run of chemistry loop
C..............................................................................
      t = 0.0d0 ! current time
      OPEN (unit = 20, file = trim(out)//'_2d.out', status = 'unknown', access = 'append')
      WRITe(20,*) '# ns =', ns
      WRITE(20,*) '# time = ', t/year
      WRITE(20,*) '# z au = ', (zAU(i), i = 1,Nz)
      CLOSE(20)
      DO k = 1, ns
        CALL write_out_speczone(t,s(k),abunzone(k,:),out)
      END DO
      tout = tfirst ! next time
      tstep = 1.0D+01**(dlog10(tlast/tfirst) / (nstep-1.0d0)) !multiplicative factor for next timestep
      print *, 'initial setup done: starting main loop at t=0'
      CALL cpu_time(start_time)
      DO j = 1, nstep
        timestep = j
        times(j) = tout
        DO i = 1, Nr*Nz
          zone = i
          IF (runzone(zone) .eqv. .TRUE.) THEN
            gdens = rho(zone) / aMp / amu
            DO k = 1, nre
              ak(k) = akzone(k,zone)
            END DO
c            DO k = 1, ns
c              abszone(k) = MAX(abserr(k)*gdens, MINABUN)
c            END DO
            CALL run_chemistry_dt(ns, yy, t, tout, relerr, abserr,abunzone)
          END IF ! end of running the zone
        END DO ! all zones done 
        t = tout
        IF (mod(j-1,nstep/10) .EQ. 0 .OR. j .EQ. nstep) THEN
            OPEN (unit = 20, file = trim(out)//'_2d.out', status = 'unknown', access = 'append')
            WRITE(20,*) '# time = ', t/year
            WRITE(20,*) '# z au = ', (zAU(i), i = 1,Nz)
            CLOSE(20)
            DO k = 1, ns
              CALL write_out_speczone(t,s(k),abunzone(k,:),out) !writeoutput
            END DO
        END IF
        print *, 'timestep =',j,'all zones done!--------'
        CALL cpu_time(finish_time)
        print '("elapsed time = ",f6.3," minutes.")',(finish_time-start_time)/60.0
        tout = tout*tstep ! update new timestep
      END DO


C Close all file(s):
      CLOSE (01)

C Exit:
      END

C..............................................................................
C
C This subroutine writes the abundances by zone to file
C
C..............................................................................
C
C Input parameter(s):
C
C ns             == number of species
C
C nz             == number of zones
C
C s              == species names
C
C rho             == gas density for the zone
C
C abundance      == array of zone by zone species abundances
C
C outname        == name for the file
C
C..............................................................................
C
C..............................................................................
      SUBROUTINE write_out_speczone(t,s,abundk,outname)
      IMPLICIT NONE

      INCLUDE "rates.h"
      INCLUDE "environ.h"
      INCLUDE "constants.h"

C Global variable(s):
      INTEGER ns
      CHARACTER*13 outname, s
      REAL*8  abundk, t
      DIMENSION abundk(nz)

C Local variable(s):
      INTEGER i, j, k
      
C write output file
      OPEN (unit = 20, file = trim(outname)//'_2d.out', status = 'unknown', access = 'append')      
      WRITE(20,66) s, (abundk(j)/ (rho(j) / (aMp*amu) ) ,j=1,Nz)
66    FORMAT(1a13,50(1pe11.3))
      CLOSE(20)
      RETURN
      END



C..............................................................................
C
C This subroutine initializes  abundances of chemical species for t=0.
C
C..............................................................................
C
C Input parameter(s):
C
C ny             == number of species
C
C y(ny)          == the names of non-conserved species,
C
C n                 == number of species with initial abundances
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
      DO i = 1, ny
        y0(i) = 0.0D0
      END DO
      DO i = 1, n
        k(i) = 0
      END DO

C Search of species 'yr' among species 'y':
      i = 0

C A loop by 'yr':
10    CONTINUE
       i = i + 1
      IF (i.gt.n)  GOTO 40  ! checked all 'yr', exit,

C A loop by 'y':
       j = 0
20     CONTINUE
        j = j + 1
        IF (j.gt.ny) GOTO 10

C Found a necessary species:
        IF (yr(i).eq.y(j)) THEN
          k(i) = j
          GOTO 10  ! switch to the next species 'yr',
        END IF

        GOTO 20  ! switch to the next species 'y',

40      CONTINUE

C Calculate the initial E abundance from the initial ion abundances
        eabun = 0
        DO i = 1,n
          namelen = len_trim(yr(i))
          IF (yr(i)(namelen:namelen) .EQ. '+') THEN
            eabun = eabun + fraction(i)
          END IF
        END DO

        DO i = 1,ny
          IF (y(i) .EQ. 'E') THEN
            y0(i) = eabun * density
            print *, 'E abundance calculated: ', y0(i)
          END IF
        END DO

C Calculate the initial abundance of species on grains and the initial
C grain abundance
      n_ice_init = 0
      DO i = 1,n
        IF ('(gr)' .EQ. yr(i)(len_trim(yr(i))-3:)) THEN
          n_ice_init = n_ice_init + fraction(i)*density
        END IF
      END DO

C Calculation of the initial abundances:
      DO i = 1, n
        IF (k(i).NE.0) THEN
            print *,'Initially present is ',y(k(i))
            IF (y(k(i)) .EQ. 'GRAIN' .OR. y(k(i)) .EQ. 'GRAIN0') THEN
              IF (zone .NE. Nz) THEN
                y0(k(i)) = fraction(i)*density*ndust
              ELSE
                y0(k(i)) = fraction(i)*density
              END IF
            ELSE
              y0(k(i)) = fraction(i)*density
            END IF
        END IF
      END DO
      print *, 'initial abundances for zone', zone, 'done!'
      RETURN
      END

C..............................................................................
C
C This subroutine initializes  abundances of chemical species for t!=0.
C
C..............................................................................
C
C Input parameter(s):
C
C ny             == number of species
C
C y(ny)          == the names of non-conserved species,
C
C n                 == number of species with initial abundances
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
      DO i = 1, ny
        y0(i) = 0.0D0
      END DO
      DO i = 1, n
        k(i) = 0
      END DO
      print *,"line 1053"

C Search of species 'yr' among species 'y':
      i = 0

C A loop by 'yr':
10    CONTINUE
        i = i + 1
        IF (i.gt.n)  GO TO 40  ! checked all 'yr', exit,

C A loop by 'y':
        j = 0
20    CONTINUE
        j = j + 1
        IF (j.gt.ny) GO TO 10

C Found a necessary species:
        IF (yr(i).eq.y(j)) THEN
          k(i) = j
          GO TO 10  ! switch to the next species 'yr',
        END IF
        GO TO 20  ! switch to the next species 'y',
40    CONTINUE


C Calculate the initial abundance of species on grains and the initial
C grain abundance
      n_ice_init = 0
      print *,"line 1083"
      DO i = 1,n
        IF ('(gr)' .EQ. yr(i)(len_trim(yr(i))-3:)) THEN
          n_ice_init = n_ice_init + fraction(i)*density
        END IF
C scale grain abundance by dust settling parameter (except for midplane zone)
        IF ('GRAIN' .EQ. trim(yr(i))) THEN
          IF (zone .NE. Nz) THEN
            ngr_init = fraction(i)*density*ndust
          ELSE
            ngr_init = fraction(i)*density
          END IF
        END IF
      END DO

C Calculation of the initial abundances:
      print *,"line 1099"
      print *,n
      DO i = 1, n
        print *,i,k(1)
        IF (k(i).NE.0) THEN
C          print *,'Initially present is ',y(k(i))
          IF (y(k(i)) .EQ. 'GRAIN' .OR. y(k(i)) .EQ. 'GRAIN0') THEN
            IF (zone .NE. Nz) THEN
              y0(k(i)) = fraction(i)*density*ndust
            ELSE
              y0(k(i)) = fraction(i)*density
            END IF
          ELSE
              y0(k(i)) = fraction(i)*density
          END IF
            print *,y0(k(i))
        END IF
      END DO
      RETURN
      END

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
C na           == number of species
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
10    CONTINUE
      IF (i.gt.na) GOTO 20      ! No 'species' in 'array', stop,
      IF (species.eq.array(i)) GOTO 30 ! Found 'species', exit,
         i = i + 1
         GOTO 10
C 'species' was not found in 'array', stop:
20    CONTINUE
      index = 0
      RETURN
C 'Species' was found in 'array', initialize 'index':
30    CONTINUE
      index = i
      RETURN
      END
