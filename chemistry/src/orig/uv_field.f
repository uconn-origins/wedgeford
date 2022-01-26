C..............................................................................
C
C This function reads in the uv field from the photons file.  Mostly taken from
C Ted's readuv.f file.  
C
C..............................................................................
C
C Input parameter(s):
C
C inpfile			== photons.dat filename
C
C..............................................................................
	subroutine readuv(inpfile)

C Include the environmental variables and xsect variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

	CHARACTER*80 inpfile

	double precision findrad, radius, UVtemp
	double precision radius_low, radius_high, radius_match
	double precision nolya_scale, LyA_interp
	character*700 buff
	character*20 words(90)
	integer nwords, read_flg, j, i, inw, nj, sflg
	integer header, nlen, LyAj, fend

c External Functions
	double precision atod
	
C Format(s)
 20		format(a700)

c Use the first radius 
	radius = Rs(zone)

	maxzone = 0
	read_flg = 0
	sflg = 0
	header = 0
	i = 1

	open(unit = 99, file = trim(inpfile))
	print *, 'READUV: opening photon file:'
	read(99,20,end=100) buff
	write(*,'(a70)') buff

C JF 12/8/10 - Read in until found radius is larger than given radius.  Then 
C pick the radius that is closest to the given radius
	fend = 0
	radius_low = 0
	radius_high = 0
	do while (fend .EQ. 0)
		read(99, '(A700)', iostat=fend) buff
		if (index(buff, 'Radius(AU)')) then
			call gwordslong(buff, nwords, words)
			findrad = atod(words(2))
			if (findrad .GE. radius) then
				radius_high = findrad
				exit
			else
				radius_low = findrad
			endif
		else
			cycle
		endif
	enddo
	close(99)

C Find radius closest to given radius
	if ((radius - radius_low) .LT. (radius_high - radius)) then
		radius_match = radius_low
	else
		radius_match = radius_high
	endif

c Match last radius if within 1%
	if ( (radius_high .EQ. 0) .AND. 
     +		((radius - radius_low)/radius .LE. 1e-2) ) then
		radius_match = radius_low
c Match first radius if within 1%
	else if ( (radius_low .EQ. 0) .AND. 
     +		((radius_high - radius)/radius .LE. 1e-2) ) then
		radius_match = radius_high
	endif

C Now read in the UV field for the found radius	
	open(unit = 99, file = trim(inpfile))
	read(99,20,end=100) buff

	fend = 0
	do while (fend .EQ. 0)
	   read(99, '(A700)', iostat=fend) buff
	   call gwordslong(buff, nwords, words)

c  loop exit construct
		if ((words(1) .eq. 'Radius(AU)').and.(sflg .eq. 1)) then
		    print *, 'readuv: done reading in photons'
		    nlam = i
		    exit
		end if

		if (words(1) .eq. 'Radius(AU)') then
			read_flg = read_flg + 1 
			findrad = atod(words(2))
			if (findrad .EQ. radius_match) then
				print *, 'readuv: isolated radius in file'
				sflg = 1 
			end if 
		end if

		if (sflg .eq. 1) header = header + 1

		if ((sflg .eq. 1).and.(header .gt. 3)) then
			lambda(i) = atod(words(1))
			UVtemp = 0
			nj = nwords
			Itemp = 0.0
		    do inw = 2, nj
				j = inw - 1
				uvfield(j,i) = atod(words(inw))

c set the last zone uv field to be equal to the zone above it to avoid
c extrapolation errors.  
c JF, 12/8/10 - Corrected so that it only is used if the last zone has a higher 
c flux than the one immediately above it.
				if (inw .eq. nj .and. uvfield(j-1,i) .LT. uvfield(j,i)) then
					if (ndust .eq. 1) then
						uvfield(j,i) = uvfield(j-1,i)
					else
						uvfield(j,i) = 0.0
					endif
				endif
				if (uvfield(j,i) .lt. 0.0) uvfield(j,i) = 0.0

C Find the zone with the maximum flux at the given wavelength (1000 A)
c j = zone, i = wavelength
				if (uvfield(j,i) .GT. UVtemp .AND. zAU(j) .LE. ZMAX) then
					UVtemp = uvfield(j,i)
					UVmaxzone(i) = j
				endif
				
		    end do
		    i = i + 1
		end if
	  end do
100	continue
	if (sflg .eq. 0) then
		print *, 'ERROR in reading photons.dat file'
		print *, 'unable to find radius in photons.dat file.  Exiting'
		stop
	end if

	nlam = nlam - 1
	nj = nj - 1

	close(99)

C calculate the integrated flux at each zone
C each wavelength bin is 10 angstroms wide
C only used for normal photodesorption, so remove LyA flux
	do i=1,nj
		totflux(i) = 0.0
		do j=1,nlam
			if (lambda(j) .NE. 1210) then
				totflux(i) = totflux(i) + uvfield(i,j) * 10
			endif
		enddo
	enddo


C if include_lya is false, remove LyA from the uv field and scale the 
C UV field so that the total flux is constant
C 8/28/09 - Remove the UV scaling, adding UV flux that we don't want
C i = zone, j = wavelength
	IF (.NOT. include_lya) then
		do j=1,nlam
			if (lambda(j) .EQ. 1210) then
				LyAj = j
				exit
			endif
		enddo
		
		do i=1,nj
c Ignore zones with no flux (outside of 400 x 400 AU box for Tom's
c photons.dat files		
			if (totflux(i) .EQ. 0) then
				cycle
			endif

			LyAinterp = (uvfield(i,LyAj-1) + uvfield(i,LyAj+1))/2.0
			uvfield(i,LyAj) = LyAinterp
			totflux(i) = totflux(i) + 10 * uvfield(i, LyAj)

c			nolya_scale = (totflux(i) + 
c     +			10.0*(uvfield(i,LyAj)-LyAinterp)) / totflux(i)
c			do j=1,nlam
c				uvfield(i,j) = uvfield(i,j) * nolya_scale
c			enddo

		enddo
	ENDIF

	return
	end


C..............................................................................
C
C PHOTORATE_NOCALC
C This function returns the photorate for the reactions that DON'T calculate
C it from the UV field and cross section files.
C
C..............................................................................
C
C Input parameter(s):
C
C r1			== First reactant name
C
C alpha			== alpha value from rreacs file
C
C..............................................................................
	DOUBLE PRECISION FUNCTION photorate_nocalc(r1, alpha)
	IMPLICIT NONE

C Include the environmental variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

	CHARACTER*13 r1
	DOUBLE PRECISION alpha
	
	CHARACTER*80 post, pre, fname
	CHARACTER*80 radstring, radd
	DOUBLE PRECISION radius
	LOGICAL fexist
	
	DOUBLE PRECISION maxlam, etau, prate

c Use the first radius 
	radius = Rs(zone)

	maxlam = 1500.0
	call find_etau(maxlam, etau, zone)
c	print *, 'PRATE NOCALC: ', r1, g0_100, radius, alpha, etau
	photorate_nocalc = g0_100*(100.0/radius)**2*alpha*etau

C Output results to .prate files
c	pre = '/Nirgal1/fogel/astrochem/semenov/prate/'
c	post = '.'//trim(dust)//'.nocalc'

c	IF(radius >= 100.0) THEN
c        write(radstring, '(f7.3)') radius
c        radd = radstring(1:3)
c	ELSE IF(radius >= 10.0) THEN
c        write(radstring, '(f6.3)') radius
c        radd = radstring(1:4)
c	ELSE
c        write(radstring, '(f5.3)') radius
c        radd = radstring(1:5)
c	END IF

c	fname = trim(pre)//trim(r1)//'.'//trim(radd)//trim(post)

c	INQUIRE(file=fname, exist=fexist)
c    IF (.not. fexist) THEN
c		OPEN(unit=11, file=fname, status='new')
c	ELSE
c		OPEN(unit=11, file=fname, status='old', access='append')
c	ENDIF

c	write(11,900) zone, photorate_nocalc
c 900	format(1x, 'nocalc photorate for zone ', i3, ' = ', e24.18)

c	CLOSE(11)

	RETURN
	END


C..............................................................................
C
C PHOTORATE
C This function returns the photorate, either calculated or read in from a 
C file with previously calculated rates
C
C..............................................................................
C
C Created 04/14/2008
C
C..............................................................................
C
C Input parameter(s):
C
C r1			== First reactant name
C
C..............................................................................
	DOUBLE PRECISION FUNCTION photorate(r1, rtype)
	IMPLICIT NONE

C Include the environmental variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

	CHARACTER*13 r1
	INTEGER rtype
	
	CHARACTER*80 add, add2, preadd, radstring, radd
	CHARACTER*80 fnamexs, fnamecalc, path
	DOUBLE PRECISION radius
	LOGICAL xsfexist
	INTEGER i

c Use the first radius 
	radius = Rs(zone)
	write(radstring, '(f7.3)') radius

c	preadd = '/Nirgal1/fogel/astrochem/semenov/xsect/'
	preadd = 'xsect/'
	add = '.photoxs'
	add2 = '.'//trim(dust)//'.calc'
	
C Add _# to file name if there are multiple paths.  Path # = rtype - 29
	IF(rtype.GT.30) THEN
		write(path, '(I1)') rtype-29
		path = '_'//trim(path)
	ELSE
		path = ''
	END IF

	fnamexs = trim(preadd)//trim(r1)//trim(path)//trim(add)

	fnamecalc = trim(preadd)//trim(r1)//trim(path)//'.'
	IF(radius >= 100.0) THEN
        write(radstring, '(f9.5)') radius
	ELSE IF(radius >= 10.0) THEN
        write(radstring, '(f8.5)') radius
	ELSE
        write(radstring, '(f7.5)') radius
	END IF
	radd = radstring
	
	fnamecalc = trim(fnamecalc)//trim(radd)//add2

C See if file of pre-calculated photorates exists.  If not, read in the xsect
C file, calculate and save the photorates to a file.  
	INQUIRE(file=fnamecalc, exist=xsfexist)
	IF (.not. xsfexist) THEN
		CALL readxs(fnamexs)
		CALL calc_prate(fnamecalc, r1)
	END IF
	
C Open file and read in appropriate rate
	OPEN(unit=10, file=fnamecalc, status='old')
	DO i=1,zone
		read(10,100) photorate
 100	format(26x, e24.18) 
	END DO

	close(10)

	RETURN
	END

C..............................................................................
C
C READXS
C Read in the data from the .xsect files and save them in the appropriate
C arrays
C
C..............................................................................
C Input parameter(s):
C
C fname			== file name of xsect file
C
C..............................................................................
	SUBROUTINE readxs(fname)

C Include common variables for xsect calculations
	INCLUDE "rates.h"

	character*80 fname, line
	double precision a, b, c
	integer ib, i, fend

C Initialize all variables to 0
	nband_read = 0
	do i=1,NBANDS
		siglam_l(i) = 0
		siglam_u(i) = 0
		sigma(i) = 0
	end do
	fend = 0
	ib = 0

	open(unit = 100, file = fname)
	
C Skip lines with a "c" as the first character (or if the line is blank)
	do while (fend .EQ. 0)
		read(100, '(a80)', iostat=fend) line
		if (len_trim(line) .GT. 0 .AND. line(1:1) .NE. 'c') then
			read(line,*) a,b,c
			siglam_l(ib) = a
			siglam_u(ib) = b 
			sigma(ib) = c*1.e-18
			ib = ib+1
		endif
	enddo

C number of wavelength bands read in
	nband_read = ib-1
	
	close(100)
 	
	RETURN
	END

C..............................................................................
C CALC_PRATE
C From the cross section file, calculate the rates for every zone
C
C..............................................................................
C Input parameter(s):
C
C fname			== file name to save rates to
C
C..............................................................................
	SUBROUTINE calc_prate(fname, r1)

C Include the environmental variables and xsect variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

	DOUBLE PRECISION photorate(Nz), prate
	EXTERNAL prate
	double precision ulim,llim, sim
	integer specno, j
	character*80 fname
	CHARACTER*13 r1

	open(unit = 100, file = fname, status='new')

	llim = 925.0
	ulim = 1590.0

C if the previous photorate was less than 1e-25, set the current rate to 0.  
C This saves computation time by not needing to calculate very small rates.
C 9/4/09 - Removed because causing errors when first zone had no UV field
C because it was outside the 400x400 AU box (no photorate was ever
C being calculated)
	do j = 1, Nz
c		if(j.gt.1 .and. photorate(j-1).le.1e-25) then
c			photorate(j) = 0
c		else
			call qsimp(prate,llim,ulim,sim,j)
			photorate(j) = sim
c		endif

c		print *, 'photorate (', j, ') = ', photorate(j)
		write(100,900) j, photorate(j)
	end do

 900	format(1x, 'photorate for zone ', i3, ' = ', e24.18)
	close(100)
		
	return
	end 

C..............................................................................
C PRATE - used in trapzd
C Returns the value of the cross section at the input point x
C..............................................................................
	DOUBLE PRECISION FUNCTION prate(x,jzone)

C Include the environmental variables and xsect variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

	double precision x, photuv, dphotuv, dlam
	integer ilam, jzone, i

	prate = 0

C Use the first radius 
	radius = Rs(zone)

C  first determine the UV photon intensity at x, then interpolate
      do i = 1, nlam
		if(lambda(i) .gt. x) then
			ilam = i
			exit
		end if
	end do

	dlam = lambda(i) - lambda(i-1)
	dphotuv = ((x-lambda(i-1))/dlam) * 
     +		  (uvfield(jzone,i)-uvfield(jzone,i-1))
	photuv = uvfield(jzone,i-1) + dphotuv
	
c	print *, i, jzone
c	print *, lambda(i-1), uvfield(jzone,i), uvfield(jzone,i-1)

C  now determine the vaule of the cross-section within the wavelength
C  range in angstroms
	do i = 0, nband_read
		if ((x .gt. siglam_l(i)) 
     +  .and. (x .le. siglam_u(i))) then
			exit
		end if
	end do

	prate = photuv*sigma(i)*((100.0/radius)**2.0)

	RETURN
	END


C..............................................................................
C
C PHOTODESORP_RATE
C This function returns the photodesorption rate (both continuum and lyman alpha)
C
C Calculations from Hassel et al. 2008 in prep
C..............................................................................
C
C Input parameter(s):
C
C yield = number of adsorbed particles ejected per incident photon
C theta = fractional surface coverage
C r1	= reactant
C rtype	= reaction type (24 = lyman alpha, or 23 = continuum)
C
C..............................................................................
C H2O, CO2, CO, CH4, NH3: theta = 1
C Other: theta = 1e-2
C
C Continuum: yield = beta
C Lyman alpha: yield from file
C..............................................................................
	DOUBLE PRECISION FUNCTION photodesorp_rate(specname, pdtype, yield)
	IMPLICIT NONE

C Include the environmental variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"
	INCLUDE "Fcn.h"
	INCLUDE "constants.h"
	
	CHARACTER*13 specname, nogrspec
	INTEGER pdtype, i, LyAwave, grspecindex, grindex
	DOUBLE PRECISION sigmapd, yield, flux
	DOUBLE PRECISION sigmafac, sigmagr
	DOUBLE PRECISION Rdust, Nsites, sitedens, ngrspec
	DOUBLE PRECISION ngr, Mlayers, Mlayersabun
	
c surface density of sites - Barlow p. 411
c	sigmafac = 1.0e-15 / 4.0
c	sigmagr = 3.0e-10
	Rdust = 1000e-8
	Nsites = 1.0e6
	sigmagr = pi * Rdust**2.0 
	LyAwave = 1216

	n_ice = 0.0
	photodesorp_rate = 0.0
	
C remove '(gr)' from the species name)
	nogrspec = specname(1:len_trim(specname)-4)

C Calculate the abundance of on-grain species
	i = 0
	do while (grspecs(i) .NE. -1)
		IF (timestep .EQ. 1) THEN
			n_ice = n_ice + abundances(grspecs(i), timestep)
		ELSE
			n_ice = n_ice + abundances(grspecs(i), timestep-1)
		ENDIF
		i = i+1
	end do

c Get the current abundance of grains
	CALL ispecies('GRAIN        ', ns, s, grindex)
	IF (timestep .EQ. 1) THEN
		ngr = abundances(grindex, timestep)
	ELSE
		ngr = abundances(grindex, timestep-1)
	END IF

C get the current abundance of the species we're looking at
	CALL ispecies(specname, ns, s, grspecindex)
	IF (timestep .EQ. 1) THEN
		ngrspec = abundances(grspecindex, timestep)
	ELSE
		ngrspec = abundances(grspecindex, timestep-1)
	END IF

c set the abundance for the first zone (when abundances = 0)
	IF (n_ice .LT. MINABUN) THEN
		n_ice = n_ice_init
	ENDIF
	IF (ngr .LT. MINABUN) THEN
		ngr = ngr_init
	ENDIF
	IF (ngrspec .LT. MINABUN) THEN
		ngrspec = MINABUN
	ENDIF

c Number of monolayers is nice / ngr / Nsites
c Abundance of monolayers is then nice / Mlayers = ngr * Nsites
	Mlayers = n_ice / ngr / Nsites

	photodesorp_rate = yield * sigmagr/Nsites * 2.0 * ngrspec / n_ice

C Continuum rate
	IF (pdtype .EQ. -23) THEN
        photodesorp_rate = photodesorp_rate * totflux(zone)
c		photodesorp_rate = totflux(zone) * yield * sigmagr/Nsites
c		photodesorp_rate = photodesorp_rate * 2.0 * ngrspec / n_ice
c		print *, 'totflux(zone): ', zone, totflux(zone)
C Lyman alpha rate
	ELSE IF (pdtype .EQ. -24) THEN
c find the flux at lyman alpha wavelength (1216 A) (subtract the interpolation off)
		do i=1,nwavl
			if (lambda(i) .EQ. 1210) then
				flux = 10.0 * uvfield(zone,i)
				exit
			endif
		enddo

        photodesorp_rate = photodesorp_rate * flux
c		photodesorp_rate = flux * yield * sigmagr/Nsites
c		photodesorp_rate = photodesorp_rate * 2.0 * ngrspec / n_ice

C ERROR
	ELSE
		write(*,*) 'ERROR: illegal call to photodesorp_rate function'
		write(*,*) 'Reaction type must be -23 or -24'
		write(*,*) 'Reaction Type Given: ', pdtype
		photodesorp_rate = -1.0
		stop
	END IF

c Set maximum monolayers limit
	IF (Mlayers .GE. 2) THEN
		photodesorp_rate = photodesorp_rate * Nsites * ngr / ngrspec 
	ENDIF

	RETURN
	END

C..............................................................................
C Helper Functions
C..............................................................................
C..............................................................................
c GWORDS subdivides a string into words separated by spaces and tabs.
c This version should work for very long lines (~700
c characters long) in photons.dat with 60 zones.
C..............................................................................
        subroutine gwordslong(buff,nwords,word)
cimplicit undefined (a-z)
        implicit none
        character*700 buff
        character*20 word(90), q
        integer nwords
        integer blflg, i, ilet
        character*1 qd

        nwords = 0
        blflg = 1
        do 100 i = 1, 700
			qd = buff(i:i)
			if (qd .eq. ' ' .or. qd .eq. '  ') then
				if (nwords .gt. 0 .and. blflg .eq. 0) then
					word(nwords) = q(1:ilet)
				endif
				blflg = 1
				go to 100
			endif
			if (blflg .eq. 1) then
				nwords = nwords + 1
				blflg = 0
				ilet = 1
			else
				ilet = ilet + 1
                endif
                q(ilet:ilet) = qd
100     continue
        if (blflg .eq. 0) word(nwords) = q(1:ilet)

        return
        end

C..............................................................................
C ATOD converts a given string to a double precision
C..............................................................................
        double precision function atod(string)
        implicit none
        character*80 string
        character*80 dummy

        write(dummy,90) string
90      format(a80)
        read(dummy,*) atod

        return
        end
C..............................................................................
C qsimp routine taken from NR - Integrate a function using Simpson's Rule
C "Returns the integral of the function from a to b.  The constants EPS can be 
C  set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1
C  is the maximum allowed number of steps.  Integration is performed by 
C  Simpson's rule."
C..............................................................................
      SUBROUTINE qsimp(func,a,b,s,jzone)
      INTEGER JMAX, jzone
      DOUBLE PRECISION a,b,func,s,EPS
	  EXTERNAL func
      PARAMETER (EPS=1.e-4, JMAX=25)
      INTEGER j
      DOUBLE PRECISION os,ost,st
      ost=0.
      os=0.
      do 11 j=1,JMAX
        call trapzd(func,a,b,st,j,jzone)
        s=(4.*st-ost)/3.
        if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or.(s.eq.0..and.os.eq.0.)) return
        endif
        os=s
        ost=st
11    continue
	  write(*,*) 'ERROR: too many steps in qsimp'
	  STOP
      END

C..............................................................................
C trapzd routine taken from NR
C "Routine implementing the extended trapezoidal rule"
C..............................................................................
      SUBROUTINE trapzd(func,a,b,s,n,jzone)
      INTEGER n, jzone
      DOUBLE PRECISION a,b,s,func
	  EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x

      if (n.eq.1) then
        s=0.5*(b-a)*(func(a,jzone)+func(b,jzone))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x,jzone)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif

      return
      END
