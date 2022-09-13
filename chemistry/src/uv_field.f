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
	double precision ISM_field, I0(3)
	character*900 buff
	character*20 words(90)
	integer nwords, read_flg, j, i, inw, nj, sflg
	integer header, nlen, LyAj, fend
	integer testnum

c External Functions
	double precision atod

C Format(s)
 20		format(a900)

c Use the first radius
	radius = Rs(zone)

c Factors for analytical ISM field (Habing field)
	I0(1) = 3.2028e15
	I0(2) = -5.15442e18
	I0(3) = 2.0546e21

	maxzone = 0
	read_flg = 0
	sflg = 0
	header = 0
	i = 1

	open(unit = 99, file = trim(inpfile))
c	print *, 'READUV: opening photon file:'
	read(99,20,end=100) buff
	write(*,'(a70)') buff

C JF 12/8/10 - Read in until found radius is larger than given radius.  Then
C pick the radius that is closest to the given radius
	fend = 0
	radius_low = 0
	radius_high = 0
	do while (fend .EQ. 0)
		read(99, '(A900)', iostat=fend) buff
		if (index(buff, 'Radius(AU)') .ne. 0) then
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
	   read(99, '(A900)', iostat=fend) buff
	   call gwordslong(buff, nwords, words)

c  loop exit construct
		if ((words(1) .eq. 'Radius(AU)').and.(sflg .eq. 1)) then
c		    print *, 'readuv: done reading in photons'
		    nlam = i
		    exit
		end if

		if (words(1) .eq. 'Radius(AU)') then
			read_flg = read_flg + 1
			findrad = atod(words(2))
			if (findrad .EQ. radius_match) then
c				print *, 'readuv: isolated radius in file'
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

c Include cosmic-ray induced UV field (Shen 2004)
c This is the Draine field, G0 used elsewhere is from Habing which is a factor
c of 1.7 smaller.

c  Commenting out the following lines Mar 9, 2013 ( LIC)
c				uvfield(j,i) = uvfield(j,i) + 1e-4 * ISM_field
c     &					* CRatten(j)

c Set erroneous values to 0
				if (uvfield(j,i) .lt. 0.0) uvfield(j,i) = 0.0

C Find the zone with the maximum flux at the given wavelength (1000 A)
c j = zone, i = wavelength
				if (uvfield(j,i) .GT. UVtemp .AND. zAU(j)
     &					.LE. ZMAX) then
					UVtemp = uvfield(j,i)
					UVmaxzone(i) = j
				endif

                ! zeroing UV, remove later:
c                uvfield(j,i) = 0.0

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
	DOUBLE PRECISION G0eval, G0_old

	radius = Rs(zone)
	maxlam = 1500.0

c  OLD METHOD:
c	call find_etau(maxlam, etau, zone)
c	G0_old = g0_100*(100.0/radius)**2*etau
c	photorate_nocalc = g0_100*(100.0/radius)**2*alpha*etau
c	print *, 'PRATE NOCALC: ', r1, g0_100, radius, alpha, etau

C New method:
	call find_G0(maxlam, G0eval, zone, 1)
	photorate_nocalc = G0eval*alpha

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

c	print *,fnamecalc
C See if file of pre-calculated photorates exists.  If not, read in the xsect
C file, calculate and save the photorates to a file.

c commented out check if exist by LIC 3/9/2013
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
	double PRECISION gb
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

c	print *, 'Here now, ils'
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
c	if(jzone>40) then
c	 print *, 'Photuv=',photuv,jzone,lambda(i)
c	end if

C original version:  Commented Mar 9 , 2013
c	prate = photuv*sigma(i)*((100.0/radius)**2.0)
	prate = photuv*sigma(i)

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
	INTEGER pdtype, i, grspecindex, grindex
	DOUBLE PRECISION sigmapd, yield, flux
	DOUBLE PRECISION sigmafac, sigmagr
	DOUBLE PRECISION Rdust, Nsites, sitedens, ngrspec
	DOUBLE PRECISION Mlayers, Mlayersabun
	DOUBLE PRECISION fr,molarea,rgr

c	sigmafac = 1.0e-15 / 4.0
c	sigmagr = 3.0e-10
c	Rdust = 1e-5 * freezeeffic**(-1.0/1.5) ! In Angstrom
c	sigmagr = pi * Rdust**2.0  ! In Angstrom^2
c	sigmagr = sigmagr * freezeeffic  ! Total surface area of grains go down with grain growth.
C	Nsites = 1.0e6 * (Rdust/1e-5)**2  ! Number of sites scale with size of grain^2



	n_ice = 0.0
	photodesorp_rate = 0.0

C   Remove '(gr)' from the species name)
	nogrspec = specname(1:len_trim(specname)-4)

C   Calculate the abundance of on-grain species
	i = 0
	do while (grspecs(i) .NE. -1)
		IF (timestep .EQ. 1) THEN
			n_ice = n_ice + abundances(grspecs(i), timestep)
		ELSE
			n_ice = n_ice + abundances(grspecs(i), timestep-1)
		ENDIF
		i = i+1
	end do

c   Get the current abundance of grains
	CALL ispecies('GRAIN        ', ns, s, grindex)
	IF (timestep .EQ. 1) THEN
		ngr(zone) = abundances(grindex, timestep)
	ELSE
		ngr(zone) = abundances(grindex, timestep-1)
	END IF

C   Get the current abundance of the species we're looking at
	CALL ispecies(specname, ns, s, grspecindex)
	IF (timestep .EQ. 1) THEN
		ngrspec = abundances(grspecindex, timestep)
	ELSE
		ngrspec = abundances(grspecindex, timestep-1)
	END IF

c   Set the abundance for the first zone (when abundances = 0)
	IF (n_ice .LT. MINABUN) THEN
		n_ice = n_ice_init
	ENDIF
	IF (ngr(zone) .LT. MINABUN) THEN
		ngr(zone) = ngr_init
	ENDIF
	IF (ngrspec .LT. MINABUN) THEN
		ngrspec = MINABUN
	ENDIF

C   Default grain size = 0.1 micron or 1e5 angstrom.
C   As small grains are removed via epsilon (ndust), more ice per grain.

	Nsites = 1.0e6
	if (ndust .ne. 1.0) Nsites = Nsites*ndust

C   ngr propto rgr**-3.5
C   sig_gr propto rgr**2
C   Increasing rgr by 2.0 changes sigma by 4.0
C   but ngr goes down by a factor of 0.088
C   total change, freezeeffic = rgr**-1.5

	if (freezeeffic .ne. 1.0) then
c        sigmagr = pi*(rgr*fr
        Nsites = Nsites * locdust(zone) !sigmagr/(pi*rgr**2)  ! scale number of sites by increase/decrease in surface area per grain
c        Nsites = Nsites * (freezeeffic**(3.5/1.5))
	end if

c   Number of monolayers is nice / ngr / Nsites
c   Abundance of monolayers is then nice / Mlayers = ngr * Nsites
	Mlayers = n_ice / ngr(zone) / Nsites

	molarea = 3.14 * (1e-8/2.0) ** 2  ! 1 angstrom sites
	photodesorp_rate = yield * molarea !* freezeeffic !* 2.0  ! * ngrspec / n_ice ! temporary comment LIC

C Continuum rate
	IF (pdtype .EQ. -23) THEN
		photodesorp_rate = photodesorp_rate * totflux(zone)
C Lyman alpha rate
	ELSE IF (pdtype .EQ. -24) THEN
		do i=1,nwavl
			if (lambda(i) .EQ. 1210) then
				flux = 10.0 * uvfield(zone,i)
				exit
			endif
		enddo
		photodesorp_rate = photodesorp_rate * flux
	ELSE
		write(*,*) 'ERROR: illegal call to photodesorp_rate function'
		write(*,*) 'Reaction type must be -23 or -24'
		write(*,*) 'Reaction Type Given: ', pdtype
		photodesorp_rate = -1.0
		stop
	END IF

c Set maximum monolayers limit
	IF (Mlayers .GE. 2.0) THEN
		photodesorp_rate = photodesorp_rate * 2.0 / Mlayers
	ENDIF

	RETURN
	END

C..............................................................................
C Helper Functions
C..............................................................................
C..............................................................................
c GWORDS subdivides a string into words separated by spaces and tabs.
c This version should work for very long lines (~900
c characters long) in photons.dat with 60 zones.
C..............................................................................
        subroutine gwordslong(buff,nwords,word)
cimplicit undefined (a-z)
        implicit none
        character*900 buff
        character*20 word(90), q
        integer nwords
        integer blflg, i, ilet
        character*1 qd

        nwords = 0
        blflg = 1
        do 100 i = 1, 900
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

C..............................................................................
C Carbon Monolayer Calculation
C carbonmlayer
C..............................................................................
	SUBROUTINE carbonmlayer(fr)
	IMPLICIT NONE

C Common Blocks
	INCLUDE "Fcn.h"
	INCLUDE "environ.h"
	INCLUDE "constants.h"
	INCLUDE "rates.h"

	INTEGER n_mols,ml_n
	PARAMETER (n_mols = 5)
	CHARACTER*13 species
	CHARACTER*13 specarr(n_mols)
	INTEGER imol
	INTEGER r1index, grindex,p
	DOUBLE PRECISION nr1, Nummonol(n_mols)
	DOUBLE PRECISION polint,Nsites,c_ice_lyr,fr
	real nummonotot,l_diff

C Define Carbon Goo:
      specarr(1) ='CO(gr)       '
      specarr(2) ='CO2(gr)      '
      specarr(3) ='H2CO(gr)     '
      specarr(4) ='CH3OH(gr)    '
      specarr(5) ='CH4(gr)      '

c      print *,itime,zonedefT
c
c
	  fr = 1.0
      nummonotot = 0.0
	  Nsites = 1.0e+06
C Calculate grain species:
	  CALL ispecies('GRAIN        ', ns, s, grindex)
c	  write(*,*) grindex
      do imol = 1, n_mols
		species = specarr(imol)
		CALL ispecies(species, ns, s, r1index)

		IF (timestep .eq. 2) THEN
			nr1 = abundances(r1index, timestep)
			ngr(zone) = abundances(grindex, timestep)
		ELSE
			nr1 = abundances(r1index, timestep-1)
			ngr(zone) = abundances(grindex, timestep-1)
		END IF

		Nummonol(imol) = nr1/ngr(zone)/Nsites
		nummonotot = nummonotot + Nummonol(imol)
c		dataarr(imol,3) = Numr1(imol)
c		write(*,*) species,timestep,Nummonol(imol),ngr
	  enddo
	  c_ice_lyr = ANint(nummonotot)
c      if (timestep.ge. 2) write(*,*) 'Total N_c: ',c_ice_lyr

C	  Ice diffusion length from Oberg 2009:
	  l_diff = 0.6+0.024*Td(zone)
	  if (c_ice_lyr.gt.0.0) then
		fr = 0.0

		do ml_n = 1,4
			fr = fr + exp(-(c_ice_lyr+(ml_n-1))/l_diff)-
     +			exp(-(c_ice_lyr+1+(ml_n-1))/l_diff)
		enddo
	  endif
c	  if (timestep.ge. 2) write(*,*) 'fraction desorb:',fr
c
c
      return
      END


C..............................................................................
C
C This function reads in the isrf field from Kamber's photons file.  Mostly taken from
C Ted's readuv.f file.
C
C 5/14/13 -- KRS & LIC
C
C..............................................................................
C
C Input parameter(s):
C
C inpfile			== photons.dat filename for ISRF
C
C..............................................................................

	subroutine readisrf(inpfile)

C Include the environmental variables and xsect variables
	INCLUDE "environ.h"
	INCLUDE "rates.h"

  	CHARACTER*80 inpfile

	double precision findrad, radius, UVtemp
	double precision radius_low, radius_high, radius_match
	double precision nolya_scale, LyA_interp
	double precision ISM_field, I0(3)
	character*900 buff
	character*20 words(90)
	integer nwords, read_flg, j, i, inw, nj, sflg
	integer header, nlen, LyAj, fend
	integer testnum

c External Functions
	double precision atod

 20		format(a900)

	radius = Rs(zone)
	maxzone = 0
	read_flg = 0
	sflg = 0
	header = 0
	i = 1

	open(unit = 99, file = trim(inpfile))
c	print *, 'READUV: opening photon file:'
	read(99,20,end=101) buff
	write(*,'(a70)') buff

C JF 12/8/10 - Read in until found radius is larger than given radius.  Then
C pick the radius that is closest to the given radius
	fend = 0
	radius_low = 0
	radius_high = 0
	do while (fend .EQ. 0)
		read(99, '(A900)', iostat=fend) buff
		if (index(buff, 'Radius(AU)') .ne. 0) then
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

c	print *, 'rmatch!',radius_match,radius

C Now read in the UV field for the found radius
	open(unit = 99, file = trim(inpfile))
	read(99,20,end=101) buff


	fend = 0
	do while (fend .EQ. 0)
	   read(99, '(A900)', iostat=fend) buff
	   call gwordslong(buff, nwords, words)

c  loop exit construct
		if ((words(1) .eq. 'Radius(AU)').and.(sflg .eq. 1)) then
c		    print *, 'readuv: done reading in photons'
		    nlam = i
		    exit
		end if

		if (words(1) .eq. 'Radius(AU)') then
			read_flg = read_flg + 1
			findrad = atod(words(2))
			if (findrad .EQ. radius_match) then
c				print *, 'readuv: isolated radius in file'
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
				isrffield(j,i) = atod(words(inw))

c set the last zone uv field to be equal to the zone above it to avoid
c extrapolation errors.
c JF, 12/8/10 - Corrected so that it only is used if the last zone has a higher
c flux than the one immediately above it.
				if (inw .eq. nj .and. isrffield(j-1,i) .LT. isrffield(j,i)) then
					if (ndust .eq. 1) then
						isrffield(j,i) = isrffield(j-1,i)
					else
						isrffield(j,i) = 0.0
					endif
				endif

c Set erroneous values to 0
				if (isrffield(j,i) .lt. 0.0) isrffield(j,i) = 0.0
C  Now add the ISRF component to the stellar UVfield:
				uvfield(j,i) = uvfield(j,i)+isrffield(j,i)
		    end do

		    i = i + 1
		end if

	  end do
101	continue

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
		do j=1,nlam
			if (lambda(j) .NE. 1210) then
				totflux(i) = totflux(i) + isrffield(i,j) * 10
			endif
		enddo
	enddo

	return
	end
