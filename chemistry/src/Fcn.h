C..............................................................................
C
C Global parameter(s):
C
C nreac       == maximal amount of chemical reactions to be read,
C
C nspec       == maximal amount of chemical species to be considered,
C
C nspec2      == maximal amount of initially present chemical species,
C
C ntime       == maximal amount of taken time steps,
C
C nss         x nspec == ~ maximal amount of non-zero Jacobian elements
C
C MINABUN     == minimum abunance allowed for network
C
C..............................................................................
C
C Common block(s):
C
C..............................................................................
C
C B1:
C
C ns           == amount of species,
C
C s(ns)        == a species set,
C
C nr           == amount of reactions in the input file,
C
C reacindex(nr)    == indices of the corresponding chemical reactions,
C
C r1(nr)       == names of the first reactants,
C
C ir1(nr)      == array of indices of 'r1' species in extended species set 's',
C
C r2(nr)       == names of the second reactants,
C
C ir2(nr)      == array of indices of 'r2' species in extended species set 's',
C
C p1(nr)       == names of the first products,
C
C ip1(nr)      == array of indices of 'p1' species in extended species set 's',
C
C p2(nr)       == names of the second products,
C
C ip2(nr)      == array of indices of 'p2' species in extended species set 's',
C
C p3(nr)       == names of the third products,
C
C ip3(nr)      == array of indices of 'p3' species in extended species set 's',
C
C p4(nr)       == names of the fourth products,
C
C ip4(nr)      == array of indices of 'p4' species in extended species set 's',
C
C p5(nr)       == names of the fifth products,
C
C ip5(nr)      == array of indices of 'p5' species in extended species set 's',
C
C alpha(nr)    == first components of the rate coeffs.,
C
C beta(nr)     == second components of the rate coeffs.,
C
C gamma(nr)    == third components of the rate coeffs.,
C
C ak(nr)       == rate coefficients, computed from 'alpha', 'beta' and 'gamma',
C
C tprint	   == testing variable so that rates are only printed once per Fcn
C
C rtype(nr)	   == type of reaction, from input file
C
C testspec	   == Species name to print out formation/destruction rates of
C
C ratetest	   == Whether to print out ratetest.out file with formation/destruction rates
C
C shieldtest   == Whether to print out selfshield.test file with test date from self shielding code
C
C..............................................................................
C Global parameter(s):
      INTEGER nreac, nspec, nspec2, ntime, nss
      PARAMETER (nreac=60000,nspec=5000,nspec2=5000,ntime=999,nss=nspec)

	  DOUBLE PRECISION MINABUN
	  PARAMETER (MINABUN=1.0e-30)

C Common block (rates, abundances, etc.):
      INTEGER ns, nre, reacindex, ir1, ir2, ip1, ip2, ip3, ip4,
     1  ip5, nJ, tprint, rtype, timestep, grspecs
      REAL*8 alpha, beta, gamma, ak, times, abundances
      CHARACTER*13 s, r1, r2, p1, p2, p3, p4, p5
      DIMENSION s(nspec), r1(nreac), r2(nreac), p1(nreac), p2(nreac),
     1  p3(nreac), p4(nreac), p5(nreac), reacindex(nreac), alpha(nreac),
     2  beta(nreac), gamma(nreac), ir1(nreac), ir2(nreac), ip1(nreac),
     3  ip2(nreac), ip3(nreac), ip4(nreac), ip5(nreac), times(ntime),
     4  ak(nreac),rtype(nreac),abundances(nspec,ntime),grspecs(nspec)



      COMMON /BL1/ ns, s, grspecs, nre, reacindex, ir1, ir2, ip1, ip2,
     1  ip3, ip4, ip5, r1, r2, p1, p2, p3, p4, p5, alpha, beta, gamma,
     2  ak, times, nJ, tprint, rtype, abundances, timestep


C..............................................................................
C ratetest		== whether rout files are created or not
C
C shieldtest	== whether shieldtest files are created or not
C
C maxratetest	== Maximum number of species to produce rout files for
C
C Nratetest		== Number of species in flags file to create rout files for
C
C fnameRT		== filename of rout file
C
C RTspec(Nratetest)		== Array of species to print rout files for
C
C RToutput(Nratetest)	== array of lines for rout files (rewritten per time step)
C
C..............................................................................
	INTEGER maxratetest
	PARAMETER (maxratetest=25)

	INTEGER Nratetest, NRTreacs(maxratetest)
	INTEGER RTreacs(maxratetest, nreac)
	DOUBLE PRECISION RTrates(maxratetest, nreac)
	LOGICAL ratetest, shieldtest
	CHARACTER*80 fnameRT(maxratetest)
	CHARACTER*13 RTspec(maxratetest)

	COMMON /RTEST/ RTrates, RTreacs, NRTreacs, Nratetest, ratetest,
     1  shieldtest, fnameRT, RTspec
