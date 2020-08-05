
	subroutine inparam()

	include 'vpdefs.h'


C***************** input parameters for autovp ********************************
C 1.0	fsigma -- number of noise vectors subtracted from smooth to get minflux
C 1.d15	NHImaxns -- NHI maximum for non-saturated line
C 1.d17	NHImaxsat -- NHI maximum for saturated line (MUST HAVE NO WINGS)
C 100.	bparmax -- b param maximum
C 10.	bminsat -- b param minimum for saturated line
C 0.98	freduce -- factor by which b and N are reduced during fitting
C 0	fitflag -- =0 for Gaussian, =1 for Voigt profiles, =2 for doublet Gauss.
C 2.0	sigmasm -- width of Gaussian smoothing in pixels (smooth)
C 0.5	bfact -- factor which sets max vel diff to combine lines (combine)
C 10	nsplitmin -- # of lines above which to split region
C 3	nlinemin -- min # of lines in split region
C 8.0	N_sigma -- detection threshold
C 100	maxiter -- maximum number of iterations for Marquardt minimization
C 1.d8	alamdamax -- maximum alamda in Marquardt minimization
C 2.0	chisqbad -- chisq above which it's a "bad fit", i.e. try 1-par min
C 1.2	chisqgood -- chisq below which it's a "good fit"
C 0.001	chisqtol -- if chisq reduces by less than this factor, it's converged
C 0.1	faccept -- for tossing lines, if chisq goes up by less than this, accept
C 0.4	fdoublet -- if doublet pos. has less than fdoublet as much flux, discard

c	open(UNIT=3,FILE='/disks/virgo2/oppenheimer/ionfiles/autovp.par',
	open(UNIT=3,FILE='./autovp.par',
     & STATUS='OLD')
	read(3,*)
	read(3,*) fsigma
	read(3,*) NHImaxns
	read(3,*) NHImaxsat
	read(3,*) bparmax
	read(3,*) bminsat
	read(3,*) freduce
	read(3,*) fitflag
	read(3,*) sigmasm
	read(3,*) bfact 
	read(3,*) nsplitmin
	read(3,*) nlinemin
	read(3,*) N_sigma
	read(3,*) maxiter
	read(3,*) alamdamax
	read(3,*) chisqbad
	read(3,*) chisqgood
	read(3,*) chisqtol
	read(3,*) faccept
	read(3,*) fdoublet
	CLOSE(3)

	write(6,*) 'fitflag= ',fitflag
	npixcusp = 1+INT(20./(vel(2)-vel(1)))	! 20 km/s cusp region
	nwidth = 1+INT(10./(vel(2)-vel(1)))	! 20 km/s detection region
	nsmfact = 4
	nint = 8
	minint = 4
	fvarymin = 0.01

	return
	end
