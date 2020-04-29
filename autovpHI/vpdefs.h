c	include file vpdefs.h

c	parameter declarations file for autofit
c
c	ndata 		# of data points (pixels) in spectrum
c	nregions	# of regions in which features are detected
c	waveln		wavelength of each pixel in data set
c	vel		velocity of each pixel in data set
c	flux		normalized flux of each pixel in data set
c	sigma		noise of each pixel in data set
c	minflux		minimum acceptable smoothed flux for model at each pixel
c	workflux	model flux 
c	resid		residual flux = flux - workflux at each pixel
c	noise		sigma on resid flux
c	dfdv		derivative of minflux wrt v (when minflux is polyn. fit)
c	d2fdv2		second derivative of minflux wrt v
c	region		list of min,max of pixels demarcating detection regions
c	vmin,vmax	min,max velocity of region
c	NHI,dNHI	column densities for components, errors
c	bpar,dbpar	b parameters for compenents, errors
c	vline,dvline	velocity of line center for compenents, errors
c	centpix		central pixel value for components
c	nlines		number of total components in model
c	velfit		velocity array within fitting region (vmin->vmax)
c	lambda0		rest wavelength of transition
c	fosc		oscillator strength of transition
c	gam8		gamma of transition (in 10^8 units)
c	redz		redshift of data
c	con1,con2	constants used in modelling Voigt profiles
c	ion_name	name of ion
c	cln_file,sig_file	filenames

      implicit double precision (a-h,o-z)
      double precision	NHImin,bparmin,tolderiv
      parameter		(nmax = 80000, maxregions = 5000, maxlines = 50)
      parameter		(maxparm = 3*maxlines)
      parameter		(NHImin = 0.005, bparmin = 1.0d0)
      parameter		(ntrials = 20)
      parameter		(tolderiv = 0.01)
      integer		ndata,nregions
      double precision  waveln,vel,flux,noise,sigma
      double precision  temp,rho,met
      double precision  minflux,workflux,resid,dfdv,d2fdv2
      integer		region
      double precision	NHI,bpar,vline,dNHI,dbpar,dvline,ewline
      double precision  tempout, rhoout, metout
      integer		centpix,lineinfo,nlines
      double precision	lambda0,fosc,gam8,redz,con1,con2,lambda1
      double precision	fsigma,avenoise
      double precision	velfit,parm0,deltap
      double precision	vmin,vmax
      character*80      cln_file,sig_file
      character*10      ion_name

	COMMON/datacom/waveln(nmax),vel(nmax),flux(nmax),sigma(nmax),
     &			temp(nmax),rho(nmax),met(nmax),ndata
	COMMON/workcom/minflux(nmax),workflux(nmax),resid(nmax),
     &			dfdv(nmax),noise(nmax),d2fdv2(nmax),
     &			fsigma,avenoise
	COMMON/regioncom/region(maxregions,2),vmin,vmax,nregions
	COMMON/ioncom/lambda0,lambda1,fosc,gam8,redz,con1,con2,ion_name
	COMMON/filecom/cln_file,sig_file
	COMMON/linecom/NHI(maxregions),bpar(maxregions),
     &			vline(maxregions),dNHI(maxregions),
     &			dbpar(maxregions),dvline(maxregions),
     &			centpix(maxregions),lineinfo(maxregions),
     &			ewline(maxregions),
     &                  tempout(maxregions),rhoout(maxregions),
     &                  metout(maxregions),nlines
	COMMON/fitcom/velfit(nmax),deltap(maxparm),parm0(maxparm)

C***************** input parameters for autovp ********************************
C 1.0	fsigma -- number of noise vectors subtracted from smooth to get minflux
C 1.d15	NHImaxns -- NHI maximum for non-saturated line
C 1.d17	NHImaxsat -- NHI maximum for saturated line (MUST HAVE NO WINGS)
C 100.	bparmax -- b param maximum
C 10.	bminsat -- b param minimum for saturated line
C 0.98	freduce -- factor by which b and N are reduced during fitting
C 0	fitflag -- =0 for Gaussian, =1 for Voigt profiles, =2 for doublet Gauss
C 10	npixcusp -- number of pixels for cusp region
C 3	nwidth -- half-width of detection window (also # of overlap pix)
C 3.0	sigmasm -- width of Gaussian smoothing in pixels (smooth)
C 4	nsmfact -- number of smoothing factors over which to smooth (smooth)
C 0.5	bfact -- factor which sets max vel diff to combine lines (combine)
C 6	nint -- # of pixels in 1 interval of piecewise polyn fit (polyfit)
C 3	minint -- minimum nint; also overlap between intervals
C 10	nsplitmin -- # of lines above which to split region
C 3	nlinemin -- min # of lines in split region
C 0.01	fvarymin -- sets min factor by which parm is varied (oneparmin,nparmin)
C 8.0	N_sigma -- detection threshold
C 100   maxiter -- maximum number of iterations for Marquardt minimization
C 1.d8  alamdamax -- maximum alamda in Marquardt minimization
C 2.0   chisqbad -- chisq above which it's a "bad fit", i.e. try 1-par min
C 1.2   chisqgood -- chisq below which it's a "good fit"
C 0.001 chisqtol -- if chisq reduces by less than this factor, it's converged
C 0.1   faccept -- for tossing lines, if chisq goes up by less than this, acceptC 0.4   fdoublet -- if doublet pos. has less than fdoublet as much flux, discard

	double precision fisgma,NHImaxns,NHImaxsat,bparmax,bminsat
	double precision freduce,sigmasm,bfact,fvarymin,N_sigma
	double precision alamdamax,chisqbad,chisqgood,
     &			chisqtol,faccept,fdoublet
	integer nwidth,nsmfact,nint,minint,nsplitmin,nlinemin
	integer npixcusp,maxiter,fitflag

	COMMON/paramcom/fisgma,NHImaxns,NHImaxsat,bparmax,bminsat,
     &			freduce,sigmasm,bfact,fvarymin,N_sigma,dfdvsat,
     &			alamdamax,chisqbad,chisqgood,chisqtol,faccept,
     &			fdoublet,ipivint,nwidth,nsmfact,nint,minint,
     &			nsplitmin,nlinemin,npixcusp,maxiter,fitflag

