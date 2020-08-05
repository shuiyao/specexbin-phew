	program tipsy2profit
C
C	Converts artificial spectra from TIPSY format (list of 
C	vel & optical depths at given redshift) to format readable by
C	cwc's PROFIT program.
C	Also shifts spectrum so that boundary is not within a
C	spectral feature (since simulation has periodic BC).
C
	parameter(c=2.99792458e5,pi=3.14159265,npixmax=12000)
	parameter (rt2piinv=0.398942280, one=1.0)
	character*80 infile,outfile
	integer i,j,istart,nskip,npix,iseed
	real Intensity(2*npixmax),v(2*npixmax)
	real lambda(2*npixmax),Noise(2*npixmax)
	real tau,Ism(npixmax),vsm(npixmax)
	real sn,z,counts,DeltaInt,select,Prob,maxint
	real lambda0,sigmasm,totwt,totI
	real meanflux,meanforig,ftau
	real RAN1
	logical start

	sn = 20.		! Default Signal/noise ratio, 10^6 -> infinity
	z = 3.0			! default redshift of simulation box
	taufact = 1.0		! default normalization factor for optical depth
	deltal = 0.0		! default pixel size (<0 -> in km/s)
	lambda0 = 1215.6701	! default transition rest wavelength
        sigmasm = 0.0		! Default Gaussian FWHM smoothing, in km/s
	iflag = 0		! 0=Do nothing, 1=shift data to max, 2=truncate
				!	flux at max, 3=shift & truncate,
				!	4= shift, truncate, and restore flux
	iseed = -1
	OPEN(UNIT=3,FILE='t2p.dat',STATUS='OLD',ERR=2)
	read(3,*,ERR=2) sn,z,taufact,deltal,lambda0,sigmasm,iflag,iseed
	CLOSE(3)
 2	continue
	if( deltal.eq.0.0 ) then
	if( z.ge.2.0 ) then
	    deltal=0.06	! set wavelength spacing: HIRES
	else if( z.ge.1.5 ) then
	    deltal=0.06	! set wavelength spacing: intermediate
	else
	    deltal=0.06 ! set wavelength spacing: FOS
	end if
	end if
	rnfact = 1./1.4 ! parametrized by Q0000 data
	rn = rnfact/sn
	write(6,'(6f10.4,i5)') sn,z,taufact,deltal,lambda0,sigmasm,iflag
	infile = 't2p.tmp'
	outfile = 'H1216N.raw'
	OPEN(UNIT=3,FILE=infile,STATUS='OLD')
	OPEN(UNIT=4,FILE=outfile,STATUS='UNKNOWN')
	maxint = 0.
	tau = RAN1(iseed)	! seeds random number generator
	DO 10 i=1,npixmax
	    READ(3,*,END=12) v(i),tau
	    if( i.eq.1 ) v0 = v(1)
	    v(i) = v(i)-v0
	    Intensity(i) = EXP(-tau*taufact)
	    lambda(i) = lambda0*(1+z)*(1+v(i)/c)
 10	CONTINUE
 12	CONTINUE
	npix = i-1
	np = npix
C
C  Interpolate so that spacing is approximately deltal
	dv0 = v(2)-v(1)
	vbox = v(npix)+dv0
	if( deltal.lt.0.0 ) deltav = -deltal
	if( deltal.gt.0.0 ) deltav = deltal*c/(1+z)/lambda0
	np = 1+INT(vbox/deltav)
	deltav = vbox/np
C
C  Bin to desired velocity spacing
	if (deltav.gt.dv0 ) then
C	write(6,*) deltav,deltal*c/(1+z)/lambda0,dv0,vbox,np
          nsm = 1+INT(0.5*deltav/dv0)
          do 45 i=1,np
            vsm(i) = deltav*(i-1)
	    ibeg = 1+INT((vsm(i)-0.5*deltav-v(1))/dv0)
	    if( ibeg.le.1 ) ibeg = ibeg-1
	    iend = 1+INT((vsm(i)+0.5*deltav-v(1))/dv0)
	    Ism(i) = 0.0
            do 50 k=ibeg,iend
              j = k
              if (j.lt.1) j=j+npix
              if (j.gt.npix) j=j-npix
C	      dv = ABS(v(j)-vsm(i))
C	      if( dv.gt.0.5*vbox ) dv = vbox-dv
C	      dv = dv-0.5*deltav
C	      if( dv.gt.0.0 ) dv = dv0-dv
C	      if( dv.le.0.0 ) dv = 0.0
              if( k.eq.ibeg) then
		dv = vsm(i)-0.5*deltav
		if( dv.lt.0.0 ) dv = dv+vbox
		dv = ABS(v(j+1)-dv)
		if( dv.gt.0.5*vbox ) dv = vbox-dv
		fraclo = dv/dv0
		Ism(i) = Ism(i) + fraclo*Intensity(j)
	      else if( k.eq.iend ) then
		dv = vsm(i)+0.5*deltav
		if( dv.gt.vbox ) dv = dv-vbox
		dv = ABS(dv-v(j))
		if( dv.gt.0.5*vbox ) dv = vbox-dv
		frachi = dv/dv0
                Ism(i) = Ism(i) + frachi*Intensity(j)
	      else
		Ism(i) = Ism(i) + Intensity(j)
	      end if
C	if(i.eq.14) write(6,*) k,Ism(i),v(j),Intensity(j)
 50       continue
	    Ism(i) = Ism(i)/(fraclo+frachi+iend-ibeg-1)
C	write(6,*) i,Ism(i),fraclo,frachi,(fraclo+frachi+iend-ibeg-1)
 45       continue
	else
	  write(6,*) 'WARNING: t2p.tmp has insufficient resolution',
     &	    dv0,deltav
	  np = npix
	  do 55 i=1,np
	    vsm(i) = v(i)
	    Ism(i) = Intensity(i)
 55	  continue
	end if
C
C  Smooth with Gaussian resolution element if desired
	if (sigmasm.ne.0.0) then
        nsm = INT(sigmasm*3.0/deltav)
        do 60 i=1,np
          do 70 k=i-nsm,i+nsm
	    j = k
	    if( k.lt.1 ) j = k+np
	    if( k.gt.np ) j = k-np
            dx = vsm(j)-vsm(i)
	    wt = EXP(-0.50*(dx/sigmasm)**2)
	    totwt = totwt+wt
            totI = totI+Intensity(j)*wt
 70       continue
	  Ism(i) = totI/totwt
 60     continue
	end if
C
C  Fit "continuum" by chopping all flux above highest pixel, rescaling fluxes
	meanforig = 0.
	do i=1,np
	  meanforig = meanforig + Ism(i)
	end do
	meanforig = meanforig / np
	do i=1,np
	  if( maxint.lt.Ism(i) ) maxint = Ism(i)
	end do
	if (iflag.eq.0.OR.iflag.eq.1) maxint = 1.	! turn this off
	do i=1,np
	  Ism(i) = Ism(i) / maxint
	end do

C  Add back in chopped flux, by adjusting optical depths (if desired)
	if( iflag.eq.4 ) then
	  do i=1,100
	    meanflux = 0.
	    do j=1,np
	      meanflux = meanflux + Ism(j)
	    end do
	    meanflux = meanflux / np
	    ftau = 1.+(meanflux-meanforig)/meanforig
	    if( ftau.lt.0.5 ) ftau = 0.5
	    if( ftau.gt.2.0 ) ftau = 2.0
	    write(6,*) "meanf,orig,dmean=",meanflux,meanforig,abs(ftau-1.)
	    if( abs(ftau-1.).lt.1.e-3 ) goto 95
	    do j=1,np
	      Ism(j) = exp(ftau*LOG(Ism(j)))
	    end do
	  end do
 95	  continue
	end if
C  Add noise
	DO 100 i=1,np
	    shotnoise = 0.0
	    if( Ism(i).gt.0.0 ) shotnoise = sqrt(Ism(i))/sn
C	    Noise(i) = sqrt(shotnoise*shotnoise+rn*rn)
C     &		/sqrt(1.+rnfact*rnfact)
	    Noise(i) = 1./sn
	    IF (sn.lt.1e6) THEN
 5		DeltaInt = (RAN1(0)-0.5)*10.*Noise(i)
	    	Prob = 1./(Noise(i)*sqrt(2.*pi)) * exp (-0.5*
     &		    (DeltaInt/Noise(i))**2.)
		select = RAN1(0)/(Noise(i)*sqrt(2.*pi))
		if (Prob.gt.select) Ism(i) = 
     &		    Ism(i)+DeltaInt
		if (Prob.le.select) GOTO 5
C	write(6,*) Ism(i)-DeltaInt,Ism(i),Noise(i),DeltaInt
	    END IF
 100	CONTINUE
C
C  Shift data to have max at ends
	istart = 0
	maxint = 0.0
	if (iflag.eq.1.OR.iflag.eq.3.OR.iflag.eq.4) then
	DO i=1,np
	    IF (Ism(i).gt.maxint) THEN
		maxint = Ism(i)
		istart = i
	    END IF
	END DO
	DO i=1,np
	  Ism(i+np) = Ism(i)
	  Noise(i+np) = Noise(i)
	END DO
	DO i=1,np
	  Ism(i) = Ism(i+istart-1)
	  Noise(i) = Noise(i+istart-1)
	END DO
	endif
	meanflux = 0.
	do j=1,np
	  meanflux = meanflux + Ism(j)
	end do
	meanflux = meanflux / np
C
C  compute wavelengths and output result
	write(6,'(a,2f8.2,2i7,2f8.3)') 
     &	  'l2c: z,sn,istart,np,fmax,<flux>=',z,sn,istart,np,maxint,meanflux
	DO 120 i=1,np
	    lambda(i) = lambda0*(1+z)*(1+vsm(i)/c)
	    WRITE(4,'(2f12.4,2f12.6)') 
     &		lambda(i),vsm(i),Ism(i),Noise(i)
 120	CONTINUE
C
C  Add one wrapped pixel for better endpoint constraint
	vsm(np+1) = vsm(np)+vsm(2)-vsm(1)
	lambda(np+1) = lambda0*(1+z)*(1+vsm(np+1)/c)
	WRITE(4,'(2f12.4,2f12.6)') 
     &		lambda(np+1),vsm(np+1),Ism(1),Noise(1)
C
C  If too few pixels, pad end with continuum to help AutoVP
	if( np.lt.10 ) then
	  DO 130 i=np+2,10
	    vsm(i) = vsm(i-1)+vsm(2)-vsm(1)
	    lambda(i) = lambda0*(1+z)*(1+vsm(i)/c)
	    WRITE(4,'(2f12.4,2f12.6)') lambda(i),vsm(i),one,Noise(1)
 130	  CONTINUE
	endif

	CLOSE(4)
C
	STOP
	END
	
