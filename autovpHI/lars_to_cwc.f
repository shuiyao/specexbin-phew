	program lars_to_cwc
C
C	Converts artificial spectra from Lars' format (list of 
C	optical depths at given redshift) to format readable by
C	cwc's PROFIT program.
C	Also shifts spectrum so that boundary is not within a
C	spectral feature (since simulation has periodic BC).
C
	parameter(c=300000.0,pi=3.14159265,npix=1000,rn=4.3)
	parameter (rt2piinv = 0.398942280)
	character*80 infile,outfile
	integer i,istart,nskip
	real Intensity(2*npix),v(2*npix),lambda(2*npix),Noise(2*npix)
	real tau,smooth(npix)
	real sn,z,counts,DeltaInt,select,Prob,maxint
	real RAN1
	logical start

	sn = 30.		! Signal to noise ratio, 10^6 -> infinity
	z = 3.0			! redshift of simulation box
	taufact = 0.725		! normalization factor for optical depth
	OPEN(UNIT=3,FILE='l2c.dat',STATUS='OLD',ERR=2)
	read(3,*,ERR=2) sn,z,taufact	! input values if input file exists
	CLOSE(3)
 2	continue
	nskip = 1
	write(6,*) sn,z,taufact
	np = npix
	deltav = 50.*sqrt(1+z)*22.2222/npix	! velocity across box
	infile = 'l2c.tmp'
	outfile = 'H1216N.raw'
	OPEN(UNIT=3,FILE=infile,STATUS='OLD')
	OPEN(UNIT=4,FILE=outfile,STATUS='UNKNOWN')
	maxint = 0.
	tau = RAN1(-1)	! seeds random number generator
	DO 10 i=1,npix
	    READ(3,*) tau
	    Intensity(i) = EXP(-tau/taufact)
	    v(i) = deltav*(i-1)
	    lambda(i) = 1215.6701*(1+z)*(1+v(i)/c)
	    counts = 0.5*sn*sn*(1+sqrt(1+4.*rn*rn/(sn*sn)))
	    Noise(i) = sqrt(counts*Intensity(i) + rn*rn)/counts
c            IF (Intensity(i).ge.0.95) start = .TRUE.
c            IF (.NOT.start) istart = i
	    IF (sn.lt.1e6) THEN
 5		DeltaInt = (RAN1(0)-0.5)*10.*Noise(i)
	    	Prob = 1./(Noise(i)*sqrt(2.*pi)) * exp (-0.5*
     &		    (DeltaInt/Noise(i))**2.)
		select = RAN1(0)/(Noise(i)*sqrt(2.*pi))
		if (Prob.gt.select) Intensity(i) = 
     &		    MAX(Intensity(i)+DeltaInt,0.)
		if (Prob.le.select) GOTO 5
	    END IF
	    IF (Intensity(i).gt.maxint) THEN
		maxint = Intensity(i)
		istart = i
	    END IF
 10	CONTINUE
C
C  Shift data to have max at ends
	DO 15 i=1,npix
	  Intensity(i+npix) = Intensity(i)
	  Noise(i+npix) = Noise(i)
 15	continue
	DO 20 i=1,npix
	  Intensity(i) = Intensity(i+istart-1)
	  Noise(i) = Noise(i+istart-1)
 20	CONTINUE
C
C  Interpolate so that spacing is deltav
	deltal = 0.06	! HIRES resolution (Hu et al)
	deltav = deltal*c/(1+z)/1215.6701
	dv0 = v(2)-v(1)
	if (deltav.gt.dv0) then
	  vnext = v(1)
	  iv = 1
	  vmax = v(npix)
	  do 40 while (vnext.le.vmax-deltav)
	    vnext = vnext + deltav
	    iv = iv+1
	    ibelow = INT(vnext/dv0)+1
	    Intensity(iv) = ((vnext-v(ibelow))*Intensity(ibelow+1) +
     &		(v(ibelow+1)-vnext)*Intensity(ibelow))/
     &		(v(ibelow+1)-v(ibelow))
	    Noise(iv) = ((vnext-v(ibelow))*Noise(ibelow+1) +
     &		(v(ibelow+1)-vnext)*Noise(ibelow))/
     &		(v(ibelow+1)-v(ibelow))
	    lambda(iv) = ((vnext-v(ibelow))*lambda(ibelow+1) +
     &		(v(ibelow+1)-vnext)*lambda(ibelow))/
     &		(v(ibelow+1)-v(ibelow))
C	write(6,*) ibelow,lambda(ibelow),lambda(iv),lambda(ibelow+1)
	    v(iv) = vnext
 40	  continue
	  np = iv+1
	  Intensity(np) = Intensity(1)
	  Noise(np) = Noise(1)
	  v(np) = v(iv)+deltav
	  lambda(np) = lambda(iv)+deltal
	end if
C
C  Smooth if desired
        sigmasm = 0.0*nskip
	if (sigmasm.ne.0.0) then
        nsm = INT(sigmasm*3.0)
        do 60 i=1,np
          smooth(i) = 0.
          do 70 k=i-nsm,i+nsm
            j = k
            if (j.lt.1) j=j+np
            if (j.gt.np) j=j-np
            dx = 1.0*(j-i)
            smooth(i) = smooth(i) +
     &          rt2piinv*EXP(-0.50*(dx/sigmasm)**2)*Intensity(j)/sigmasm
 70       continue
 60     continue
	else
	  do 80 i=1,np
	    smooth(i) = Intensity(i)
 80	  continue
	end if
C
C  output result
	write(6,'(a,2f10.3,5i10)') 
     &	  'l2c: z,sn,istart,nskip=',z,sn,istart,nskip,np
	DO 120 i=1,np,nskip
	    WRITE(4,'(2f10.4,2f10.6)') 
     &		lambda(i),v(i),smooth(i),Noise(i)
 120	CONTINUE
	CLOSE(4)
C
	STOP
	END
	
