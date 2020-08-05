	program clipsy2profit
C
C	Converts artificial spectra from FOGGY format (list of 
C	vel & optical depths at given redshift) to format readable by
C	cwc's PROFIT program.
C	Also shifts spectrum so that boundary is not within a
C	spectral feature (since simulation has periodic BC).
C
	parameter(c=300000.0,pi=3.14159265,npix=1000)
	parameter (rt2piinv = 0.398942280)
	character*80 infile,outfile,ionfile
	integer i,istart,nskip
	real Intensity(2*npix),v(2*npix),lambda(2*npix),Noise(2*npix)
	real tauHI,IntHI(2*npix)
	real tau,smooth(npix)
	real sn,z,counts,DeltaInt,select,Prob,maxint
	real lambda0
	real RAN1
	character*10 ion_name,ion_str
	logical start

	ion_name = 'H1216N'	! default ion name
	sn = 30.		! Default Signal/noise ratio, 10^6 -> infinity
	z = 3.0			! default redshift of simulation box
	taufact = 1.0		! default normalization factor for optical depth
	OPEN(UNIT=3,FILE='c2p.dat',STATUS='OLD',ERR=2)
	read(3,'(a10)',ERR=2) ion_name	! input values if file exists
	read(3,*,ERR=2) sn,z,taufact
	CLOSE(3)
 2	continue
	deltal=0.044	! HIRES value
	rn = 1./(1.4*sn)
	nskip = 1
	np = npix

C  Get ion rest wavelength
c	ionfile = '/home/rad/lyalpha/autovp/ions.dat'
	ionfile = './ions.dat'
	OPEN(UNIT=2,FILE=ionfile,STATUS='OLD')
        do 4 j=1,1000
          read(2,'(a10,f9.4)',END=5) ion_str,lambda0
          if (ion_str.eq.ion_name) goto 6
 4      continue
 5      continue
        write(6,*) 'Could not find ion',ion_name
	stop
 6      continue
	CLOSE(2)

C  Open files for input/output
	infile = 'c2p.tmp'
	j=0
	do 7 i=1,10
	  if (ion_name(i:i).ne.' ') j=j+1
 7	continue
	outfile = ion_name(1:j)//'.raw'
	OPEN(UNIT=3,FILE=infile,STATUS='OLD')
	OPEN(UNIT=4,FILE=outfile,STATUS='UNKNOWN')
	write(6,'(a20,4f10.4)') outfile,sn,z,taufact,deltal
	outfile = ion_name(1:j)//'.nf'
	OPEN(UNIT=5,FILE=outfile,STATUS='UNKNOWN')
C
C  Generate intensities from optical depths, adding noise if desired
	maxint = 0.
	tau = RAN1(-1)	! seeds random number generator
	DO 10 i=1,npix
	    READ(3,*) v(i),tau,tauHI
	    Intensity(i) = EXP(-tau/taufact)
	    IntHI(i) = EXP(-tauHI/taufact)
	    IF (IntHI(i).gt.maxint) THEN
		maxint = IntHI(i)
		istart = i
	    END IF
	    IntHI(i) = Intensity(i)
	    lambda(i) = lambda0*(1+z)*(1+v(i)/c)
C	    counts = 0.5*sn*sn*(1+sqrt(1+4.*rn*rn/(sn*sn)))
C	    Noise(i) = sqrt(counts*Intensity(i) + rn*rn)/counts
            shotnoise = sqrt(Intensity(i))/sn
            Noise(i) = sqrt(shotnoise*shotnoise+rn*rn)
	    IF (sn.lt.1e6) THEN
 8		DeltaInt = (RAN1(0)-0.5)*10.*Noise(i)
	    	Prob = 1./(Noise(i)*sqrt(2.*pi)) * exp (-0.5*
     &		    (DeltaInt/Noise(i))**2.)
		select = RAN1(0)/(Noise(i)*sqrt(2.*pi))
		if (Prob.gt.select) Intensity(i) = 
     &		    MAX(Intensity(i)+DeltaInt,0.)
		if (Prob.le.select) GOTO 8
	    END IF
 10	CONTINUE
C
C  Shift data to have max at ends
	DO 15 i=1,npix
	  Intensity(i+npix) = Intensity(i)
	  IntHI(i+npix) = IntHI(i)
	  Noise(i+npix) = Noise(i)
 15	continue
	DO 20 i=1,npix
	  Intensity(i) = Intensity(i+istart-1)
	  IntHI(i) = IntHI(i+istart-1)
	  Noise(i) = Noise(i+istart-1)
 20	CONTINUE

C  Interpolate so that spacing is deltal in wavelength
	deltav = deltal*c/(1+z)/lambda0
	dv0 = v(2)-v(1)
	if (deltav.gt.dv0) then
	  vnext = v(1)
	  iv = 1
	  vmax = v(npix)
	  do 40 while (vnext.le.vmax-deltav)
	    vnext = vnext + deltav
	    iv = iv+1
	    ibelow = INT(vnext/dv0)
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
	    WRITE(5,'(2f10.4,2f10.6)') 
     &		lambda(i),v(i),IntHI(i),Noise(i)
 120	CONTINUE
	CLOSE(4)
C
	STOP
	END
	
