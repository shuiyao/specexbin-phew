
	SUBROUTINE readdata

c  reads in data in the following format (compatible with PROFIT):
c
c  spectral info (file called ion_name//'.cln'):
c      wavelength   velocity   flux   noise
c
c  also sets workflux (the model) to be 1
c  and resid (the residual flux = flux-workflux+1) to be = flux
c 
        include 'vpdefs.h'
        include 'const.dek'
	integer i,j
        double precision gamma
	character*10 ion_str
c
c  get ion info from ions.dat
c        open(unit=1,file='/home/rad/lyalpha/autovp/ions.dat',
	open(unit=1,file='./ions.dat',
     &	  status='old')
	do 10 j=1,1000
          read(1,'(a10,f9.4,f7.4,g12.3,e9.4)',END=15)
     &      ion_str,lambda0,fosc,gamma,lambda1
	  if (ion_str.eq.ion_name) goto 20
 10	continue
 15	continue
	write(6,*) 'Could not find ion',ion_name
 20	continue
	gam8 = gamma/1.d8
        con1 = 1.0d5*fosc*sqrt(pi)*e*e*lambda0**2 / (me*c*c)
        con2 = (gam8*lambda0**2)/(4.0d0*pi*c)
	close(1)

c  read in data
      call fappend(ion_name,'cln',cln_file)
      open(unit=1,file=cln_file,status='old')
	avenoise = 0.d0
      ndata = 0
      do 5 i=1,nmax
       read(1,*,end=6) waveln(i),vel(i),flux(i),sigma(i),rho(i),
     & temp(i),met(i)
       ndata = ndata + 1
       resid(i) = flux(i)
       workflux(i) = 1.d0
       noise(i) = sigma(i)
       avenoise = avenoise+noise(i)
  5   continue
  6   close(unit=1)
      avenoise = avenoise/ndata
C      do 90 i=1,ndata
C	vel(i) = ckms*(waveln(i)/waveln(1)-1.d0)
C	if(i.le.10) write(6,*) i,waveln(i),vel(i)
C 90   continue
	redz = (waveln(1)/(lambda0*(1.d0+vel(1)/ckms)))-1.d0
        write(6,'(a,a10,4f12.4)') 'ion: ',ion_name,redz,lambda0,
     &	waveln(1),vel(1)
c
c  get user-set parameters
	CALL inparam
	if( chisqbad.eq.0.0 ) then
	  chisqbad = ((1./avenoise)*0.5)**0.25-0.25
	  if( chisqbad.gt.2.5 ) chisqbad = 2.5
	  if( chisqbad.lt.1.0 ) chisqbad = 1.0
	  chisqgood = 0.5*chisqbad
	  write(6,*) 'New chisqbad = ',chisqbad,avenoise
	end if
	if( ion_name(1:1).ne.'H' ) then
	  bparmax = bparmax/4.
	  bminsat = bminsat/4.
	end if

	return
	end

	integer function len1(a)
	character*10 a
	len1 = 0
	do while(a(len1:len1).ne.' ')
	  len1 = len1+1
	end do
	return
	end
