      program maxcont
c
c  fits a continuum level by recursive iteration fit of 4th-order
c  spline polynomial.
c
	implicit double precision (a-h,o-z)
      integer            nmax
      parameter          (nmax = 5000)
      integer            i,j,k,ndata
      double precision   wave(nmax),vel(nmax),splinefit(nmax),
     &                   signorm(nmax),fluxnorm(nmax)
      double precision   lowlim,hilim,chisq,chisqold
      double precision	 sigmax,fluxmax,velmax,flux1
      character*80       ion_name,sig_file,cln_file
c
c  get the ion info for this run
      call getarg(1,ion_name)
      call fappend(ion_name,'raw',sig_file)
      open(unit=1,file=sig_file,err=1000,status='old')
	write(6,*) 'FITTING CONTINUUM: ',sig_file(1:40)
c  read in data
      ndata = 0
      sigmax = 0.d0
	fluxmax = 0.d0
      do 5 i=1,nmax
       read(1,*,end=6) wave(i),vel(i),fluxnorm(i),signorm(i)
       sigmax = MAX(sigmax,signorm(i))
       fluxmax = MAX(fluxmax,fluxnorm(i))
	if (fluxmax.eq.fluxnorm(i)) ifmax = i
       ndata = ndata + 1
  5   continue
  6   close(unit=1)
      velmax = vel(ndata)
	flux1 = 0.0
	do 10 i=1,5
	  flux1 = flux1+fluxnorm(i)
 10	continue
	flux1 = flux1/5-sigmax
	flux1 = fluxmax-2.5*sigmax
	if(flux1.gt.1.d0) flux1 = 1.d0
	fluxdec = 0.0
	do 15 i=1,ndata
	  splinefit(i) = flux1
	  fluxdec = fluxdec+(1.-splinefit(i))
 15	continue
	fluxdec = fluxdec/ndata
	fluxmax = fluxmax-1.5*sigmax
	if(fluxmax.gt.1.d0) fluxmax = 1.d0
c
c  normalize spectra and noise level to continuum, and output
 50   call fappend(ion_name,'cln',cln_file)
      open(unit=2,file=cln_file,status='unknown')
      do 20 i=1,ndata
	fluxnorm(i) = fluxnorm(i)/splinefit(i)
	signorm(i) = signorm(i)/splinefit(i)
	write(2,'(2f12.4,2f12.6)') 
     &	 wave(i),vel(i),fluxnorm(i),signorm(i)
 20   continue
      close(unit=2)  
	write(6,'(3x,a,f12.4)') 'Flux decrement: ',fluxdec
c
 1000 stop
c
      end

c=======================================================================
c
c  
c
c
      subroutine          fappend(infile,delim,outfile)   
c
c
c
c-----------------------------------------------------------------------
c
c  replace the delimeter on infile with delim and output the
c  result in outfile
c  if the infile has no delimeter, then append delim prefixed 
c  with a "." to infile and output the result in outfile
c
c=======================================================================
      integer             i,k,lend      
      character*(*)       infile,delim,outfile
c
       lend = len(infile)
       k    = 0
       do 09 i=1,lend
        k = i
        if ((infile(i:i).eq.'.').or.(infile(i:i).eq.' ')) goto 10
 09   continue
c
 10   if (delim.eq.' ') then
       outfile = infile(1:k-1)
      else
       if (infile(k:k).eq.'.') outfile = infile(1:k)//delim
       if (infile(k:k).eq.' ') outfile = infile(1:k-1)//"."//delim
      end if
c
c
c  normal exit
       return
c
       end
c
c=======================================================================
