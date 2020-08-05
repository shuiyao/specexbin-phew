      integer function ew_finder(minpix,maxpix,detect,line,flx)
c
c  find the regions in the spectrum that are objective features
c  between pixel values minpax and maxpix,
c  above a confidence limit set by N_sigma (see vpdefs.h)
c
c  returns the min,max pixels of each detection region (detect)
c  and the number of regions detected (line, also = ew_finder)
c
c  modified: added ew limits spectrum calculations: Mon Feb 26 1996
c
	include 	'vpdefs.h'
	integer		 detect(maxregions,2)
	double precision  flx(*)
      integer            i,j,k,m,line,id,minpix,maxpix
      double precision   zabs
      double precision   lile(nmax),bige(nmax),lils(nmax),bigs(nmax)
      logical		 flag
c
c  compute the detection spectrum
      call fappend(ion_name,'ewlim',sig_file)
      open(unit=2,file=sig_file,status='unknown')
        lile(minpix) = 0.
        lils(minpix) = 0.
        lile(maxpix) = 0.
        lils(maxpix) = 0.
      do 13 i=minpix+1,maxpix-1
       lile(i) = 0.5*abs(waveln(i-1)-waveln(i+1))*(1.0d0-flx(i))
       lils(i) = 0.5*abs(waveln(i-1)-waveln(i+1))*noise(i)
       bige(i) = 0.0d0
       bigs(i) = 0.0d0
 13   continue
      do 14 i=minpix+nwidth,maxpix-nwidth
       do 17 j=i-nwidth,i+nwidth
        bige(i) = bige(i) + lile(j)
        bigs(i) = bigs(i) + lils(j)**2
 17    continue
       bigs(i) = sqrt(bigs(i))
       zabs = waveln(i)/lambda0
       write(2,201) waveln(i),vel(i),flx(i),noise(i),
     &              bige(i)/bigs(i),N_sigma*bigs(i)/zabs,
     &              0.5*N_sigma*bigs(i)/zabs
 14   continue
      close(unit=2)  
c
	do 60 m=1,maxregions
	  detect(m,1) = 0
	  detect(m,2) = 0
 60	continue
c
c  we are go for ion processing
c  scan for the detected regions
      flag = .false.
      line = 0
      do 27 i=minpix+nwidth,maxpix-nwidth
        if (.not.flag) then
         if (bige(i)/bigs(i).ge.N_sigma) then
          flag = .true.
          line = line + 1     
          detect(line,1) = i - nwidth
         end if
        else
         if (bige(i)/bigs(i).le.N_sigma) then
          flag = .false.
          detect(line,2) = i + nwidth
         end if
        end if
 27   continue
c
c  Handle case where detection region extends beyond boundary of spectrum
	do 50 m=1,line
	  if (detect(m,1).lt.1) detect(m,1) = 1
	  if (detect(m,2).lt.1) detect(m,2) = ndata
C	 write(6,*) 'pre',m,detect(m,1),detect(m,2)
 50	continue
c
c check regions for minimum size
      do 40 k=1,line
	if( detect(k,2)-detect(k,1).le.2 ) then
	  if( detect(k,1).gt.1) detect(k,1) = detect(k,1)-1
	  if( detect(k,2).lt.ndata) detect(k,2) = detect(k,2)+1
	else if( detect(k,2)-detect(k,1).eq.3 ) then
	  if( detect(k,1).gt.1) detect(k,1) = detect(k,1)-1
	endif
 40   continue
c
c  OK... now, the way this is written, some regions can overlap,
c  which is silly and redundant.  So, let's find 'em and combine 'em
c  while book-keeping properly
      m = 1
      do 31 while (m.lt.line)
        if (detect(m,2)+nwidth.ge.detect(m+1,1)) then  ! overlap=nwidth pixels
         detect(m,2) = detect(m+1,2)
         do 33 k=m+1,line
          detect(k,1) = detect(k+1,1)
          detect(k,2) = detect(k+1,2)
 33      continue
	 line = line-1
	else
	 m = m+1
        end if
 31   continue
c
c  now write the regions for record keeping and further processing
      open(unit=3,file="ew_regions.dat",status="unknown")
      do 29 i=1,line
C       write(3,301) i,1.25,waveln(detect(i,1)),waveln(detect(i,2)),
C     &              0.5*(waveln(detect(i,1))+waveln(detect(i,2)))
       write(3,301) i,detect(i,1),detect(i,2),
     &		waveln(detect(i,1)),waveln(detect(i,2)),
     &		    vel(detect(i,1)),vel(detect(i,2))
c	write(6,*) i,detect(i,1),detect(i,2)
 29   continue
      close(unit=3)
c
 201  format(1x,2f12.4,6f9.4)
 301  format(1x,i5,2i6,4f12.4)
c
	ew_finder = line

 1000 return
c
      end
c
