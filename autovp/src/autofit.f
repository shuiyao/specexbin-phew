      program autofit
c
c  automated routine to approximately fit Voigt profiles to data
c
	include 'vpdefs.h'
	include 'const.dek'
	integer i
c
c  get the ion for this run
        call getarg(1,ion_name)
	CALL readdata
c
c  fit Voigt profiles in each detected region
	CALL fitregion
c
c  combine small, nearby lines; also eliminate non-ion lines and doublet pairs
	CALL combine
C
        nelems = 1
        vhelio = 0.d0
        call fappend(ion_name,'pro',sig_file)
        open(unit=2,file=sig_file,status='unknown')
	write(2,'(1x,i5,5f14.6)') nelems,redz,vhelio
	write(2,'(1x,a10,i5)') ion_name,nlines
        do 70 i=1,nlines
          write(2,'(i6,f14.6,f16.6,f14.6,f10.3,1x,f10.3,1x,f10.3)')
     &      i,NHI(i)/1.d13,vel(centpix(i)),bpar(i),
     &      rho(centpix(i)),temp(centpix(i)),met(centpix(i))
c          write(2,'(i6,f14.6,f16.6,f14.6)')
c     &      i,NHI(i)/1.d13,vel(centpix(i)),bpar(i)
 70     continue
        close(2)

 1000	stop
	end
c
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
