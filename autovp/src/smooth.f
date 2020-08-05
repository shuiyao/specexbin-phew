
      subroutine smooth(imin,imax)
c
c  smooths the spectrum, finds minima, and subtracts noise to
c  to get minimum acceptable fit spectrum
c
	include 	'vpdefs.h'
	parameter (rt2piinv = 0.39894228d0)
	integer imin,imax
	double precision dx,fluxreal
	integer i,j,nsm
c
	do 10 i=imin,imax
	  minflux(i) = 0.
	  nsm = INT(sigmasm*nsmfact)
	  do 20 j=i-nsm,i+nsm
	    dx = 1.d0*(j-i)
	    if (j.lt.1) then
	      fluxreal = resid(j+ndata)
	    else if (j.gt.ndata) then
	      fluxreal = resid(j-ndata)
	    else
	      fluxreal = resid(j)
	    endif
c  Averaging
C	    minflux(i) = minflux(i) + fluxreal/(2*nsm+1)
c  Gaussian smoothing
	    minflux(i) = minflux(i) +
     &		rt2piinv*EXP(-0.5d0*(dx/sigmasm)**2)*fluxreal/sigmasm
 20	  continue
	  if (minflux(i).gt.1.d0) minflux(i) = 1.d0
 10	continue

        call fappend(ion_name,'min',sig_file)
        open(unit=2,file=sig_file,status='unknown')
        do 60 i=imin,imax
          write(2,'(i6,2f12.4,4f12.6)') 
     &	    i,waveln(i),vel(i),flux(i),
     &      minflux(i),resid(i),noise(i)
 60     continue
        close(2)

	return
	end

c=======================================================================

      subroutine smoothsat(imin,imax)
c
c  smooths the spectrum, finds minima, and subtracts noise to
c  to get minimum acceptable fit spectrum
c
	include 	'vpdefs.h'
	parameter (rt2piinv = 0.39894228d0)
	integer imin,imax
	double precision dx,fluxreal,maxnoise
	integer i,j,nsm
c
	do 5 i=imin,imax
	  maxnoise = MAX(maxnoise,noise(i))
 5	continue
	do 10 i=imin,imax
	  minflux(i) = 0.
	  sigmasm = 3.d0*(noise(i)/maxnoise)**2
	  nsm = INT(4.d0*sigmasm)
	  do 20 j=i-nsm,i+nsm
	    dx = 1.d0*(j-i)
	    if (j.lt.1) then
	      fluxreal = resid(j+ndata)
	    else if (j.gt.ndata) then
	      fluxreal = resid(j-ndata)
	    else
	      fluxreal = resid(j)
	    endif
c  Averaging
	    minflux(i) = minflux(i) + fluxreal/(2*nsm+1)
c  Gaussian smoothing
c	    minflux(i) = minflux(i) +
c     &		rt2piinv*EXP(-0.5d0*(dx/sigmasm)**2)*fluxreal/sigmasm
 20	  continue
	  if (minflux(i).gt.1.d0) minflux(i) = 1.d0
 10	continue

C        call fappend(ion_name,'min',sig_file)
C        open(unit=2,file=sig_file,status='unknown')
C        do 60 i=1,ndata
C          write(2,'(2f12.4,2f12.6)') waveln(i),vel(i),minflux(i),
C     &      noise(i)
C 60     continue
C        close(2)

	return
	end
