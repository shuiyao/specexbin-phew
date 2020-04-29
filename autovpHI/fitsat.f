	subroutine fitsat(n,imin,imax)

c  fit a saturated line

	include 'vpdefs.h'
	integer n,imin,imax,ilow,ihi,i,k,ideconv
	double precision NHItrials(ntrials,2)
	double precision minchisqr,chisqr,bestchisqr
c
        do 7 i=imin,imax
          minflux(i) = minflux(i) - fsigma*noise(i)
 7      continue
c  deconvolve blended saturated lines into separate regions, if necessary
	ideconv = imax
	CALL deconv(n,imin,imax,ilow,ihi,ideconv)

c  try various b to see which has lowest chi-square
	bestchisqr = 1.d30
	do 2 k=1,4*ntrials
	  bpar(n) = bminsat + (bparmax-bminsat)*k/(4*ntrials)

	  CALL NHItest(n,imin,ideconv,ilow,ihi,NHItrials)

	  NHI(n) = NHItrials(ntrials/2,1)
	  CALL model(n,n,imin,ideconv,1)	! determines workflux
	  minchisqr = chisqr(MAX(imin,ilow-npixcusp),ilow) +
     &	    chisqr(ihi,MIN(ideconv,ihi+npixcusp))
C	write(6,'(2g12.5,5i10)') bpar(n),NHI(n),imin,ideconv,
C     &	ilow,ihi,npixcusp
	  if (minchisqr.le.bestchisqr) then
	    bestchisqr = minchisqr
	    bestbpar = bpar(n)
	  end if
 2	continue

c  determine best values
	bpar(n) = bestbpar
	CALL NHItest(n,imin,ideconv,ilow,ihi,NHItrials)
	NHI(n) = NHItrials(ntrials/2,1)

c  compute residual flux
	CALL model(1,n,1,ndata,1)
	  do 50 i=imin,ideconv
	    residold = resid(i)
	    resid(i) = 1.d0 + MAX(flux(i),0.d0) - workflux(i)
            if (resid(i).gt.1.d0) resid(i) = 1.d0
	    noise(i) = noise(i)*sqrt(MAX(resid(i),noise(i))/
     &	      (MAX(residold,noise(i))))
 50	  continue

	return 
	end

C=============================================================================

	subroutine NHItest(n,imin,imax,ilow,ihi,NHItrials)

c  brackets best NHI value to match cusp region for a given b param

	include 'vpdefs.h'
	integer n,imin,imax,ilow,ihi,jtries
	double precision NHItrials(ntrials,2),minchisqr
	double precision NHIlow,NHIhi,dlogNHI
	double precision chisqr

	NHIlow = NHImin
	NHIhi = NHImaxsat
	jtries = 0

c  try various NHI to see which has lowest chi-square
	do 5 j=1,5
C	write(6,*) log10(NHIlow),log10(NHIhi)
	    do 10 i=1,ntrials
	      dlogNHI = (log10(NHIhi) - log10(NHIlow))/ntrials
	      NHItrials(i,1) = 10**(log10(NHIlow) + dlogNHI*i)
	      NHI(n) = NHItrials(i,1)
	      CALL model(n,n,imin,imax,1)	! determines workflux
	      NHItrials(i,2) = 
     &		chisqr(MAX(imin,ilow-npixcusp/2),ilow+npixcusp/2) +
     &	        chisqr(ihi-npixcusp/2,MIN(imax,ihi+npixcusp/2))
c  disallow trial if the model goes 5-sigma below data anywhere
	  do 7 k=imin,imax
	    if(workflux(k).lt.minflux(k)-5.*noise(k)) NHItrials(i,2)=1.e30
 7	  continue
C	if(n.eq.4) write(6,*) i,log10(NHItrials(i,1)),NHItrials(i,2)

 10	    continue

c  bracket NHI which give lowest chisqr, going from largest to smallest NHI
	    NHIhi = NHItrials(ntrials,1)
	    NHIlow = NHItrials(ntrials-2,1)
	    minchisqr = NHItrials(ntrials,2)
	    do 20 i=ntrials-1,2,-1
	      if (NHItrials(i,2).gt.1.5*minchisqr) goto 5
	      if (NHItrials(i,2).lt.minchisqr) then
	        minchisqr = NHItrials(i,2)
	        NHIlow = NHItrials(i-1,1)
	        NHIhi = NHItrials(i+1,1)
	      end if
 20	    continue

 5	continue

	return 
	end

C=============================================================================

	double precision function chisqr(ilow,ihi)
	include 'vpdefs.h'
	integer i,ilow,ihi

        chisqr = 0.d0
        do 10 i=ilow,ihi
          chisqr = chisqr + ((resid(i)-workflux(i))/noise(i))**2
 10     continue

        return
	end

C=============================================================================

	subroutine deconv(n,imin,imax,idown,iup,ideconv)

	include 'vpdefs.h'
	integer imin,imax,n,i,iup,idown,j,nsatreg,ideconv
	double precision satflux(nmax)
	integer satreg(maxregions,2)
	integer ew_finder

c  identify region of saturation
	do 10 i=imin,centpix(n)
	  if (minflux(i).lt.fsigma*sigma(i)) then
	    idown = i
	    goto 12
	  end if
 10	continue
 12	continue
	do 15 i=imax,centpix(n),-1
	  if (minflux(i).lt.fsigma*sigma(i)) then
	    iup = i
	    goto 17
	  end if
 15	continue
 17	continue

c  identify if more than one saturated feature in detection region
	do 20 i=1,ndata
	  satflux(i) = 1.d0-resid(i)
 20	continue
	j = ew_finder(idown,iup,satreg,nsatreg,satflux)
c  if none found, do more foolproof test to be sure
	if (nsatreg.le.1) then
	  do 25 i=idown,iup
	    if( minflux(i).gt.4.*sigma(i) ) then
	      nsatreg = nsatreg+1
	      satreg(1,1) = idown+1
	      satreg(1,2) = i
	      satreg(2,1) = i+1
	      satreg(2,2) = iup-1
	      j = j+1
	      goto 27
	    end if
 25	  continue
 27	  continue
	end if

	nsatreg = 0
	do 30 i=1,j
	  if (satreg(i,1).gt.idown.AND.
     &	    satreg(i,2).lt.iup) nsatreg = nsatreg+1
 30	continue

	centpix(n) = (iup+idown)/2
	if (nsatreg.eq.0) return

c  found more than one crossing: deconvolve
	centpix(n) = (idown+satreg(1,1))/2
	ideconv = (satreg(1,1)+satreg(1,2))/2
	write(6,*) 'nregions:',j,vel(satreg(1,1)),vel(satreg(1,2))
	write(6,*) 'sregions:',nsatreg,vel(idown),vel(iup)
	write(6,'(a,3i6,3f10.1)') 'deconv: ',imin,ideconv,centpix(n),
     &	vel(imin),vel(ideconv),vel(centpix(n))

	return
	end

