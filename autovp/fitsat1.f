	subroutine fitsat(n,imin,imax)

c  fit a saturated line

	include 'vpdefs.h'
	integer n,imin,imax,ilow,ihi
	logical shifted,saturated,deconv
	integer decreaseN,modelbelow,adjustb

c  deconvolve blended saturated lines into separate regions
	if (.NOT.deconv(n,imin,imax,ilow,ihi,ipivot))
     &	  CALL fitlim(n,imin,imax,ilow,ihi,ipivot)
c  first guess is HUGE
	NHI(n) = NHImaxsat
	bpar(n) = MIN(bparmax,vel(imax)-vel(imin))
	CALL model(n,n,imin,imax,1)	! determines workflux
c  decrease NHI while adjusting b to match at pivot points
	do 10 while (decreaseN(n,imin,imax,ilow,ihi).ne.0)
	  NHI(n) = freduce*NHI(n)
	  CALL model(n,n,imin,imax,1)	! determines workflux
	  if (adjustb(n,imin,imax,ilow,ihi,ipivot).eq.-1) GOTO 20
 10	continue
 20	continue
c    	write(6,*) 'adjusted NHI:',NHI(n),bpar(n)

c  compute residual flux
	CALL model(1,n,1,ndata,1)
	if (bpar(n).lt.bminsat.AND.NHI(n).gt.NHImaxsat) then
     	  write(6,*) '*** unfittable line ***',imin,imax
	  do 40 i=imin,imax
	    resid(i) = 1.d0
 40	  continue
	else 
	  do 50 i=imin,imax
	    residold = resid(i)
	    resid(i) = 1.d0 + MAX(flux(i),0.d0) - workflux(i)
            if (resid(i).gt.1.d0) resid(i) = 1.d0
	    noise(i) = noise(i)*sqrt(MAX(resid(i),noise(i))/
     &	      (MAX(residold,noise(i))))
 50	  continue
	end if

	return 
	end

C=============================================================================

	logical function deconv(n,imin,imax,idown,iup,ipivot)

	include 'vpdefs.h'
	integer imin,imax,n,i,iup,idown,j,ipivot,ic
	double precision maxflux

c  identify if more than one saturated feature in detection region
	deconv = .FALSE.
	iup = 0
	idown = 0
	ic = centpix(n)
	scl = 1.d0
 5	continue
c	write(6,*) scl,ic,noise(ic)
	do 10 i=imin,imax-1
	  if (minflux(i).gt.scl*fsigma*noise(ic).AND.minflux(i+1).lt.
     &	    scl*fsigma*noise(ic)) then
	    if ((i-idown).gt.nwidth.AND.idown.gt.0) goto 20
	    if (idown.eq.0) idown = i
	  end if
	  if (minflux(i).lt.scl*fsigma*noise(ic).AND.minflux(i+1).gt.
     &	    scl*fsigma*noise(ic)) then
	    iup = i
	  end if
 10	continue
c  found only one up crossing and one down crossing
	if (iup-idown.lt.nwidth.AND.scl.lt.avenoise/noise(ic)) then
	  scl = scl + 1.d0
	  goto 5
	end if
	centpix(n) = (idown+iup)/2
	return
c  found more than one crossing: deconvolve
c  idown is 1st down crossing, iup is 1st up crossing, i is 2nd down crossing
 20	continue
	if (iup-idown.lt.nwidth.OR.i-iup.lt.nwidth
     &	  .AND.scl.lt.avenoise/noise(ic)) then
	  scl = scl + 1.d0
	  goto 5
	end if
	deconv = .TRUE.
c	write(6,*) 'imin, imax was',idown,iup,i
	maxflux = minflux(iup)
	imax = iup
	do 30 j=iup+1,i
	  if (minflux(j).gt.maxflux) then
	    maxflux = minflux(j)
	    imax = j
	  end if
 30	continue
c  new guess at central pixel value
	centpix(n) = (idown+iup)/2
	ipivot = imax
	if (imax-centpix(n).gt.centpix(n)-imin) ipivot = imin
	write(6,*) 'deconv: new imin,imax,centpix =', imin,imax,centpix(n)

	return
	end

C=============================================================================

	subroutine fitlim(n,imin,imax,ilow,ihi,ipivot)

	include 'vpdefs.h'
	integer imin,imax,n,i,ihi,ilow,ipivot,jhi,jlow

c	write(6,*) 'fitlim:',centpix(n),vel(centpix(n)),dfdv(centpix(n))
	if (ABS(dfdv(centpix(n))).gt.dfdvsat) then
     	  write(6,*) 'deriv at center too large',dfdv(centpix(n))
	  dfdvsat = 1.5d0*ABS(dfdv(centpix(n)))
	end if
	i = centpix(n)
	do 10 while (ABS(dfdv(i)).le.dfdvsat.AND.i.ge.imin+npixcusp)
	  ilow = i
	  i = i-1
 10	continue

	i = centpix(n)
	do 20 while (ABS(dfdv(i)).le.dfdvsat.AND.i.le.imax-npixcusp)
	  ihi = i
	  i = i+1
 20	continue

	centpix(n) = (ilow+ihi)/2
	ipivot = imax
	if (imax-centpix(n).gt.centpix(n)-imin) ipivot = imin
c	write(6,*) 'new',centpix(n),ilow,ihi

	return
	end

C=============================================================================

        integer function decreaseN(n,imin,imax,ilow,ihi)

        include 'vpdefs.h'
        integer imin,imax,n,i,ilow,ihi

c  decrease NHI if model is below minimum anywhere in cusp region
	decreaseN = 0
	do 10 i=ilow-npixcusp,ilow-1
c	write(6,*) i,workflux(i),minflux(i),decreaseN
	  if (workflux(i).lt.minflux(i)-fsigma*noise(i)) then
	    decreaseN = decreaseN+1
	    goto 15
	  end if
 10	continue
 15	continue
	do 20 i=ihi+1,ihi+npixcusp
	  if (workflux(i).lt.minflux(i)-fsigma*noise(i)) then
	    decreaseN = decreaseN+2
	    goto 25
	  end if
 20	continue
 25	continue

	do 30 i=ilow,ihi
	  if (workflux(i).gt.minflux(i)+fsigma*noise(i)) decreaseN=0
 30	continue

	return
	end

C=============================================================================

	integer function adjustb(n,imin,imax,ilow,ihi,ipivot)

	include 'vpdefs.h'
	integer imin,imax,n,i,ihi,ilow,ipivot
	double precision fminpiv

	adjustb = 0
	fminpiv = minflux(ipivot)-fsigma*noise(ipivot)
	if (workflux(ipivot).gt.fminpiv) then
	  if (ipivot.ge.ihi+npixcusp) ipivot = ipivot-1
	  if (ipivot.le.ilow-npixcusp) ipivot = ipivot+1
	end if
	fminpiv = minflux(ipivot)-fsigma*noise(ipivot)
c	write(6,*) 'b',ipivot,bpar(n),NHI(n)
	if (workflux(ipivot).lt.fminpiv) then
	  do 10 while (workflux(ipivot).lt.fminpiv)
	    if (bpar(n).le.bminsat) return
	    bpar(n)=freduce*bpar(n)
c	write(6,*) 'b',bpar(n),workflux(ipivot),fminpiv
	    CALL model(n,n,imin,imax,1)
 10	  continue
	else
	  do 20 while (workflux(ipivot).gt.fminpiv
     &	    .AND.workflux(ilow).lt.minflux(ilow)+fsigma*noise(ilow)
     &	    .AND.workflux(ihi).lt.minflux(ihi)+fsigma*noise(ihi))
	    if (bpar(n).ge.bparmax) return
	    bpar(n)=bpar(n)/freduce
	    CALL model(n,n,imin,imax,1)
 20	  continue
	end if

	if (ipivot.ge.ihi) then
	do 35 ipiv = ipivot,imax,ipivint
	  fminpiv = minflux(ipiv)-fsigma*noise(ipiv)
	  do 30 while (workflux(ipiv).lt.fminpiv)
	    bpar(n)=bpar(n)*freduce
	    CALL model(n,n,imin,imax,1)
 30	  continue
 35	continue
	else
	do 45 ipiv = ipivot,imin,ipivint
	  fminpiv = minflux(ipiv)-fsigma*noise(ipiv)
	  do 40 while (workflux(ipiv).lt.fminpiv)
	    bpar(n)=bpar(n)*freduce
	    CALL model(n,n,imin,imax,1)
 40	  continue
 45	continue
	end if

	return
	end

