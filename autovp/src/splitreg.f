
	integer function splitreg(ireg)
c
c  if region reg has too many lines, split at a location somewhere
c  in the middle where the flux is highest.
c  truncates the region ireg, and adds a new region to end of list.
c  
	include 'vpdefs.h'

	integer ireg,nl,j,l,jsplit,ilow,ihi
	double precision vl(maxregions)
	double precision fmax

 1	continue
	nl = 0
        vmin = vel(region(ireg,1))
        vmax = vel(region(ireg,2))
        do 20 j=1,nlines
          if (vline(j).ge.vmin.AND.vline(j).le.vmax) nl = nl+1
 20	continue
	splitreg = nl
	if (nl.le.nsplitmin) return
c
c  find region within which split point lies: betn lines nlinemin & nl-nlinemin
	nl = 0
        do 25 j=1,nlines
          if (vline(j).ge.vmin.AND.vline(j).le.vmax) then
	    nl = nl+1
	    vl(nl) = (vline(j)-vel(1))*ndata/(vel(ndata)-vel(1))
	  end if
 25	continue
	CALL PIKSRT(nl,vl)
	ilow = 1+INT(vl(nlinemin))
	ihi = INT(vl(nl-nlinemin))
c
c  find maximum flux pixel where deriv is approx 0
	CALL polyfit(region(ireg,1),region(ireg,2))
	fmax = 0.d0
	jsplit = ilow-1
	do 30 j=ilow,ihi
	  if (dfdv(j)*dfdv(j+1).lt.0.d0.AND.
     &	    dfdv(j-nint/2)*dfdv(j+nint/2+1).lt.0.d0) then
	    if (flux(j).gt.fmax) then
	      fmax = flux(j)
	      jsplit = j
	    end if
	  end if
 30	continue
c
c  if that didn't work, split at max flux
	if (jsplit.eq.ilow-1.OR.flux(jsplit).lt.0.1) then
	  fmax = 0.d0
	  do 40 j=ilow,ihi
	    if (flux(j).gt.fmax) then
	      fmax = flux(j)
	      jsplit = j
	    end if
 40	  continue
	end if
c
c  split region into two at pixel jsplit
	write(6,*) 'SPLIT REGION at',jsplit,vel(jsplit)
	nregions = nregions+1
	region(nregions,1) = jsplit+1
	region(nregions,2) = region(ireg,2)
	region(ireg,2) = jsplit
c  try further splitting if necessary
	goto 1

	end

