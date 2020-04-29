      subroutine combine
c
c  combine small very nearby lines
c  combine if both lines within bfact*b of each other, and have N<NHImaxns
c
	include 	'vpdefs.h'
	double precision dv,dvmin
	integer i,j,k,ncomb

	ncomb = 1
	if( nlines.le.1 ) return 	! nothing to combine
c  keep combining until none left to combine
 	do 5 while (ncomb.ne.0)
	ncomb = 0
	do 10 i=1,nlines
c  find closest line to i, call it line j
	  dvmin = 10000000.
	  do 20 k=1,nlines
	    dv = ABS(vel(centpix(i))-vel(centpix(k)))
	    if (dv.lt.dvmin.AND.i.ne.k) then
	      dvmin = dv
	      j = k
	    endif
 20	  continue
c  if i and j are close and small enough, combine into line i
	  dv = ABS(vel(centpix(i))-vel(centpix(j)))
	  if (bfact*bpar(i).gt.dv.AND.bfact*bpar(j).gt.dv.AND.
     &	    NHI(i).lt.NHImaxns.AND.NHI(j).lt.NHImaxns.AND.i.ne.j) then
	    NHI(i) = NHI(i) + NHI(j)			! add column densities
	    bpar(i) = sqrt(bpar(i)**2+bpar(j)**2) + dv  ! expand width by dv
	    centpix(i) = INT((NHI(i)*centpix(i)+
     &		NHI(j)*centpix(j))/(NHI(i)+NHI(j)))
c  set values of line j so that it will not be combined elsewhere
	    NHI(j) = NHImaxns
	    bpar(j) = 0.d0
	    ncomb = ncomb + 1
	write(6,*) 'combining',ncomb,i,j,dv
	  end if
 10	continue
c	do 23 i=1,nlines
c	write(6,*) i,NHI(i),bpar(i)
c 23	continue
c  remove lines which have been combined into others
	i = 0
	do 40 k=1,nlines
	  if (bpar(k).ne.0.) then
	    i = i+1
	    NHI(i) = NHI(k)
	    bpar(i) = bpar(k)
	    centpix(i) = centpix(k)
	  end if
 40	continue
	nlines = i
c
 5	continue

	return 
	end
