
	integer function tossline(parm,parerr,np,index)
c
c  Decides if any lines are worth trying to get rid of.
c  First, if any line is very small, toss it immediately.
c  Otherwise, find the weakest line (i.e. has greatest error in NHI and b)
c  and try tossing that line (if errors > parameter values).
c  Tosses line by switching tossed line to the end of the line
c  list and reducing the number of lines by 1.  
c  If it is deemed that the line CANNOT be tossed without substantial
c  increase in chisq, we can then recover line info easily.
c
	include 'vpdefs.h'
	double precision parm(*),parerr(*)
	double precision weakest,colmax,weakness
	integer np,iweak
	integer index(*)

	tossline = 0
c  don't toss if it's the last line
	if (np.le.3) return
c
c  if any line is VERY small, toss it
	do 20 i=1,np,3
	  iline = index(i)
	  if (lineinfo(iline).eq.-2) goto 20
	  if (NHI(iline).le.2.*NHImin) then
	    CALL switchlines(iline,nlines)
	    nlines = nlines-1
	    write(6,*) 'unconditionally tossing small line',iline
	    tossline = -iline		! negative signals unconditional toss
	    return
	  end if
 20	continue
	do 22 i=1,np,3
	  iline = index(i)
	  if (lineinfo(iline).eq.-2) goto 22
	  if (bpar(iline).le.2.*bparmin) then
	    CALL switchlines(iline,nlines)
	    nlines = nlines-1
	    write(6,*) 'unconditionally tossing small line',iline
	    tossline = -iline		! negative signals unconditional toss
	    return
	  end if
 22	continue
c
c  if a line is predetermined to be tossed, skip over criteria and just toss
	do 5 i=1,nlines
	  if (lineinfo(i).eq.2) then
	    iweak = i
	    lineinfo(i) = 0
	    goto 80
	  end if
 5	continue
c  determine weakest line using weakness criterion 
	iweak = 0
	do 10 i=1,np,3
	  iline = index(i)
	  if (lineinfo(iline).lt.0) goto 10
	  weakness = sqrt((parerr(i+1)/parm(i+1))**2+
     &	    (parerr(i+2)/parm(i+2))**2)
	  if (i.eq.1) weakest = weakness
	  if (weakest.le.weakness) then
	    weakest = weakness
	    iweak = iline
	  end if
 10	continue
	if (weakest.le.1.d0) return

c  find any lines that have NHI factor of 10 below largest; try tossing those
	if (weakest.lt.0.d0) then	! NOT USED (weakest never lt 0)
	  colmax = 0.d0
	  do 30 i=1,np,3
	    colmax = MAX(colmax,parm(i+1))
 30	  continue
	  do 40 i=1,np,3
	    if (lineinfo(index(i)).lt.0) goto 40
	    weakness = colmax/parm(i+1)
	    if (weakness.gt.20.d0.AND.weakest.le.weakness) then
	      iweak = index(i)
	      weakest = weakness
	    end if
 40	  continue
	  if (weakest.le.1.d0) return
	end if

c  toss weakest line by filling last line into slot of weakest
c  retain info of weakest line in last line's slot
 80	continue
	CALL switchlines(iweak,nlines)
	tossline = iweak
	if( tossline.ne.0 ) nlines = nlines-1

	return
	end

	logical function acctoss(chisq,ochisqtoss,parm,np,nd)

	include 'vpdefs.h'
	double precision parm(*)
	double precision chisq,ochisqtoss
	integer np,nd

	acctoss = .FALSE.
C	if (chisq.lt.ochisqtoss.OR.ochisqtoss.eq.0.d0) acctoss = .TRUE.
	if (chisq-ochisqtoss.lt.faccept*chisq.OR.ochisqtoss.eq.0.d0
     &    .OR.(chisq/nd.lt.chisqgood.AND.chisq.lt.2.*ochisqtoss))
     &	  acctoss = .TRUE.

	return
	end

	subroutine switchlines(i,j)

	include 'vpdefs.h'
	integer i,j

	CALL switch(NHI(i),NHI(j))
	CALL switch(bpar(i),bpar(j))
	CALL switch(vline(i),vline(j))
	CALL switch(dNHI(i),dNHI(j))
	CALL switch(dbpar(i),dbpar(j))
	CALL switch(dvline(i),dvline(j))
	itemp = lineinfo(i)
	lineinfo(i) = lineinfo(j)
	lineinfo(j) = itemp

	return
	end

	subroutine switch(x,y)

	double precision x,y,temp

	temp = x
	x = y
	y = temp

	return
	end
