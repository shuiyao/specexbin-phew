      program minfit
c
c
	include 'vpdefs.h'
	include 'const.dek'

        double precision fluxfit(nmax),error(nmax)
	double precision parm(maxparm),parerr(maxparm)
	double precision parmold(maxparm)
	double precision covar(maxparm,maxparm),alpha(maxparm,maxparm)
	double precision alamda,chisq,oldchisq
	integer listparm(maxparm),index(maxparm)
	integer i,j,k,kk,np,npv
	integer tossline,itoss
	double precision ochisqtoss,ochisqadd
	double precision tau,dew
        double precision vhi,vlow
	integer splitreg
	logical acctoss,equivwidth
	character*80 infile,outfile
	integer ew_finder
c
c  get the input data: cln file and first-guess lines in PROFIT format
        call getarg(1,infile)
	CALL getlines(infile)
	CALL readdata
c
c  find detection regions
	i = ew_finder(1,ndata,region,nregions,resid)
	write(6,*) '# of regions, lines:',nregions,nlines
c
c  fit Voigt profiles in each detected region
	i = 0
	do 40 while (i.lt.nregions)
	  i = i+1
c  if too many lines, split region intelligently; new region added at end
c  if no lines, ignore region
	  vmin = vel(region(i,1))
	  vmax = vel(region(i,2))
	  if( splitreg(i).eq.0 ) then
	    write(6,*) ' -----------------------------------------'
	    write(6,*) ' SKIPPING REGION ',i,vmin,vmax
	    goto 40
	  end if
	  vmin = vel(region(i,1))
	  vmax = vel(region(i,2))
	  nd = region(i,2)-region(i,1)+1
c  load spectral data for detection region i
	  do 50 j=region(i,1),region(i,2)
	    k = j - region(i,1) + 1
	    velfit(k) = vel(j)
	    fluxfit(k) = 1.d0 - workflux(j) + flux(j)
	    error(k) = noise(j)
 50	  continue
	  write(6,*) ' -----------------------------------------'
	  write(6,*) ' DOING REGION',i,vmin,vmax
c  initialize line tossing/adding parameters
	  do 52 j=1,nlines
	    lineinfo(j) = 1	! active line, which hasn't been tossed
 52	  continue
	  ochisqtoss = 0.d0
	  itoss = 0
	  ochisqadd = 1.d6*chisqbad*ndata/faccept  ! start huge
c  loop for tossing lines
 55	  continue
c  load line parameters
	  CALL initparm(parm,parerr,index,listparm,np,npv)
	  if (np.eq.0) goto 40

c  BEGIN CHI-SQ ITERATION LOOP
	  chisq = 0.d0
	  alamda = -1.d0
     	  write(6,*) 'starting chisq = ',chi2(velfit,fluxfit,error,
     &        nd,parm,np)
	  do 70 niter = 1,maxiter
	    oldchisq = chisq
c  do one interation of Marquardt minimization
	    CALL MRQMIN(velfit,fluxfit,error,nd,parm,np,listparm,
     &		npv,covar,alpha,maxparm,chisq,alamda)
c	do kk=1,nd
c	write(6,'(3f10.4)') velfit(kk),fluxfit(kk),error(kk)
c	end do
c  output if iteration was successful (reduced chisq)
	    if (chisq.lt.oldchisq.OR.alamda.ge.alamdamax) then
     	      write(6,*) 'successful iter:',niter,
     &	        chisq,oldchisq
	    end if
c  check covariance matrix to make sure no lines are irrelevant to fit
	    do 77 j=1,npv
	      nzero = 0
	      do 78 k=1,npv
	        if (covar(j,k).eq.0.0) nzero = nzero+1
 78	      continue
c	write(6,*) '**** covariance matrix ****'
c	do kk=1,npv
c	write(6,'(15g10.1)') (covar(kk,k),k=1,npv)
c	end do
c  if an entire row of covar matrix is 0, freeze that parameter and restart
	      if (nzero.eq.npv) then
		CALL fixparm(j,parm,parerr,np,npv,listparm,index)
		alamda = -1.d0		! signal to restart minimization
		goto 80
	      end if
 77	    continue
	    if (ABS(chisq-oldchisq)/chisq.lt.chisqtol.AND.
     &	      chisq.lt.oldchisq) goto 80
	    if (alamda.ge.alamdamax) goto 80
 70	  continue 
	  write(6,*) 'uh-oh...overiterated!!!',chisq
 80	  continue

c  DONE CHI-SQ ITERATION LOOP -- compute covariance matrix
	  alamda = 0.d0		! signal to compute final covar matrix
	  CALL MRQMIN(velfit,fluxfit,error,nd,parm,np,listparm,
     &	    npv,covar,alpha,maxparm,chisq,alamda)
c  output fit values
	  write(6,'(a,f10.3,i8,3f11.3)') 'FIT RESULTS:',
     &	    chi2(velfit,fluxfit,error,nd,parm,np)/nd,niter,vmin,vmax,
     &	    lambda0*(1.+redz)*(1.+0.5*(vmin+vmax)/ckms)
	  write(6,*) ' i     NHI         dNHI            v          dv        
     &    b          db'
	  do 85 j=1,np,3
c  approximate errors as sqrt of diag of covariance matrix
	    parerr(j) = sqrt(ABS(covar(j,j)))
	    parerr(j+1) = sqrt(ABS(covar(j+1,j+1)))
	    parerr(j+2) = sqrt(ABS(covar(j+2,j+2)))
	    write(6,'(i3,3(2f12.4,1x))') index(j),parm(j+1),parerr(j+1),
     &		parm(j),parerr(j),
     &		parm(j+2),parerr(j+2)
 85	  continue

c  TRY TOSSING UNNECCESARY LINE AND REFITTING
c  if a line was tossed, accept only if new fit comparable to old fit
	  if (acctoss(chisq,ochisqtoss,parm,np,nd)) then
c  load new values into line parameters
	    CALL loadline(parm,parerr,np,index) ! accept latest fit
	    if (itoss.ne.0) lineinfo(ABS(itoss)) = 0	! tag tossed line
c  see if we can toss another line
	    itoss = tossline(parm,parerr,np,index)
	    if (itoss.ne.0) then
	      write(6,*) 'tossing line ',ABS(itoss),' and refitting'
	      ochisqtoss = chisq
	      if (itoss.lt.0) ochisqtoss = 0.d0
	      goto 55
	    end if
c  if line NOT tossable, add line back in and try tossing another
	  else 
 86	    continue
	    nlines = nlines+1		! add line back in (at end)
	    write(6,*) 'Toss NOT accepted ',itoss,chisq,nlines
	    lineinfo(nlines) = -1	! tag line as not tossable
	    CALL initparm(parm,parerr,index,listparm,np,npv)  ! reload all lines
	    chisq = ochisqtoss
	    itoss = tossline(parm,parerr,np,index)
	    if (itoss.ne.0) then
	      write(6,*) 'tossing line ',itoss,' and refitting'
	      goto 55
	    end if
	    CALL loadline(parm,parerr,np,index) ! accept previous fit
	  end if
c  IF NOT GOOD FIT, TRY MINIMIZING ONE PARAMETER AT A TIME
	  if (chisq/nd.gt.chisqbad) then
	    write(6,*) 'One-parameter minimization',chisq
	    CALL initparm(parm,parerr,index,listparm,np,npv) ! reload all lines
	    ochisqtoss = chisq
	    chisqold = 2.*chisq
c  do one-param minimization as long as it keeps reducing chisq
	    do 94 while (chisq-chisqold.lt.-chisqtol*chisq)
	      chisqold = chisq
     	      CALL oneparmin(velfit,fluxfit,error,nd,parm,np,
     &	        listparm,npv,chisq)
    	      CALL nparmin(velfit,fluxfit,error,nd,parm,np,
     &	        listparm,npv,chisq)
 94	    continue
c  if one param minimization helped, then go back and try Marquardt again
	    if (chisq-ochisqtoss.lt.-chisqtol*chisq) then
	      CALL loadline(parm,parerr,np,index)  ! accept new parms
	      write(6,*) 'Reduced chisq to',chisq
	      ochisqtoss = 0.d0
	      goto 55
	    end if
	  end if
c  IF STILL NOT GOOD FIT, TRY ADDING LINE
	  if (chisq/nd.gt.chisqbad) then
	    if (chisq.lt.(1.d0-faccept)*ochisqadd) then
	      ochisqadd = chisq
	      ochisqtoss = 0.d0
	      do 91 j=1,np
		parmold(j) = parm(j)
		if (lineinfo(index(j)).eq.-2) lineinfo(index(j)) = -1
 91	      continue
	      npold = np
	      CALL addline(i,velfit,fluxfit,error,nd,parm,np)
	      goto 55
	    else		! set lines back to old fit
	      do 92 j=1,np
		parm(j) = parmold(j)
		if (lineinfo(index(j)).eq.-2) then
		  lineinfo(index(j)) = 2
		  itoss = tossline(parm,parerr,np,index)
		end if
 92	      continue
	      write(6,*) 'trial line NOT accepted',itoss
	      chisq = ochisqadd
	      np = npold
	      CALL loadline(parm,parerr,np,index)
	    end if
	    CALL initparm(parm,parerr,index,listparm,np,npv) ! reload all lines
	  end if

c  THE END: computes working model (workflux) with new lines inserted
	  if (chisq/nd.gt.chisqbad) write(6,'(a,a12,3f12.4)') 
     &	    'CHECK FIT ',infile,chisq/nd,vmin,vmax
	  do 89 j=1,np,3
	    CALL model(index(j),index(j),1,ndata,2) 
 89	  continue
 40	continue
c
c  output all lines in PROFIT format
	i = 0
	do 95 while (infile(i:i).ne.'.')
	  i = i+1
 95	continue
        call fappend(infile(1:i),'vpm',outfile)
	write(6,*) 'DONE FIT -- OUTPUTTING TO ',outfile(1:20)

c  compute mean optical depth
	do 96 j=1,ndata
	  workflux(j) = 1.d0
 96	continue
	CALL model(1,nlines,1,ndata,2) 
	tau = 0.d0
	k=0
	do 87 j=1,ndata
C	  if( j.gt.240.AND.j.lt.320 ) write(6,*) 'TAU! ',j,workflux(j-1),
C     &	vhelio,-log(vhelio)
C	  if( workflux(j).ne.0.d0 ) vhelio = workflux(j)
	  if( workflux(j).gt.0.1 ) then
	    k = k+1
	    tau = tau - log(workflux(j))
	  end if
 87	continue
	tau = tau/k

c  Check for doublet pairs
        if( fdoublet.ne.0.d0.AND.lambda1.ne.0.d0 ) then
	  CALL doublet2
	endif

c  Calculate equiv widths of lines
	i = ew_finder(1,ndata,region,nregions,flux)
	if( fitflag.eq.2 ) equivwidth = .TRUE.
	if( equivwidth ) then
	  do 202 i=1,nlines
	    do 204 j=1,ndata
	      workflux(j) = 1.d0
 204	    continue
	    CALL model(i,i,1,ndata,2) 
	    ewline(i) = 0.d0
C	    if( bpar(i).lt.bparmax ) then
	    do 206 j=2,ndata
	      dew = (waveln(j)-waveln(j-1))*(1.-workflux(j))
     &	/(waveln(j)/lambda0)		! convert to rest equiv width
	      ewline(i)=ewline(i)+dew
 206	    continue
C	    else
C	    do 208 j=1,nregions
C	      if( vline(i).ge.vel(region(j,1)).AND.
C     &		  vline(i).lt.vel(region(j,2)) ) ireg = j
C 208	    continue
C	    do 207 j=2,ndata
C	      if(j.ge.region(ireg,1).AND.j.le.region(ireg,2)) 
C     &	      ewline(i)=ewline(i)+
C     &		(waveln(j)-waveln(j-1))*MAX(0.d0,1.-flux(j)-sigma(j))
C 207	    continue
C	    end if
 202	  continue
	end if
	
c  Output final model spectrum
	errmax = 999999.d0
	do 98 j=1,ndata
	  workflux(j) = 1.d0
 98	continue
	do j=1,nlines
	    NHI(j) = MIN(errmax,NHI(j))
	    bpar(j) = 0.01*(INT(bpar(j)*100.+0.5)/1)
	    if( ewline(j).gt.0.45 .AND. ewline(j).lt.0.75 ) 
     &		CALL model(j,j,1,ndata,2) 
	end do
	do j=1,ndata
	  resid(j) = workflux(j)
	  workflux(j) = 1.d0
	end do
	do j=1,nlines
	    NHI(j) = MIN(errmax,NHI(j))
	    bpar(j) = 0.01*(INT(bpar(j)*100.+0.5)/1)
	    if( ewline(j).lt.10000.) CALL model(j,j,1,ndata,2) 
	end do
        call fappend(ion_name,'mod',sig_file)
        open(unit=2,file=sig_file,status='unknown')
        do 88 j=1,ndata
          write(2,'(f12.4,f16.4,3f12.6)') waveln(j),vel(j),workflux(j),
     &       noise(j),resid(j)
 88     continue
        close(2)

	nelems = 2
	vhelio = 0.d0
	do 97 j=2,ndata
	   if( flux(j).ne.1.0 ) wavehi = waveln(j-1)
  97	continue
        open(unit=2,file=outfile,status='unknown')
c	write(2,'(1x,i5,4f14.6)') nelems,redz,vhelio,
c     &	  (wavehi-waveln(1))/lambda0,tau
c	write(2,'(1x,a10,i5)') ion_name,nlines
        do 90 j=1,nlines
            i = 0
            do 901 while(vel(i)<vline(j))
c                write(6,*) vline(j),vel(i),j,i
                vlow = vel(i)
                i = i+1 
  901       continue
            vhi = vel(i)
            write(6,*) vlow, vhi, vel(i),vline(j)
            write(6,*) rho(i-1),rho(i),met(i-1),met(i)
            rhoout(j) = ((vel(i)-vline(j))*rho(i-1)+ 
     &          (vline(j)-vel(i-1))*rho(i))/(vel(i)-vel(i-1))
            tempout(j) = ((vel(i)-vline(j))*temp(i-1)+
     &          (vline(j)-vel(i-1))*temp(i))/(vel(i)-vel(i-1))
            metout(j) = ((vel(i)-vline(j))*met(i-1)+
     &          (vline(j)-vel(i-1))*met(i))/(vel(i)-vel(i-1))
            write(6,*) rhoout(j),metout(j) 
	  if( equivwidth ) 
     &    write(2,'(i4,1x,f11.2,1x,f11.1,1x,f10.2,1x,
     &     f10.2,1x,f10.2,1x,f10.2,1x,f9.3,1x,f10.3,1x,f10.3,1x,
     &      f10.3,1x,f12.5)')
     &      j,MIN(errmax,NHI(j)),vline(j),bpar(j),MIN(errmax,dNHI(j)),
     &      MIN(errmax,dvline(j)),MIN(errmax,dbpar(j))
     &      ,ewline(j),lambda0*exp(vline(j)/2.99792458e5),
     &      rhoout(j),tempout(j),metout(j)
	  if( .NOT.equivwidth ) 
     &    write(2,'(i4,1x,f12.3,1x,f11.3,1x,f10.3,1x,f10.3,1x,
     &     f10.3,1x,f10.3,1x,f10.3,1x,f10.3,1x,f12.5)')
     &      j,MIN(errmax,NHI(j)),vline(j),bpar(j),MIN(errmax,dNHI(j)),
     &      MIN(errmax,dvline(j)),MIN(errmax,dbpar(j)),rhoout(j),
     &      tempout(j),metout(j)
 90     continue
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
