	subroutine doublet
c
c  Initial (crude & conservative) check for doublets
c
	include 	'vpdefs.h'
	include 	'const.dek'
	integer i,j,k,idata,npixrej
	double precision dwave,wave2,z,flux2,tau1,tau2
c
	write(6,*) 'Looking for doublets...'
	dwave = (waveln(ndata)-waveln(1))/ndata
	npixrej = INT(bminsat/(vel(2)-vel(1)))
	do 5 i=1,nlines
	  lineinfo(i) = 1
 5	continue
c  loop until all detection regions processed
	do 10 i=1,nlines
	  z = (1.d0+redz)*(1.d0+vel(centpix(i))/ckms) - 1.d0
	  wave2 = (1.d0+z)*lambda1
C	  write(6,*) redz,z,vel(centpix(i))
	  idata = (wave2-waveln(1))/dwave
	  do 20 while(waveln(idata).LT.wave2.AND.waveln(idata+1).LT.wave2
     &	    .AND.idata.le.ndata)
	    idata = idata + 1
 20	  continue
	  do 30 while(waveln(idata).GT.wave2.AND.waveln(idata+1).GT.wave2
     &	    .AND.idata.ge.1)
	    idata = idata - 1
 30	  continue
	  if( idata.le.1.OR.idata.ge.ndata ) then
	    lineinfo(i) = 0
	    write(6,*) 'doublet out of range: ',i,idata
	    goto 10
	  endif
	  tau1 = -log(fsigma*noise(centpix(i)))
	  if( flux(centpix(i)).GT.fsigma*noise(centpix(i)) ) 
     &	    tau1 = -log(flux(centpix(i)))
	  flux2 = MIN(flux(idata),flux(idata+1))
	  tau2 = -log(fsigma*noise(idata))
	  if( flux2.GT.fsigma*noise(centpix(i)) ) 
     &	    tau2 = -log(flux2)
	  if( tau2.LT.fdoublet*tau1 ) lineinfo(i) = 0
	  write(6,'(2i6,5f11.4)') i,lineinfo(i),
     &	    waveln(centpix(i)),wave2,
     &	    tau1,tau2
	  if( lineinfo(i).eq.1 ) then
	    do 40 j=1,nlines
	      do 40 k=-npixrej,npixrej
	        if( idata+k.eq.centpix(j) ) then
		  lineinfo(j) = 0
		end if
 40	    continue
	  end if
 10	continue
	do 50 i=1,nlines
C	  if( lineinfo(i).EQ.1 ) write(6,*) i,centpix(i),
C     &	waveln(centpix(i)),flux(centpix(i))
	  if( lineinfo(i).EQ.0 ) bpar(i) = 0.d0
 50	continue

	return
	end

	subroutine doublet2
c
c  Second more stringent check for doublets after line minimization
c
	include 	'vpdefs.h'
	include 	'const.dek'
	integer i,j,k,idat1,idat2,npixrej
	double precision dwave,wave1,wave2,z,flux1,flux2,tau1,tau2

	write(6,*) 'Looking for doublets...'
	dwave = (waveln(ndata)-waveln(1))/ndata
	npixrej = INT(bminsat/(vel(2)-vel(1)))
	do 5 i=1,nlines
	  lineinfo(i) = 1
 5	continue
c  loop until all detection regions processed
	do 10 i=1,nlines
	  z = (1.d0+redz)*(1.d0+vline(i)/ckms) - 1.d0
	  wave1 = (1.d0+z)*lambda0
	  wave2 = (1.d0+z)*lambda1
	  idat1 = (wave1-waveln(1))/dwave
	  idat2 = (wave2-waveln(1))/dwave
C	  write(6,*) redz,z,vline(i)
	  do 18 while(waveln(idat1).LT.wave1.AND.waveln(idat1+1)
     & 		.LT.wave1)
	    idat1 = idat1 + 1
 18	  continue
	  do 28 while(waveln(idat1).GT.wave1.AND.waveln(idat1+1)
     &		.GT.wave1)
	    idat1 = idat1 - 1
 28	  continue
	  centpix(i) = idat1
	  do 20 while(waveln(idat2).LT.wave2.AND.waveln(idat2+1).LT.
     &		wave2)
	    idat2 = idat2 + 1
 20	  continue
	  do 30 while(waveln(idat2).GT.wave2.AND.waveln(idat2+1).GT.
     &		wave2)
	    idat2 = idat2 - 1
 30	  continue
	  if( idat2.le.1.OR.idat2.ge.ndata ) then
	    lineinfo(i) = 0
	    write(6,'(2i6,4f11.4,a15)') i,lineinfo(i),
     &		waveln(centpix(i)),wave2,tau1,tau2,'*out of range'
	    goto 10
	  endif
	  flux1 = MIN(workflux(idat1),workflux(idat1+1))
	  tau1 = -log(fsigma*noise(idat1))
	  if( flux1.GT.fsigma*noise(idat1) ) tau1=-log(flux1)
	  flux2 = MIN(workflux(idat2),workflux(idat2+1))
	  flux2 = flux2 - noise(idat2)
	  tau2 = -log(fsigma*noise(idat2))
	  if( flux2.GT.fsigma*noise(idat2) ) tau2=-log(flux2)
C	  if( flux1.gt.(1.-0.5*N_sigma*noise(idat1))) then
C	    lineinfo(i) = 0
C	    write(6,'(2i6,4f11.4,a10)') i,lineinfo(i),
C     &		waveln(centpix(i)),wave2,tau1,tau2,'*too small'
C	    goto 10
C	  endif
	  if( tau2.LT.fdoublet*tau1 ) lineinfo(i) = 0
	  write(6,'(2i6,2f10.3,3f10.5,f10.2)') i,lineinfo(i),
     &	    waveln(centpix(i)),wave2,
     &	    tau1,tau2,NHI(i),vline(i)
	  if( lineinfo(i).eq.1 ) then
	    do 40 j=1,nlines
	      do 40 k=-npixrej,npixrej
	        if(idat2+k.eq.centpix(j).AND.NHI(i).gt.NHI(j)) then
		  lineinfo(j) = 0
	    write(6,'(a11,2i6,4f11.4)') 'eliminating:',j,i,
     &		waveln(centpix(i)),waveln(centpix(j)),tau1,tau2
		end if
 40	    continue
	  end if
 10	continue
	do 50 i=1,nlines
C	  if( lineinfo(i).EQ.1 ) write(6,*) i,centpix(i),
C     &	waveln(centpix(i)),flux(centpix(i))
	  if( lineinfo(i).EQ.0 ) bpar(i) = 0.d0
 50	continue

c  remove lines which have been combined into others
        i = 0
        do 60 k=1,nlines
          if (bpar(k).ne.0.) then
            i = i+1
            NHI(i) = NHI(k)
            bpar(i) = bpar(k)
            vline(i) = vline(k)
          end if
 60     continue
        nlines = i

	return
	end


