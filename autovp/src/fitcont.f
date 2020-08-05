      program fitcont
c
c  fits a continuum level by recursive iteration fit of 4th-order
c  spline polynomial.
c
	implicit double precision (a-h,o-z)
      integer            nmax,fitorder
      parameter          (nmax = 5000, fitorder=3)
      integer            i,j,k,ndata,niter
      double precision   wave(nmax),vel(nmax),
     &                   signorm(nmax),splinefit(nmax)
      logical            fitdata(nmax)
      double precision   velfit(nmax),sigfit(nmax),fluxfit(nmax)
      double precision   coeff(fitorder)
      double precision	 alpha(fitorder,fitorder),
     &			 covar(fitorder,fitorder)
      integer		 listcoeff(fitorder),nfit,nwrap
      double precision   lowlim,hilim,chisq,chisqold
      double precision	 sigmax,fluxmax,fmidmax,velmax,dvel,dwave
      character*80       ion_name,sig_file,cln_file
	double precision fluxnorm
	common/fitcom/fluxnorm(nmax),sigmax,ndata
c
c  get the ion info for this run
      call getarg(1,ion_name)
      call fappend(ion_name,'raw',sig_file)
      open(unit=1,file=sig_file,err=1000,status='old')
	write(6,*) 'FITTING CONTINUUM: ',sig_file(1:40)
c  read in data
      ndata = 0
      sigmax = 0.
	fluxmax = 0.
      do 5 i=1,nmax
       read(1,*,end=6) wave(i),vel(i),fluxnorm(i),signorm(i)
       sigmax = MAX(sigmax,signorm(i))
       fluxmax = MAX(fluxmax,fluxnorm(i))
	if (fluxmax.eq.fluxnorm(i)) ifmax = i
       fitdata(i) = .TRUE.
       ndata = ndata + 1
  5   continue
  6   close(unit=1)
      velmax = vel(ndata)
c  add wrap-around to reduce boundary effects
      nwrap = 0
      dvel = (vel(ndata)-vel(1))/ndata
      dwave = (wave(ndata)-wave(1))/ndata
c  shift data up by nwrap
      do 7 i=ndata,1,-1
	wave(i+nwrap) = wave(i)
	vel(i+nwrap) = vel(i)
	fluxnorm(i+nwrap) = fluxnorm(i)
	signorm(i+nwrap) = signorm(i)
 7    continue
c  add portions below and above, using periodic BC
      do 8 i=1,nwrap
	wave(i) = wave(nwrap+1) + (i-1-nwrap)*dwave
	vel(i) = vel(nwrap+1) + (i-1-nwrap)*dvel
	fluxnorm(i) = fluxnorm(i+ndata)
	signorm(i) = signorm(i+ndata)
 8    continue
      do 9 i=ndata+1,ndata+nwrap
	wave(i) = wave(ndata) + (i-ndata)*dwave
	vel(i) = vel(ndata) + (i-ndata)*dvel
	fluxnorm(i) = fluxnorm(i-ndata)
	signorm(i) = signorm(i-ndata)
 9    continue
      ndata = ndata+2*nwrap
c
c  initialize fit
      lowlim = 2.5
      hilim = 1000000.d0
      niter = 10
      chisq = 0.d0
c  constrain at endpoints and at max flux value
      do 30 i=3,fitorder
	coeff(i) = 0.d0
	listcoeff(i-2) = i
 30   continue
      coeff(1) = fluxnorm(1)-sigmax
c  determine coeff(2) and coeff(3) from maximum in middle of spectrum
	fmidmax = fluxnorm(ndata/5)
	vfmax = vel(ndata/5)
	do 12 i=ndata/5+1,4*ndata/5
	    if(fluxnorm(i).gt.fmidmax) then
		fmidmax = fluxnorm(i)
		ifmax = i
	    end if
 12	continue
	fmidmax = fmidmax - signorm(ifmax)
      vfmax = vel(ifmax)/velmax
      coeff(2) = (vfmax*vfmax*fluxnorm(ndata)-fmidmax-
     &	(vfmax**2-1)*coeff(1))/(vfmax**2-vfmax)
      coeff(3) = (vfmax*fluxnorm(ndata)-fmidmax-
     &	(vfmax-1)*coeff(1))/(vfmax-vfmax**2)
c  exclude anything 4*lowlim*sigmamax below first guess continuum curve
	do 16 i=1,ndata
	  splinefit(i) = 0.d0
	  do 17 k=1,fitorder
	    splinefit(i)=splinefit(i)+coeff(k)*(vel(i)/velmax)**(k-1)
 17	  continue
 16	continue
	do 32 i=1,ndata
	  fitdata(i) = .TRUE.
	  if (fluxnorm(i).lt.splinefit(i)-4.*lowlim*signorm(i))
     &	    fitdata(i) = .FALSE.
 32	continue
	alamda = -1.d0
	chisq = 0.d0
      do 33 j=1,niter
c  load appropriate data into arrays to be fitted
	nfit = 0
	do 34 i=1,ndata
	  if (fitdata(i)) then
	    nfit = nfit+1
	    velfit(nfit) = vel(i)/velmax
	    fluxfit(nfit) = fluxnorm(i)
	    sigfit(nfit) = sigmax
	  endif
 34	continue
c  fit polynomial to non-excluded data set
c	chisqold = chisq
c	write(6,'(2i7,6g12.4)') j,nfit,chisq,(coeff(k),k=1,fitorder)
        CALL MRQMIN(velfit,fluxfit,sigfit,nfit,coeff,fitorder,listcoeff,
     &    fitorder-2,covar,alpha,fitorder,chisq,alamda)
	coeff(2) = fluxnorm(ndata)-fluxnorm(1)
        do 45 i=3,fitorder
     	  coeff(2) = coeff(2) - coeff(i)
 45	continue
c  compute polynomial value at each pixel
	do 36 i=1,ndata
	  splinefit(i) = 0.d0
	  do 37 k=1,fitorder
	    splinefit(i)=splinefit(i)+coeff(k)*(vel(i)/velmax)**(k-1)
 37	  continue
 36	continue
c  if converged, exit loop
c	if (ABS(chisq-chisqold)/chisq.lt.0.001) GOTO 50
c  exclude points too far above/below
	do 38 i=1,ndata
	  fitdata(i) = .TRUE.
	  if (fluxnorm(i).lt.splinefit(i)-lowlim*signorm(i))
     &	    fitdata(i) = .FALSE.
C	if(.NOT.fitdata(i)) write(6,*) i,fluxnorm(i),splinefit(i)
C	  if (fluxnorm(i).gt.splinefit(i)+hilim*signorm(i))
C     &	    fitdata(i) = .FALSE.
 38	continue
 33   continue
c
c  normalize spectra and noise level to continuum, and output
 50   call fappend(ion_name,'cln',cln_file)
      open(unit=2,file=cln_file,status='unknown')
      fluxdec = 0.d0
      do 20 i=nwrap+1,ndata-nwrap
	if (splinefit(i).gt.fluxmax) splinefit(i) = fluxmax
C	if (fluxnorm(i)/splinefit(i).gt.fluxmax) 
C     &	  splinefit(i) = fluxnorm(i)/fluxmax
	fluxnorm(i) = fluxnorm(i)/splinefit(i)
	signorm(i) = signorm(i)/splinefit(i)
	write(2,'(2f12.4,2f12.6)') 
     &	 wave(i),vel(i),fluxnorm(i),signorm(i)
	fluxdec = fluxdec+(1.d0-MAX(splinefit(i),fluxnorm(i)/fluxmax))
 20   continue
      close(unit=2)  
	write(6,'(3x,a,f12.4)') 'Flux decrement: ',fluxdec/ndata
c
 1000 stop
c
      end
c
c-----------------------------------------------------------------------
c  third order polynomial fitting function
      subroutine FUNCS(x,coeff,y,dyda,n)
      parameter          (nmax = 5000)
	double precision fluxnorm,sigmax
	integer ndata
	common/fitcom/fluxnorm(nmax),sigmax,ndata
      double precision x,coeff(n),y,dyda(n)
      integer i,n

	coeff(2) = fluxnorm(ndata)-fluxnorm(1)
        do 5 i=3,n
     	  coeff(2) = coeff(2) - coeff(i)
 5	continue
	y = 0.d0
        do 10 i=1,n
	  y = y + coeff(i)*x**(i-1)
	  dyda(i) = x**(i-1)
 10     continue
c	write(6,*) coeff(1),coeff(2),coeff(3),coeff(4)
      return
      end

c-----------------------------------------------------------------------

c  numerical recipes routines to do Marquardt Method minimization (p 523)

      SUBROUTINE MRQMIN(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     *    COVAR,ALPHA,NCA,CHISQ,ALAMDA)
	implicit double precision (a-h,o-z)
      PARAMETER (MMAX=50)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MA),
     *  COVAR(NCA,NCA),ALPHA(NCA,NCA),ATRY(MMAX),BETA(MMAX),DA(MMAX)
      IF(ALAMDA.LT.0.)THEN
        KK=MFIT+1
        DO 12 J=1,MA
          IHIT=0
          DO 11 K=1,MFIT
            IF(LISTA(K).EQ.J)IHIT=IHIT+1
11        CONTINUE
          IF (IHIT.EQ.0) THEN
            LISTA(KK)=J
            KK=KK+1
          ELSE IF (IHIT.GT.1) THEN
            write(6,*) 'Improper permutation in LISTA'
          ENDIF
12      CONTINUE
        IF (KK.NE.(MA+1)) then
	  write(6,*) 'Improper permutation in LISTA'
	end if
        ALAMDA=0.001
        CALL MRQCOF(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NCA,CHISQ)
        OCHISQ=CHISQ
        DO 13 J=1,MA
          ATRY(J)=A(J)
13      CONTINUE
      ENDIF
      DO 15 J=1,MFIT
        DO 14 K=1,MFIT
          COVAR(J,K)=ALPHA(J,K)
14      CONTINUE
        COVAR(J,J)=ALPHA(J,J)*(1.+ALAMDA)
        DA(J)=BETA(J)
15    CONTINUE
      CALL GAUSSJ(COVAR,MFIT,NCA,DA,1,1)
c	call LUDCMP(COVAR,MFIT,NCA,INDX,D)
c	call LUBKSB(COVAR,MFIT,NCA,INDX,DA)
      IF(ALAMDA.EQ.0.)THEN
        CALL COVSRT(COVAR,NCA,MA,LISTA,MFIT)
        RETURN
      ENDIF
      DO 16 J=1,MFIT
        ATRY(LISTA(J))=A(LISTA(J))+DA(J)
16    CONTINUE
      CALL MRQCOF(X,Y,SIG,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,CHISQ)
      IF(CHISQ.LT.OCHISQ)THEN
        ALAMDA=0.1*ALAMDA
        OCHISQ=CHISQ
        DO 18 J=1,MFIT
          DO 17 K=1,MFIT
            ALPHA(J,K)=COVAR(J,K)
17        CONTINUE
          BETA(J)=DA(J)
          A(LISTA(J))=ATRY(LISTA(J))
18      CONTINUE
      ELSE
        ALAMDA=10.*ALAMDA
        CHISQ=OCHISQ
      ENDIF
      RETURN
      END


      SUBROUTINE MRQCOF(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NALP,
     *CHISQ)
	implicit double precision (a-h,o-z)
      PARAMETER (MMAX=50)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),ALPHA(NALP,NALP),BETA(MA),
     *    DYDA(MMAX),LISTA(MFIT),A(MA)
      DO 12 J=1,MFIT
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
12    CONTINUE
      CHISQ=0.
      DO 15 I=1,NDATA 
	CALL FUNCS(X(I),A,YMOD,DYDA,MA)
        SIG2I=1./(SIG(I)*SIG(I))
        DY=Y(I)-YMOD
        DO 14 J=1,MFIT
          WT=DYDA(LISTA(J))*SIG2I
          DO 13 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WT*DYDA(LISTA(K))
13        CONTINUE
          BETA(J)=BETA(J)+DY*WT
14      CONTINUE
        CHISQ=CHISQ+DY*DY*SIG2I
15    CONTINUE
      DO 17 J=2,MFIT
        DO 16 K=1,J-1
          ALPHA(K,J)=ALPHA(J,K)
16      CONTINUE
17    CONTINUE
      RETURN
      END


      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
	implicit double precision (a-h,o-z)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END


      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
	implicit double precision (a-h,o-z)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                write(6,*) 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) then
	do jj=1,n
	write(6,*) (a(jj,kk),kk=1,n)
	end do
	  write(6,*) 'Singular matrix.'
	end if
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

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
