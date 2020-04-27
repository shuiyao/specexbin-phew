      subroutine polyfit(imin,imax)
c
c  fits a piecewise 3rd-order polynomial to data in region (imin,imax)
c  Does nint pixels at a time, plus minint on either side for overlap.
c
	include 'vpdefs.h'
      integer            fitorder
      parameter          (fitorder=4)
      integer            i,j,jmin,jmax
      double precision   fluxfit(nmax),sigfit(nmax)
      double precision   coeff(fitorder)
      double precision	 work1(nmax,fitorder),work2(fitorder,fitorder),
     &			 work3(fitorder)
      integer		 listcoeff(fitorder),nfit
      double precision   chisq
      double precision	 velmax
	logical done
c
c  initialize fit:  fits all components of polynomial
      interval = MIN(nint,(imax-imin)/5)
      interval = MAX(interval,minint)
      do 30 i=1,fitorder
	listcoeff(i) = i
 30   continue
      jmin = imin
      jmax = jmin+interval-1
      done = .FALSE.
      do 20 while (.NOT.done)
c  first guess: linear
	velmax = 0.d0
	do 10 i=jmin,jmax
          coeff(1) = coeff(1) + resid(i)
	  velmax = MAX(vel(i),velmax)
 10	continue
	coeff(1) = coeff(1)/(jmax-jmin+1)
	coeff(2) = velmax*(resid(jmax)-resid(jmin))/(vel(jmax)-vel(jmin))
	coeff(3) = 0.d0
	coeff(4) = 0.d0
c  load appropriate data into arrays to be fitted
	nfit = 0
	do 34 i=jmin-minint,jmax+minint
	  if (i.gt.0.AND.i.le.ndata) then
	    nfit = nfit+1
	    velfit(nfit) = vel(i)/velmax
	    fluxfit(nfit) = resid(i)
	    sigfit(nfit) = noise(i)
	  endif
 34	continue
c  fit polynomial to data set
	CALL SVDFIT(velfit,fluxfit,sigfit,nfit,coeff,fitorder,work1,
     &	  work2,work3,nmax,fitorder,chisq)
c	write(6,'(i7,6g12.4)') jmin,chisq,(coeff(k),k=1,fitorder)
c  compute polynomial value at each pixel
	do 36 i=jmin,MIN(jmax,imax)
	  minflux(i) = 0.d0
	  dfdv(i) = 0.d0
	  d2fdv2(i) = 0.d0
	  do 38 j=1,fitorder
	    if (vel(i).eq.0.d0) goto 38
	    minflux(i)=minflux(i)+coeff(j)*(vel(i)/velmax)**(j-1)
	    dfdv(i)=dfdv(i)+(j-1)*coeff(j)*(vel(i)/velmax)**(j-2)
	    d2fdv2(i)=d2fdv2(i)+(j-1)*(j-2)*coeff(j)*
     &		(vel(i)/velmax)**(j-3)
 38	  continue
 36	continue
	jmin = jmax+1
	jmax = jmin+interval-1
	if (jmin.ge.imax) done = .TRUE.
 20   continue
c
c  set minflux to lie between 0 and 1.
      do 50 i=imin,imax
	minflux(i) = MAX(minflux(i),0.d0)
	minflux(i) = MIN(minflux(i),1.d0)
 50   continue

c        call fappend(ion_name,'min',sig_file)
c        open(unit=2,file=sig_file,status='unknown')
c        do 60 i=1,ndata
c          write(2,'(2f12.4,2f12.6)') vel(i),dfdv(i),minflux(i),
c     &      d2fdv2(i)
c 60     continue
c        close(2)
c	pause ' '

 1000 return
c
      end
c
c-----------------------------------------------------------------------
c  polynomial fitting function of order (fitorder-1)
      subroutine FUNCPOLY(x,afunc,n)
      double precision x,afunc(n)
      integer i,n

      afunc(1) = 1.d0
      do 10 i=2,n
	afunc(i) = afunc(i-1)*x
 10   continue
      return
      end

c-----------------------------------------------------------------------

      SUBROUTINE SVDFIT(X,Y,SIG,NDATA,A,MA,U,V,W,MP,NP,CHISQ)
      implicit double precision (a-h,o-z)
      PARAMETER(NMAX=1000,MMAX=50,TOL=1.E-5)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),V(NP,NP),
     *    U(MP,NP),W(NP),B(NMAX),AFUNC(MMAX)
      DO 12 I=1,NDATA
        CALL FUNCPOLY(X(I),AFUNC,MA)
        TMP=1./SIG(I)
        DO 11 J=1,MA
          U(I,J)=AFUNC(J)*TMP
11      CONTINUE
        B(I)=Y(I)*TMP
12    CONTINUE
      CALL SVDCMP(U,NDATA,MA,MP,NP,W,V)
      WMAX=0.
      DO 13 J=1,MA
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,MA
        IF(W(J).LT.THRESH)W(J)=0.
14    CONTINUE
      CALL SVBKSB(U,W,V,NDATA,MA,MP,NP,B,A)
      CHISQ=0.
      DO 16 I=1,NDATA
        CALL FUNCPOLY(X(I),AFUNC,MA)
        SUM=0.
        DO 15 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
15      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
16    CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------

      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=100)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.50) THEN
		write(6,*) 'No convergence in 50 iterations'
		RETURN
	  ENDIF
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------

      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
      implicit double precision (a-h,o-z)
      PARAMETER (NMAX=100)
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      DO 12 J=1,N
        S=0.
        IF(W(J).NE.0.)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END
