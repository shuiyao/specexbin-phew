C=============================================================================

        subroutine model(m,n,imin,imax,iflag)
c
c  computes value of flux in region (imin,imax) given the presence
c  of lines numbered m through n.
c  iflag = 1  -> called from autofit
c  iflag = 2  -> called from minfit
c
        include 'vpdefs.h'
        include 'const.dek'
        integer m,n,imin,imax,i,j,iflag
        double precision b1,b2,b3,b4,x,y,u,v
        double precision tauconst

	if (iflag.eq.1) then
          do 5 j=1,ndata
            workflux(j) = 1.d0
 5        continue
	end if
        do 10 i=m,n
          if (iflag.eq.1) b1 = NHI(i)/1.d13
          if (iflag.eq.2) b1 = NHI(i)
          if (iflag.eq.1) b2 = lambda0*((1.d0+vel(centpix(i))/ckms)/
     &	    (1.d0-vel(centpix(i))/ckms))**0.5	! wavelength at center
          if (iflag.eq.2) b2 = lambda0*((1.d0+vline(i)/ckms)/
     &	    (1.d0-vline(i)/ckms))**0.5
C	write(6,*) waveln(centpix(i)),b2
          b3 = lambda0 / ckms * bpar(i)       ! doppler width
c	if (imax.eq.ndata) write(6,*) b1,b2,b3
          y = con2/b3
	  b4 = lambda0/(b2*b3)
          tauconst = con1*b1/b3
          do 20 j=imin,imax
            x = (waveln(j)-b2)*b4    ! wavelength in units of doppler width
	    if (fitflag.eq.1) then
              CALL cpf12(x,y,u,v)
              tau = tauconst*u            ! optical depth
	    else
	      tau = tauconst*exp(-x*x)
	    end if
            workflux(j) = workflux(j) * exp(-tau)
 20       continue
 10     continue

        return
        end

C=============================================================================


        double precision function model1(parm,np,veloc)
c
c  computes value of flux at a single velocity veloc given the presence
c  of (np/3) lines with parameters stored in parm
c
        include 'vpdefs.h'
        include 'const.dek'
	integer np
        double precision parm(np),veloc,lambda
        double precision b1,b2,b3,b4,x,y,u,v,tau

        model1 = 1.d0
c	lambda = lambda0*(1.d0+redz)*(1.d0+veloc/ckms)
	lambda = lambda0*((1.d0+veloc/ckms)/(1.d0-veloc/ckms))**0.5
        do 10 i=1,np,3
          b1 = parm(i+1)
c         b2 = lambda0*(1.d0+redz)*(1.d0+parm(i)/ckms)
          b2 = lambda0*((1.d0+parm(i)/ckms)/(1.d0-parm(i)/ckms))**0.5
          b3 = lambda0/ckms * parm(i+2)       ! doppler width
	  b4 = lambda0/(b2*b3)
          y = con2/b3
C          x = (lambda-b2)/b3/(1.d0+redz) ! wavelength in units of doppler width
          x = (lambda-b2)*b4    ! wavelength in units of doppler width
	    if (fitflag.eq.1) then
              CALL cpf12(x,y,u,v)
              tau = con1*b1/b3*u   ! optical depth
	    else
              tau = con1*b1/b3*exp(-x*x)   ! optical depth
	    end if
          model1 = model1 * exp(-tau)
c	write(6,'(8g10.3)') parm(i),parm(i+1),parm(i+2),veloc,redz,
c     &  con1,con2,tau
	if(parm(i).gt.ckms.OR.model1.gt.1.d0) write(6,*) 'PROBLEM:',
     &	parm(i),veloc,b1,b2,b3,b4
 10     continue

        return
        end


c..routine for profile.f
c
c  routine cpf12  computes the complex probability function w(z)
c                 using the 2-REGION n=12 approximation
c
c              w(z)   = exp(-z**2) * erfa(-iz) = u(x,y) + iv(x,y)
c    where
c              z      = x + iy
c              u(x,y) = Voigt "function" = H(x,y)
c              v(x,y) = complex component to w(z) 
c              x      = wave number scale in Doppler width units
c              y      = ratio of Lorentz to Doppler width 
c
c  computation is valid for the upper half plane (y>0); computational
c  accuracy is claimed to be maximum relative error < 2.0e-6 for u 
c  and < 5.0e-6 for v;  code is from method of J. Humlicek (1979) J.
c  Quant. Spectrosc. Radiat. Transfer _21_ 309.  typed in and trivially
c  modified by Christopher W. Churchill, March 1994.
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c..
c..
c..
      subroutine         cpf12(x,y,u,v)
c..
c..
c..compute the complex probability function w(z); see above comments
c.......................................................................
c..
c..declarations
C      implicit undefined (a-z)
      integer            i,m
      parameter          (m=6)
      double precision   b1,b2,t(m),c(m),s(m)
      double precision   u,v,x,y,d,d1,d2,d3,d4,r,r2,y1,y2,y3         
c..
c..load the data ...
c..harmonic roots and weights of the n-point Gauss-Hermite formula,
c..chosen to keep absolute error of u and v < 1.0e-7 
      data b1/ 0.850d0/
      data  t/ 0.314240376  ,  0.947788391  ,  1.59768264   ,
     @         2.27950708   ,  3.02063703   ,  3.8897249    /
      data  c/ 1.01172805   , -0.75197147   ,  1.2557727e-2 ,
     @         1.00220082e-2, -2.42068135e-4,  5.00848061e-7/
      data  s/ 1.393237     ,  0.231152406  , -0.155351466  ,
     @         6.21836624e-3,  9.19082986e-5, -6.27525958e-7/
c..
c..initialize this pass; define the x<18.1 REGION boundary
      u  = 0.0d0
      v  = 0.0d0
      y1 = y + 1.50d0
      y2 = y1*y1
      b2 = 18.10d0*y + 1.650d0
c..
c..
      if ((y.gt.b1).or.(abs(x).lt.b2)) then
c..
c..eq.(6);  REGION 1  y > (x-1.65)/18.1  for x<=18.1 
c                     y > 0.85           for x> 18.1
       do 3 i=1,m
        r  = x - t(i)
        d  = 1.0d0/(r*r+y2)
        d1 = y1*d
        d2 = r*d
        r  = x + t(i)
        d  = 1.0d0/(r*r+y2)        
        d3 = y1*d
        d4 = r*d
        u  = u + c(i)*(d1+d3) - s(i)*(d2-d4)
        v  = v + c(i)*(d2+d4) + s(i)*(d1-d3)
 03    continue
       return
c..
      else
c..
c..eq.(11); REGION 2  y <= (x-1.65)/18.1  for x<=18.1 
c..                   y <= 0.85           for x> 18.1
       if (abs(x).lt.12.0d0) u = exp(-x*x)
       y3 = y + 3
       do 5 i=1,m
        r  = x - t(i)
        r2 = r*r
        d  = 1.0d0/(r2+y2)
        d1 = y1*d
        d2 = r*d        
        u  = u + y*(c(i)*(r*d2-1.50d0*d1)+s(i)*y3*d2)/(r2+2.250d0)
        r  = x + t(i)
        r2 = r*r
        d  = 1.0d0/(r2+y2)        
        d3 = y1*d
        d4 = r*d
        u  = u + y*(c(i)*(r*d4-1.50d0*d3)-s(i)*y3*d4)/(r2+2.250d0)
        v  = v + c(i)*(d2+d4) + s(i)*(d1-d3)
 05    continue
       return
c..
      end if
      end
c..
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c..eof
C=============================================================================
c..
c..
      subroutine         cpf15(x,y,u,v)
c..
c..this routine computes the complex probability function w(z)
c..
c..       w(z)   = exp(-z**2) * erfc(-iz) = u(x,y) + iv(x,y)
c..where
c..       z      = x + iy
c..       u(x,y) = Voigt "function" = H(x,y)
c..       v(x,y) = complex component to w(z) 
c..       x      = wave number scale in Doppler width units
c..       y      = ratio of Lorentz to Doppler width 
c..
c..computation is valid for the upper half plane (y>0); computational
c..accuracy is claimed to be maximum relative error < 2.0e-6 for u 
c..and < 5.0e-6 for v;  code is from method of J. Humlicek (1979) J.
c..Quant. Spectrosc. Radiat. Transfer _21_ 309.  typed in and trivially
c..modified by Christopher W. Churchill, March 1994.
c..
c..the 2-REGION n=12 approximation
c..
c..local declarations
      integer            i,m
      parameter          (m=6)
      double precision   b1,b2,t(m),c(m),s(m)
      double precision   u,v,x,y,d,d1,d2,d3,d4,r,r2,y1,y2,y3         
c..
c..load the data ...
c..harmonic roots and weights of the n-point Gauss-Hermite formula,
c..chosen to keep absolute error of u and v < 1.0e-7 
      data b1/ 0.850d0/
      data  t/ 0.314240376  ,  0.947788391  ,  1.59768264   ,
     @         2.27950708   ,  3.02063703   ,  3.8897249    /
      data  c/ 1.01172805   , -0.75197147   ,  1.2557727e-2 ,
     @         1.00220082e-2, -2.42068135e-4,  5.00848061e-7/
      data  s/ 1.393237     ,  0.231152406  , -0.155351466  ,
     @         6.21836624e-3,  9.19082986e-5, -6.27525958e-7/
c..
c..initialize this pass; define the x<18.1 REGION boundary
      u  = 0.0d0
      v  = 0.0d0
      y1 = y + 1.50d0
      y2 = y1*y1
      b2 = 18.10d0*y + 1.650d0
c..
c..
      if ((y.gt.b1).or.(abs(x).lt.b2)) then
c..
c..eq.(6);  REGION 1  y > (x-1.65)/18.1  for x<=18.1 
c..                   y > 0.85           for x> 18.1
       do 3 i=1,m
        r  = x - t(i)
        d  = 1.0d0/(r*r+y2)
        d1 = y1*d
        d2 = r*d
        r  = x + t(i)
        d  = 1.0d0/(r*r+y2)        
        d3 = y1*d
        d4 = r*d
        u  = u + c(i)*(d1+d3) - s(i)*(d2-d4)
        v  = v + c(i)*(d2+d4) + s(i)*(d1-d3)
 03    continue
       return
c..
      else
c..
c..eq.(11); REGION 2  y <= (x-1.65)/18.1  for x<=18.1 
c..                   y <= 0.85           for x> 18.1
       if (abs(x).lt.12.0d0) u = exp(-x*x)
       y3 = y + 3
       do 5 i=1,m
        r  = x - t(i)
        r2 = r*r
        d  = 1.0d0/(r2+y2)
        d1 = y1*d
        d2 = r*d        
        u  = u + y*(c(i)*(r*d2-1.50d0*d1)+s(i)*y3*d2)/(r2+2.250d0)
        r  = x + t(i)
        r2 = r*r
        d  = 1.0d0/(r2+y2)        
        d3 = y1*d
        d4 = r*d
        u  = u + y*(c(i)*(r*d4-1.50d0*d3)-s(i)*y3*d4)/(r2+2.250d0)
        v  = v + c(i)*(d2+d4) + s(i)*(d1-d3)
 05    continue
       return
c..
      end if
      end




