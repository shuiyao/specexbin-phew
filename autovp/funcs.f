c=======================================================================

        subroutine FUNCS(veloc,parm,f,dfdp,np)

c  find value of flux f and derivatives dfdp at veloc
c  given the presence of np/3 Voigt profile lines.  For each line,
c  parm(1) is vel of line center, parm(2) is NHI, parm(3) is b parameter.

        include 'vpdefs.h'
        include 'const.dek'
        integer np
        double precision veloc,f,dfdp(np),parm(np)
        double precision model1,dparm

c  check parameters for out of bounds
	CALL limitparm(parm,np)
c  compute model flux at veloc
        f = model1(parm,np,veloc)
c  compute derivatives at veloc wrt each parameter
        do 10 i=1,np
          dfdp(i) = dparm(parm,np,veloc,i)
 10     continue

        return
        end
c
c=======================================================================

        double precision function dparm(parm,np,veloc,ip)
c
c  computes derivative of model dparm at veloc wrt parameter ip from list parm.
c  compute derivative numerically on successively smaller intervals
c  until convergence.
c
        include 'vpdefs.h'
        include 'const.dek'
        integer np,iter
        double precision veloc,parm(np)
        double precision olddparm,deltah,parmold
        double precision modelhi,modellow,model1

c  remember original value of parameter
        parmold = parm(ip)
c  compute the amount you want to vary parm(ip)
        deltah = ABS(tolderiv*parm(ip))
	if (parm(ip).eq.0.d0) deltah = tolderiv
c  velocity deltah is tied to bparmax, not parm value
	if( 3*((ip-1)/3)+1.eq.ip ) deltah = tolderiv*bparmax
C	deltah = MIN(deltap(ip),deltah)
c  compute first guess at derivative
	iter = 0
 5      continue
        parm(ip) = parmold+deltah
        modelhi = model1(parm,np,veloc)
        parm(ip) = parmold-deltah
        modellow = model1(parm,np,veloc)
c  if difference (modelhi-modellow) killed by roundoff, increase deltah
        if (modelhi.eq.modellow.AND.deltah.lt.0.5*parmold) then
          deltah = 2.d0*deltah
          parm(ip) = parmold
	  if( 3*((ip-1)/3)+1.eq.ip .AND. deltah.gt.bparmax ) goto 7
	  iter = iter+1
          goto 5
        endif
c  compute first guess at derivative
 7	continue
        dparm = 0.5d0*(modelhi-modellow)/deltah
        olddparm = 2.d0*dparm
        parm(ip) = parmold
c  reduce deltah until convergence of derivative
        do 10 while(ABS((dparm-olddparm)/dparm).gt.tolderiv.AND.
     &    dparm.ne.0.)
          olddparm = dparm
          deltah = 0.5d0*deltah

          parm(ip) = parmold+deltah
          modelhi = model1(parm,np,veloc)
          parm(ip) = parmold-deltah
          modellow = model1(parm,np,veloc)
          dparm = 0.5d0*(modelhi-modellow)/deltah
	  parm(ip) = parmold
 10     continue
        if(dparm.eq.0.d0) then
c         write(6,*) 'deriv=0! ',ip,parm(ip),deltah,olddparm,
c     &		modelhi,modellow
          dparm = olddparm
        end if
c  reset original value of parm
        parm(ip) = parmold
	deltap(ip) = deltah

        return
        end

