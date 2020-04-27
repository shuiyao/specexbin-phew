
        subroutine oneparmin(x,y,sigfit,nd,parm,np,listparm,npv,chisq)

	include 'vpdefs.h'
        double precision parm(*),x(*),y(*),sigfit(*)
        integer np,npv,nd,listparm(*)
        double precision dchisqdp(maxparm),dfdp(maxparm),
     &		parmold(maxparm),chisqdiff(maxparm)
        double precision fact,f,chisq,chisqold
	double precision chi2

c  determines which of the parameters parm have significant effect
c  on chisqm and lists those npv parameters in listparm

c	chisqold = chi2(x,y,sigfit,nd,parm,np)
	chisqold = chisq
        do 5 j=1,np
	  parmold(j) = parm(j)
	  listparm(j) = j
	  chisqdiff(j) = 0.d0
	  dchisqdp(j) = 0.d0
 5	continue
c
        do 10 i=1,nd
          call FUNCS(x(i),parm,f,dfdp,np)
          fact = (y(i)-f)/(sigfit(i)*sigfit(i))
          do 20 j=1,np
            dchisqdp(j) = dchisqdp(j) - 2.d0*dfdp(j)*fact
 20       continue
c	write(6,'(6g12.4)') (dchisqdp(j),j=1,6)
 10     continue

	totalred = 1.d0
	npv = 1
 45	continue
	j = npv
	do 50 i=1,np
	  parm(i) = parmold(i)
 50	continue
	parm(j) = parmold(j)-
     &	    totalred*chisqold/dchisqdp(j)
	CALL limitparm(parm,np)
	chisq = chi2(x,y,sigfit,nd,parm,np)
	if (totalred.ge.fvarymin) then
	  totalred = 0.5d0*totalred
	  chisqdiff(j) = MAX(chisqdiff(j),chisqold-chisq)
	  goto 45
	end if
	if (npv.lt.np) then
	  totalred = 1.d0
	  npv = npv + 1
	  goto 45
	end if
c  all parameters tried; find best one (param which lowered chisq most)
	chisqdmax = chisqdiff(1)
	j = 1
	do 70 i=1,npv
	  if (chisqdiff(i).gt.chisqdmax) then
	    chisqdmax = chisqdiff(i)
	    j = i
	  end if
 70	continue
c	write(6,'(6g12.4)') (chisqdiff(i),i=1,5),chisqdmax
c  if chisq not improved, reset to old parameters and return
	if (chisqdmax.le.0.0) then
	  chisq = chisqold
          do 60 j=1,np
	    parm(j) = parmold(j)
 60	  continue
	  return
	end if
c  recompute chisq for best param
	do 80 i=1,np
	  parm(i) = parmold(i)
 80	continue
	totalred = 1.d0
 75	continue
	parm(j) = parmold(j)-
     &	    totalred*chisqold/dchisqdp(j)
	CALL limitparm(parm,np)
	chisq = chi2(x,y,sigfit,nd,parm,np)
	if (chisqold-chisq.lt.chisqdiff(j)) then
	  totalred = 0.5d0*totalred
c	  write(6,*) 'reducing totalred ',totalred,chisq
	  goto 75
	end if
c	write(6,*) '1-par min',j,npv,chisq

        return
        end


	double precision function chi2(x,y,sigfit,nd,parm,np)

        double precision parm(*),x(*),y(*),sigfit(*)
	double precision model1
        integer nd,np,i

	chi2 = 0.d0
	do 10 i=1,nd
          f = model1(parm,np,x(i))
	  chi2 = chi2 + ((y(i)-f)/sigfit(i))**2
 10	continue

	return
	end

	subroutine limitparm(parm,np)

	include 'vpdefs.h'
	double precision parm(*)
	integer i,np

        do 5 i=1,np,3
          if (parm(i).lt.vmin) parm(i) = vmin  ! don't let v go 
          if (parm(i).gt.vmax) parm(i) = vmax  ! outside region.
          if (parm(i+1).lt.NHImin) parm(i+1) = NHImin   ! limit NHI value.
C          if (parm(i+1).gt.NHImax) parm(i+1) = NHImaxsat ! limit NHI value.
          if (parm(i+2).lt.bparmin) parm(i+2) = bparmin ! min b value.
          if (parm(i+2).gt.bparmax) parm(i+2) = bparmax ! max b value
          if (parm(i+2).gt.1.25*parm0(i+2)) parm(i+2) = 1.25*parm0(i+2) ! max b
 5      continue

	return
	end

