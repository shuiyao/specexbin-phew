
        subroutine nparmin(x,y,sigfit,nd,parm,np,listparm,npv,chisq)

	include 'vpdefs.h'
        double precision parm(*),x(*),y(*),sigfit(*)
        integer np,npv,nd,listparm(*)
	integer i,ipv,j
        double precision dchisqdp(maxparm),dfdp(maxparm),
     &		parmold(maxparm)
        double precision fact,f,chisq,chisqold
	double precision chi2

c  determines which of the parameters parm have significant effect
c  on chisqm and lists those ipv parameters in listparm

c	chisqold = chi2(x,y,sigfit,nd,parm,np)
	chisqold = chisq
        do 5 j=1,np
	  parmold(j) = parm(j)
	  listparm(j) = j
	  dchisqdp(j) = 0.d0
 5	continue
c
        do 10 i=1,nd
          call FUNCS(x(i),parm,f,dfdp,np)
          fact = 2.d0*(y(i)-f)/(sigfit(i)*sigfit(i))
          do 20 j=1,np
            dchisqdp(j) = dchisqdp(j) - dfdp(j)*fact
 20       continue
 10     continue
c	write(6,'(6g12.4)') (dchisqdp(j),j=1,6)

	CALL PIKSR2(np,dchisqdp,listparm)
c	ipv = 1
c	do 30 while (ABS(dchisqdp(ipv)).gt.0.1*ABS(dchisqdp(1)))
c 	write(6,*) 'using ',ipv,dchisqdp(ipv)
c	   ipv = ipv+1
c 30	continue
c	ipv = ipv-1
c	write(6,*) 'starting ipv',ipv

	totalred = 1.d0
	ipv = 1
 45	continue
	freduce = totalred
	do 50 j=ipv,ipv
	  parm(listparm(j)) = parmold(listparm(j))-
     &	    freduce*chisqold/dchisqdp(j)
	  CALL limitparm(parm,np)
 50	continue
	chisq = chi2(x,y,sigfit,nd,parm,np)
	if (chisq.gt.chisqold.AND.totalred.gt.fvarymin) then
	  totalred = 0.5d0*totalred
c	  write(6,*) 'reducing totalred ',totalred,chisq
	  goto 45
	end if
	if (chisq.gt.chisqold.AND.totalred.le.fvarymin
     &	  .AND.ipv.lt.np) then
	  totalred = 1.d0
	  ipv = ipv + 1
c	  write(6,*) 'increasing ipv ',ipv,chisq
	  goto 45
	end if
	if (chisq.gt.chisqold) then
	  chisq = chisqold
          do 60 j=1,np
	    parm(j) = parmold(j)
 60	  continue
	else
c	  write(6,*) 'N-par min ',listparm(ipv),ipv,chisq
	end if

        return
        end

