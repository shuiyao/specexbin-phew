
	SUBROUTINE fixparm(j,parm,parerr,np,npv,listparm,index)
c
c  called when covariance is 0 for parameter j.
c  fixes parameter j at current value.
c
        include 'vpdefs.h'
        include 'const.dek'
	double precision parm(*),parerr(*)
	integer listparm(*),index(*),np,npv,j
	integer k

c  if line is very small, freeze ALL line parameters
        write(6,*) 'covar=0: param,line',listparm(j),index(listparm(j))
        k = (listparm(j)-1)/3
        k = 3*k+1
        if (parm(k).lt.vmin.OR.parm(k).gt.vmax
     &    .OR.parm(k+1).le.NHImin
     &    .OR.parm(k+2).le.bparmin) then
          parm(k) = 0.5*(vmin+vmax)
          parm(k+1) = NHImin
          parm(k+2) = bparmin
c  change listparm to exclude non-varying parameters
          do 74 ik=k,npv-3
            listparm(ik) = listparm(ik+3)
 74       continue
          listparm(npv-2) = k
          listparm(npv-1) = k+1
          listparm(npv) = k+2
          npv = npv-3
        else
c  if line not small, then freeze only parameter with 0 covar
          do 79 k=j,npv-1
            listparm(k) = listparm(k+1)
 79       continue
          listparm(npv) = j
          npv = npv-1
        end if

	return
	end
