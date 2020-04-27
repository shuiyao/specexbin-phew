
	SUBROUTINE initparm(parm,parerr,index,listparm,np,npv)
c
c  finds all lines contained within vel region (vmin,vmax).
c  loads parameters for those lines into parm(np), which are
c  the parameters varied by minimization routines.
c  listparm contains indices of npv parameters to be varied;
c  for now, vary all parameters (npv=np).
c  
        include 'vpdefs.h'
        include 'const.dek'
	double precision parm(*),parerr(*)
	integer index(*),listparm(*),np,npv,j

        np = 0
        do 60 j=1,nlines
c  find all lines contained in detection region
          if (vline(j).ge.vmin.AND.vline(j).le.vmax) then
            np = np + 1
            parm(np) = vline(j)
            parerr(np) = dvline(j)
            index(np) = j
            np = np + 1
            parm(np) = NHI(j)
            parerr(np) = dNHI(j)
            index(np) = j
            np = np + 1
            parm(np) = bpar(j)
            parerr(np) = dbpar(j)
            index(np) = j
          end if
 60     continue
c  fit all lines at once (for now)
        do 65 j=1,np
          listparm(j) = j
	  parm0(j) = parm(j)
	  deltap(j) = 0.01*parm(j)
	  if( 3*((j-1)/3)+1.eq.j ) deltap(j) = tolderiv*bparmax
 65     continue
	npv = np

c  output initial values
        write(6,*) '>>>',np/3,' lines:'
        do 67 j=1,np,3
          write(6,'(i3,3(2f12.4,1x))') index(j),
     &	    parm(j+1),parerr(j+1),
     &      parm(j),parerr(j),
     &      parm(j+2),parerr(j+2)
 67     continue

	return
	end
