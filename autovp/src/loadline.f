
	SUBROUTINE loadline(parm,parerr,np,index)

        include 'vpdefs.h'
        include 'const.dek'
	double precision parm(*),parerr(*)
	integer index(*),np,j

        do 84 j=1,np,3
          vline(index(j)) = parm(j)
          NHI(index(j+1)) = parm(j+1)
          bpar(index(j+2)) = parm(j+2)
          dvline(index(j)) = parerr(j)
          dNHI(index(j+1)) = parerr(j+1)
          dbpar(index(j+2)) = parerr(j+2)
	  if (lineinfo(index(j)).ne.-2) lineinfo(index(j)) = 1
 84     continue

	return
	end
