
        subroutine addline(ireg,x,y,sigfit,nd,parm,np)

        include 'vpdefs.h'
        double precision parm(*),x(*),y(*),sigfit(*)
        integer np,nd,ireg
        double precision maxchisq
	double precision xsubreg(20),ysubreg(20),sigsubreg(20)
        double precision chi2,model1

	maxchisq = 0.d0
        do 10 j=region(ireg,1)+nwidth,region(ireg,2)-nwidth
	  do 20 k=1,2*nwidth+1
	    xsubreg(k) = x(j-nwidth+k-region(ireg,1))
	    f = model1(parm,np,xsubreg(k))
	    ysubreg(k) = MIN(f,y(j-nwidth+k-region(ireg,1)))
	    sigsubreg(k) = sigfit(j-nwidth+k-region(ireg,1))
 20	  continue
	  chisqtemp = chi2(xsubreg,ysubreg,sigsubreg,2*nwidth+1,parm,np)
	  if (chisqtemp.gt.maxchisq) then
	    iadd = j
	    maxchisq = chisqtemp
	  end if
 10     continue
	nlines = nlines+1
	NHI(nlines) = 0.5d0
	bpar(nlines) = 20.d0
	vline(nlines) = vel(iadd)
	lineinfo(nlines) = -2
	write(6,*) 'adding line ',nlines,vline(nlines)

	return
	end
