
	subroutine getlines(infile)

	include 'vpdefs.h'
        character*80 infile
	integer j,k
c
c  read in fit lines
	write(6,*) 'inputting ',infile(1:50)
        open(unit=1,file=infile,status='old')
        read(1,*)
        read(1,'(1x,a10,i10)') ion_name,nlines
        do 30 j=1,nlines
      read(1,*) k,NHI(j),vline(j),bpar(j),tempout(j),rhoout(j),metout(j)
          dvline(j) = 0.d0
          dNHI(j) = 0.d0
          dbpar(j) = 0.d0
 30     continue
        close(1)

	return
	end
