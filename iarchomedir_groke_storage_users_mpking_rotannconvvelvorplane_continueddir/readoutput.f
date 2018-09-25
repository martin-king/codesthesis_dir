      implicit real*8(a-h,o-z)

      real*8 th(0:417,1,0:3584)

!mpkingdifferent        open(11,file='temp1e103e5r4.dat',form='unformatted')
!mpk2d         read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
!mpkingdifferent        read(11)(((th(i,j,k),i=0,353),j=1,1),k=0,1792)
        open(11,file='temp8e78e4pr1-1r10.dat',status='unknown')
        read(11,201)(((th(i,j,k),i=0,241),j=1,1),k=0,1792)
        close(11)

        open(11,file='plotmap.dat',status='unknown')
        write(11,201)(((th(i,j,k),i=0,241),j=1,1),k=0,1792)
        close(11)
  
!mpkdifferent 201    format(1(1x,e16.8))
  201  format(1(1x,e21.15))
        stop
        end
