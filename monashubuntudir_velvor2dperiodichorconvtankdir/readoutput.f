      implicit real*8(a-h,o-z)

      real*8 th(0:500,1,0:1000)

        open(11,file='temp.dat',form='unformatted')
!mpk2d         read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
        read(11)(((th(i,j,k),i=0,321),j=1,1),k=0,353)
        close(11)

        open(11,file='plotmap.dat',status='unknown')
        write(11,201)(((th(i,j,k),i=0,321),j=1,1),k=0,353)
        close(11)
  
 201    format(1(1x,e16.8))

        stop
        end
