c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      implicit real*8(a-h,o-z)
      parameter(kr=65,kz=65,kp=160,mgen=100000)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 ueln(0:kr,0:kz,0:kp),veln(0:kr,0:kz,0:kp)
      real*8 weln(0:kr,0:kz,0:kp)
      character*20 fname1
c
      print*,'filename for velocity u: '
      read*,fname1
      open(11,file=fname1,form='unformatted')
      read(11)nr,nz,np
      read(11)rmax,r0,zmax,zmin
      read(11)(((uel(i,j,k),i=0,nr),j=0,nz+1),k=0,np+1)
      close(11)
c
      print*,'filename for velocity v: '
      read*,fname1
      open(11,file=fname1,form='unformatted')
      read(11)nr,nz,np
      read(11)rmax,r0,zmax,zmin
      read(11)(((vel(i,j,k),i=0,nr+1),j=0,nz+1),k=0,np)
      close(11)
c
      print*,'filename for velocity w: '
      read*,fname1
      open(11,file=fname1,form='unformatted')
      read(11)nr,nz,np
      read(11)rmax,r0,zmax,zmin
      read(11)(((wel(i,j,k),i=0,nr+1),j=0,nz),k=0,np+1)
      close(11)
c
      print*,'filename for temperature: '
      read*,fname1
      open(11,file=fname1,form='unformatted')
      read(11)nr,nz,np
      read(11)rmax,r0,zmax,zmin
      read(11)(((th(i,j,k),i=0,nr+1),j=0,nz+1),k=0,np+1)
      close(11)
c   just initialising
      do 10 i=0,nr+1
      do 10 j=0,nz+1
      do 10 k=0,np+1
        ueln(i,j,k)=0.d0
        weln(i,j,k)=0.d0
        veln(i,j,k)=0.d0
        thn(i,j,k)=0.d0
   10 continue
c   transforming coordinates of u,v,w & th to (0:nr,0:nz,0:np)
      print*,nr,nz,np
      print*,rmax,r0,zmax,zmin
      do 20 i=1,nr-1
      do 20 j=1,nz-1
      do 20 k=1,np-1
        ip=i+1
        jp=j+1
        ka=k+1
!mpk         if(k.eq.np)ka=1
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
        ueln(i,j,k)=(uel(i,j,k)+uel(i,jp,k)+uel(i,j,ka)+uel(i,jp,ka))
     +             /4.d0
        weln(i,j,k)=(wel(i,j,k)+wel(ip,j,k)+wel(i,j,ka)+wel(ip,j,ka))
     +             /4.d0
        veln(i,j,k)=(vel(i,j,k)+vel(ip,j,k)+vel(i,jp,k)+vel(ip,jp,k))
     +             /4.d0
        thn(i,j,k)=(th(i,j,k)+th(ip,j,k)+th(i,jp,k)+th(ip,jp,k)
     +            +th(i,j,ka)+th(ip,j,ka)+th(i,jp,ka)+th(ip,jp,ka))/8.d0
   20 continue
c
      do 30 i=0,nr
      do 30 j=0,nz
      do 30 k=0,np
        ip=i+1
        jp=j+1
        ka=k+1
!mpk        if(k.eq.np)ka=1
c     thermal boundaries at discs and radial walls
        thn(i,0,k)=(th(i,0,k)+th(ip,0,k)+th(i,0,ka)+th(ip,0,ka))/4.d0
        thn(i,nz,k)=(th(i,nz,k)+th(ip,nz,k)+th(i,nz,ka)+th(ip,nz,ka))/
     +              4.d0
        thn(i,j,0)=(th(i,j,0)+th(ip,j,0)+th(i,jp,0)+th(ip,jp,0))/4.d0
        thn(i,j,np)=(th(i,j,np)+th(ip,j,np)+th(i,jp,np)+th(ip,jp,np)
     +              )/4.d0
        thn(0,j,k)=0.d0
        thn(nr,j,k)=1.d0
   30 continue
c
      open(11,file='graf.file',form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((ueln(i,j,k),i=0,nr),j=0,nz),k=0,np)
      write(11)(((veln(i,j,k),i=0,nr),j=0,nz),k=0,np)
      write(11)(((weln(i,j,k),i=0,nr),j=0,nz),k=0,np)
      write(11)(((thn(i,j,k),i=0,nr),j=0,nz),k=0,np)
      close(11)
c
      stop
      end
