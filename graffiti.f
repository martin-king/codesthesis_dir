      implicit none

      integer idim,jdim,kdim

      parameter (idim=33)
      parameter (jdim=33)
      parameter (kdim=129)

      real*8 dum(jdim,kdim)
      real*8 x(idim,jdim,kdim)
      real*8 y(idim,jdim,kdim)
      real*8 z(idim,jdim,kdim)
      real*8 vx(idim,jdim,kdim),vy(idim,jdim,kdim)
      real*8 t(idim,jdim,kdim),vz(idim,jdim,kdim)
      real*8 p(idim,jdim,kdim)
      real*8 dummy(idim,jdim,kdim)
      real*8 rmax,r0,zmax,zmin
      real tmax(kdim)
      real ximax,ximin,etamax,etamin

      character*4 names(9)
      character*80 title(8)

      character*11 blfile
      integer kr,kf,iframe,incr,nr,nz,np

      integer idrop,jdrop,kdrop,k3d,j3d,i3d

      integer imax,jmax,kmax, i,j,k, nrec, bound_type(3)
      integer isurfs(8,2), jsurfs(8), ksurfs(8,3)

      data names/'x','y','z','vx','vy','vz','p','t','extr'/


c------------------------------------------------- x & y
c
      open(11,file='grafdat',form='unformatted')
      read(11)nr,nz,np
      read(11)rmax,r0,zmax,zmin
      read(11)(((vx(j+1,i+1,k+1),i=0,nr),j=0,nz),k=0,np)
      read(11)(((vy(j+1,i+1,k+1),i=0,nr),j=0,nz),k=0,np)
      read(11)(((vz(j+1,i+1,k+1),i=0,nr),j=0,nz),k=0,np)
      read(11)(((t(j+1,i+1,k+1),i=0,nr),j=0,nz),k=0,np)
      close(11)
c
      idrop=1
      jdrop=1
      kdrop=1
      imax=nz+1
      jmax=nr+1
      kmax=np+1
      print*,imax,jmax,kmax
c
      i3d=0
      do 501 i=1,imax,idrop
        i3d=i3d+1
        j3d=0
      do 501 j=1,jmax,jdrop
        j3d=j3d+1
        k3d=0
      do 501 k=1,kmax,kdrop
        k3d=k3d+1
  501 continue

c------------------------------------------------- 
c
      bound_type(1) = 1
      bound_type(2) = 2
      bound_type(3) = 1

      title(1) = 'Title line 1'
      title(2) = 'Title line 2'
      title(3) = 'Title line 3'
      title(4) = 'Title line 4'
      title(5) = 'Title line 5'
      title(6) = 'Title line 6'
      title(7) = 'Title line 7'
      title(8) = 'Title line 8'

c------------------------------------------------- 
      open (unit=91,file='graft1.file',recl=20,form='unformatted',
     &      access='direct')

      call writetop('STAN', 1, 'EXAM', i3d, j3d, k3d,
     &              2, 0, 1, 1, 1, 1,
     &              0.0, 0.0, 0.0, 0.0, 0.0, title,
     &              bound_type, isurfs, jsurfs, ksurfs, nrec)
c------------------------------------------------- 
c
      do 5 i=1,i3d
      do 5 j=1,j3d
      do 5 k=1,k3d
         x(i,j,k) = x(1,j,k)
         y(i,j,k) = y(1,j,k)
         z(i,j,k) = z(1,j,k)
   5  continue
c-------------------------------------------------

      do 11 i=1,i3d
      do 11 j=1,j3d
      do 11 k=1,k3d
         x(i,j,k)=dble(i-1)/dble(i3d-1)*zmax
         y(i,j,k)=r0+dble(j-1)/dble(j3d-1)*(rmax-r0)
         z(i,j,k)=(3.141592654/0.5)*float(k-1)/float(k3d-1)
 11   continue
c

      call writearray(y,names(2),idim,jdim,kdim,i3d,j3d,k3d,
     &                dummy, nrec )
      call writearray(z,names(3),idim,jdim,kdim,i3d,j3d,k3d,
     &                dummy, nrec )
      call writearray(x,names(1),idim,jdim,kdim,i3d,j3d,k3d,
     &                dummy, nrec )
c------------------------------------------------- p,t normalised
c     open(11,file='grafdat',form='unformatted')
c     read(11)nr,nz,np
c     read(11)rmax,r0,zmax,zmin
c     read(11)(((vx(j+1,i+1,k+1),i=0,nr),j=0,nz),k=0,np)
c     read(11)(((vy(j+1,i+1,k+1),i=0,nr),j=0,nz),k=0,np)
c     read(11)(((vz(j+1,i+1,k+1),i=0,nr),j=0,nz),k=0,np)
c     read(11)(((t(j+1,i+1,k+1),i=0,nr),j=0,nz),k=0,np)
c     close(11)
c
c-------------------------------------------------
      call writearray(vx,names(4),idim,jdim,kdim,i3d,j3d,k3d,
     &                dummy, nrec )
      call writearray(vy,names(5),idim,jdim,kdim,i3d,j3d,k3d,
     &                dummy, nrec )
      call writearray(vz,names(6),idim,jdim,kdim,i3d,j3d,k3d,
     &                dummy, nrec )
      call writearray(p,names(7),idim,jdim,kdim,i3d,j3d,k3d,
     &                dummy, nrec )
      call writearray(t,names(8),idim,jdim,kdim,i3d,j3d,k3d,
     &                dummy, nrec )
c-------------------------------------------------

      close (unit=91)

      write(*,*) 'vbl1.file  has now been created'

      stop
      end
c-------------------------------------------------
      subroutine writetop(datatype, version, cprog, imax, jmax, kmax,
     &                    icor, i2d, isym, nisurfs, njsurfs, nksurfs,
     &                    ptch, rgas, omega, cp, egam, title,
     &                    bound_type, isurfs, jsurfs, ksurfs, nrec)

      implicit none

      character*4 datatype, cprog
      character*80 title(8)

      integer version, imax, jmax, kmax, icor, i2d, isym , i, j, k, l
      integer nisurfs, njsurfs, nksurfs, nrec
      integer bound_type(3), isurfs(8,nisurfs), jsurfs(8,njsurfs),
     &        ksurfs(8,nksurfs)

      real*4 ptch, rgas, omega, cp, egam

      write(91,rec=1) datatype, version, cprog
      write(91,rec=2) imax, jmax, kmax, icor, i2d, isym, nisurfs,
     &                njsurfs, nksurfs
      write(91,rec=3) ptch, rgas, omega, cp, egam
      do 10 i = 1,8
         write(91,rec=3+i) title(i)
   10 continue
      write(91,rec=12) bound_type
      nrec = 12
      do 20 i = 1,nisurfs
         nrec = nrec + 1
         write(91,rec=nrec) ( isurfs(l,i), l=1,8 )
   20 continue
      do 30 j = 1,njsurfs
         nrec = nrec + 1
         write(91,rec=nrec) ( jsurfs(l,j), l=1,8 )
   30 continue
      do 40 k = 1,nksurfs
         nrec = nrec + 1
         write(91,rec=nrec) ( ksurfs(l,k), l=1,8 )
   40 continue

      return
      end
c-------------------------------------------------
      subroutine  writearray(array, name, idim, jdim, kdim,
     &                       imax, jmax, kmax, dummy, nrec )

      implicit none

      integer idim, jdim, kdim, imax, jmax, kmax, i, j, k, count, ip 
      integer nrec

      character*4 name

      real*8 array(idim,jdim,kdim)
      real*8 dummy(*)

      count = 0
      do 10 k = 1, kmax
      do 10 j = 1, jmax
      do 10 i = 1, imax
         count = count + 1
         dummy(count) = array(i,j,k)
   10 continue

      nrec = nrec + 1
      write(91,rec=nrec) name

      do 20 ip = 1,imax*jmax*kmax,20
         nrec = nrec + 1
         write(91,rec=nrec) ( sngl(dummy(k)), k=ip,ip+19 )
   20 continue

      return
      end
c-------------------------------------------------
