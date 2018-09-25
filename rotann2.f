c----------------------------------------------------------------------c
c     Program ROTANN.F                                                 c
c     ================                                                 c
c     Mathematical modelling of flow in a rotating annulus             c
c     with heat flux applied in the radial direction.                  c
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=128,nzmax=128,mgen=100000)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax),thetan(0:nrmax,0:nzmax)
      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      real*8 omegmax(0:mgen),omegmin(0:mgen),httrans(0:mgen)
      real*8 phimax(0:mgen),phimin(0:mgen)
      character*20 fname1,fname2,fname3,fname4,fname5,fname6,fname7
      character*20 fname,ans
      common/val/cf1,cf2,cf3,cf4,Ra,Re,Ptl
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
      common/sol2/thetan,theta,thetap
      common/num/dr,dz,nr,nz,dt
      common/radius/r0
c
      print*,'filename for vorticity: '
      read*,fname1
      print*,'filename for stream function: '
      read*,fname2
      print*,'filename for axial velocity: '
      read*,fname3
      print*,'filename for temperature: '
      read*,fname4
      print*,'filename for max/min vorticity data (filename or "no"): '
      read*,fname5
      print*,'filename for max/min stream function (filename or "no"): ' 
      read*,fname6
      print*,'filename for heat transfer data (filename or "no"): '
      read*,fname7 
      print*,'time step: '
      read*,dt
      print*,'number of intervals in r-direction: '
      read*,nr
      print*,'number of intervals in z-direction: '
      read*,nz
      print*,'maximum number of time steps: '
      read*,nts
      print*,'r0,rmax,zmax'
      read*,r0,rmax,zmax
      print*,'rotational Rayleigh number: '
      read*,Ra
      print*,'rotational Reynolds number: '
      read*,Re
      print*,'Prandtl number: '
      read*,Ptl
      print*,'temperature perturbation factor'
      read*,pert
c
c----------------------------------------------------------------------c
c     Calculate constants                                              c
c----------------------------------------------------------------------c
c
      pi=dacos(-1.d0)
      dr=(rmax-r0)/dble(nr)
      dz=zmax/dble(nz)
      zmin=0.d0
      cf1=1.d0/(2.d0*Ptl*dt)+1.d0/dr**2+1.d0/dz**2
      cf2=1.d0/(2.d0*Ptl*dt)-1.d0/dr**2-1.d0/dz**2
      cf3=1.d0/(2.d0*dt)+1.d0/dr**2+1.d0/dz**2
      cf4=1.d0/(2.d0*dt)-1.d0/dr**2-1.d0/dz**2
c
c----------------------------------------------------------------------c
c     Set initial conditions                                           c
c----------------------------------------------------------------------c
c
      do 10 i=0,nr
      do 10 j=0,nz
        r=r0+dble(i)*dr
        z=dble(j)*dz
        omeg(i,j)=0.d0
        phi(i,j)=0.d0
c       vel(i,j)=-Ptl*Re*r
        vel(i,j)=0.d0
        vel(0,j)=0.d0
        vel(nr,j)=0.d0
        vel(i,0)=0.d0
        vel(i,nz)=0.d0
        theta(i,j)=(1.d0-dlog(r/rmax)/dlog(r0/rmax))
c    +            +pert*dsin(2.d0*pi*(r/rmax-0.5d0))*dcos(2.d0*pi*z)
c       theta(i,j)=0.5d0
c       theta(0,j)=0.d0
c       theta(nr,j)=1.d0
  10  continue
        theta(nr/4,nz/4)=theta(nr/4,nz/4)+0.01d0
c
      print*,'initial conditions from file?(y/n)'
      read*,ans
      if(ans.eq.'yes'.or.ans.eq.'y')then
        print*,'vorticity filename'
        read*,fname
        open(7,file=fname,form='unformatted')
        read(7)nr,nz
        read(7)rmax,r0,zmax,zmin
        read(7)((omeg(i,j),i=0,nr),j=0,nz)
        close(7)
c
        print*,'strm func. filename'
        read*,fname
        open(8,file=fname,form='unformatted')
        read(8)nr,nz
        read(8)rmax,r0,zmax,zmin
        read(8)((phi(i,j),i=0,nr),j=0,nz)
        close(8)
c
        print*,'velocity filename'
        read*,fname
        open(9,file=fname,form='unformatted')
        read(9)nr,nz
        read(9)rmax,r0,zmax,zmin
        read(9)((vel(i,j),i=0,nr),j=0,nz)
        close(9)
c
        print*,'temperature filename'
        read*,fname
        open(10,file=fname,form='unformatted')
        read(10)nr,nz
        read(10)rmax,r0,zmax,zmin
        read(10)((theta(i,j),i=0,nr),j=0,nz)
        close(10)
      endif
c
      do 15 i=0,nr
      do 15 j=0,nz
        omegp(i,j)=omeg(i,j)
        velp(i,j)=vel(i,j)
        thetap(i,j)=theta(i,j)
  15  continue
        call rufcon2(theta,phi,nr,nz,amax,amin)
c
c----------------------------------------------------------------------c
c     Solve for vorticity, stream function, axial velocity and         c
c     temperature on internal points                                   c
c----------------------------------------------------------------------c
c
      call poissoninit(nr,nz,dr,dz)
c
      do 20 k=1,nts
        call vorticity
        call velocity
        call temperature
        call poisson(phi,omegn)
        call update
        call heat(trans)
c
        if(k/20.eq.dble(k)/20.d0)call rufcon2(theta,phi,nr,nz,amax,amin)
        httrans(k)=-r0*dlog(r0/rmax)*trans/(zmax*rmax)
        print*,httrans(k)
c
c----------------------------------------------------------------------c
c     Convergence check and                                            c
c----------------------------------------------------------------------c
c
        errmax1=0.d0
        errmax2=0.d0
        errmax3=0.d0
        gmax=-10000.d0
        gmin=10000.d0
        pimax=-10000.d0
        pimin=10000.d0
c
        do 30 i=0,nr
        do 30 j=0,nz
          error1=abs(omeg(i,j)-omegp(i,j))
          error2=abs(vel(i,j)-velp(i,j))
          error3=abs(theta(i,j)-thetap(i,j))
          if(error1 .gt. errmax1)errmax1=error1 
          if(error2 .gt. errmax2)errmax2=error2
          if(error3 .gt. errmax3)errmax3=error3
          if(omeg(i,j) .gt. gmax)gmax=omeg(i,j)
          if(omeg(i,j) .lt. gmin)gmin=omeg(i,j)
          if(phi(i,j) .gt. pimax)pimax=phi(i,j)
          if(phi(i,j) .lt. pimin)pimin=phi(i,j)
  30    continue
c
        omegmax(k)=gmax
        omegmin(k)=gmin 
        phimax(k)=pimax
        phimin(k)=pimin
c
        print*,errmax1,errmax2,errmax3,k
c
        if(errmax1.lt.1e-6.and.errmax2.lt.1e-6.and.errmax3.lt.1e-6)then
          print*,'converged'
          goto 40
        endif
  20  continue
c
      print*,'not converged'
c
  40  continue
c
c---------------------------------------------------------------------c
c     Write vorticity, stream function, velocity and temperature data c
c     to file                                                         c
c---------------------------------------------------------------------c
c
      if(fname5 .ne. 'no')then
        open(11,file=fname5,status='unknown')
        write(11,100)(i*dt,omegmax(i),omegmin(i),i=1,k-1)
        close(11)
      endif
c
      if(fname6 .ne. 'no')then
        open(12,file=fname6,status='unknown')
        write(12,100)(i*dt,phimax(i),phimin(i),i=1,k-1)
        close(12)
      endif
c
      if(k .ge. nts)then
        open(13,file='WARNING.WARNING',status='unknown')
        write(13,*)fname1,fname2,fname3,fname4,'USELESS'
        write(13,*)errmax1,errmax2,errmax3,nts
        close(13)
      endif
c
      open(7,file=fname1,form='unformatted')
      write(7)nr,nz
      write(7)rmax,r0,zmax,zmin
      write(7)((omeg(i,j),i=0,nr),j=0,nz)
      close(7)
c
      open(8,file=fname2,form='unformatted')
      write(8)nr,nz
      write(8)rmax,r0,zmax,zmin
      write(8)((phi(i,j),i=0,nr),j=0,nz)
      close(8)
c
      open(9,file=fname3,form='unformatted')
      write(9)nr,nz
      write(9)rmax,r0,zmax,zmin
      write(9)((vel(i,j),i=0,nr),j=0,nz)
      close(9)
c
      open(10,file=fname4,form='unformatted')
      write(10)nr,nz
      write(10)rmax,r0,zmax,zmin
      write(10)((theta(i,j),i=0,nr),j=0,nz)
      close(10)
c
      if(fname7 .ne. 'no')then
        open(14,file=fname7,status='unknown')
        write(14,200)(i*dt,httrans(i),i=1,k-1)
        close(14)
      endif
c
 100  format(3(1x,e16.8))
 200  format(2(1x,e16.8))
c
      stop
      end
c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     Subroutine VORTICITY                                             c
c     ====================                                             c
c     Calculates values of vorticity using the duFort-Frankel method   c
c----------------------------------------------------------------------c
c
      subroutine vorticity
c    
      implicit real*8(a-h,o-z)
      parameter(nrmax=128,nzmax=128)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax),thetan(0:nrmax,0:nzmax)
      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      real*8 arak(0:nrmax,0:nzmax)
      common/val/cf1,cf2,cf3,cf4,Ra,Re,Ptl
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
      common/sol2/thetan,theta,thetap
      common/num/dr,dz,nr,nz,dt
      common/radius/r0
c
      call arakawa(arak,phi,omeg,nr,nz,dr,dz)
      RRP=Ra/Re/Ptl
!mpk      RRP=0.d0
c
      do 10 i=1,nr-1
      do 10 j=1,nz-1
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
c
        r=r0+dble(i)*dr
c
        A=cf1+1.d0/(2.d0*r**2)
        AA=cf2-1.d0/(2.d0*r**2)
        B=(1.d0/(2.d0*Ptl*r*dz))*(phi(i,j)*(omeg(i,jp)-omeg(i,jm))
     +   +omeg(i,j)*(phi(i,jp)-phi(i,jm)))
        C=(1.d0/(Ptl*r*dz))*vel(i,j)*(vel(i,jp)-vel(i,jm))
        D=Re*(vel(i,jp)-vel(i,jm))/dz
!mpk    D=0.d0
        E=r*Ra*(theta(i,jp)-theta(i,jm))/(2.d0*dz)
        EE=RRP*vel(i,j)*(theta(i,jp)-theta(i,jm))/dz
     +    +RRP*theta(i,j)*(vel(i,jp)-vel(i,jm))/dz
        F=(omeg(ip,j)+omeg(im,j))/dr**2
        G=(1.d0/(2.d0*r*dr))*(omeg(ip,j)-omeg(im,j))
        H=(omeg(i,jp)+omeg(i,jm))/dz**2
c
        omegn(i,j)=(AA*omegp(i,j)+arak(i,j)/Ptl+B+C+D-E-EE+F+G+H)/A
c
  10  continue
c
      return
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     subroutine VELOCITY                                              c
c     ===================                                              c
c     Calculates values of velocity using the duFort-Frankel method    c
c----------------------------------------------------------------------c
c
      subroutine velocity
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=128,nzmax=128)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax),thetan(0:nrmax,0:nzmax)
      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      real*8 arak(0:nrmax,0:nzmax)
      common/val/cf1,cf2,cf3,cf4,Ra,Re,Ptl
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
      common/sol2/thetan,theta,thetap
      common/num/dr,dz,nr,nz,dt
      common/radius/r0
c
      call arakawa(arak,phi,vel,nr,nz,dr,dz)
      RRP=Ra/Re/Ptl
!mpk      RRP=0.d0
c
      do 10 i=1,nr-1
      do 10 j=1,nz-1
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
c
        r=r0+dble(i)*dr
c
        A=cf1+1.d0/(2.d0*r**2)
        AA=cf2-1.d0/(2.d0*r**2)
        B=(1.d0/(2.d0*Ptl*r*dz))*(vel(i,j)*(phi(i,jp)-phi(i,jm))
     +   -phi(i,j)*(vel(i,jp)-vel(i,jm)))
        C=Re*(phi(i,jp)-phi(i,jm))/dz
!mpk        C=0.d0
        CC=-RRP*theta(i,j)*(phi(i,jp)-phi(i,jm))/dz
c       CC=0.d0
        D=(vel(ip,j)+vel(im,j))/dr**2
        E=(1.d0/(2.d0*r*dr))*(vel(ip,j)-vel(im,j))
        F=(vel(i,jp)+vel(i,jm))/dz**2
c
        veln(i,j)=(AA*velp(i,j)+arak(i,j)/Ptl-B-C-CC+D+E+F)/A
c
  10  continue
c
      return
      end

c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     subroutine TEMPERATURE                                           c
c     ======================                                           c
c     Calculates temperature using the duFort-Frankel mehod            c
c----------------------------------------------------------------------c
c
      subroutine temperature
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=128,nzmax=128)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax),thetan(0:nrmax,0:nzmax)
      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      real*8 arak(0:nrmax,0:nzmax)
      common/val/cf1,cf2,cf3,cf4,Ra,Re,Ptl
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
      common/sol2/thetan,theta,thetap
      common/num/dr,dz,nr,nz,dt
      common/radius/r0
c
      call arakawa(arak,phi,theta,nr,nz,dr,dz)
c
      do 10 i=1,nr-1
      do 10 j=0,nz
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
c
        r=r0+dble(i)*dr
        if(j.eq.0)then
          jm=jp
          arak(i,j)=0.d0
        endif
        if(j.eq.nz)then
          jp=jm
          arak(i,j)=0.d0
        endif
c
        B=(1.d0/(2.d0*r*dz))*phi(i,j)*(theta(i,jp)-theta(i,jm))
        C=(theta(ip,j)+theta(im,j))/dr**2
        D=(1.d0/(2.d0*r*dr))*(theta(ip,j)-theta(im,j))
        E=(theta(i,jp)+theta(i,jm))/dz**2
c
        thetan(i,j)=(cf4*thetap(i,j)+arak(i,j)+B+C+D+E)/cf3
c
  10  continue
c
      return
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     subroutine ARAKAWA                                               c
c     ==================                                               c
c     computes non-linear terms of the form psi_x th_y - psi_y th_x    c
c----------------------------------------------------------------------c
c
      subroutine arakawa(arak,psi,th,nx,ny,hx,hy)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=128,nzmax=128)
      real*8 arak(0:nrmax,0:nzmax),psi(0:nrmax,0:nzmax)
      real*8 th(0:nrmax,0:nzmax)
c
      do 10 i=1,nx-1
      do 10 j=1,ny-1
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        arak(i,j)=(psi(ip,j)-psi(im,j))*(th(i,j+1)-th(i,j-1))
     +           -(psi(i,j+1)-psi(i,j-1))*(th(ip,j)-th(im,j))
     +           +psi(ip,j)*(th(ip,j+1)-th(ip,j-1))
     +           -psi(im,j)*(th(im,j+1)-th(im,j-1))
     +           -psi(i,j+1)*(th(ip,j+1)-th(im,j+1))
     +           +psi(i,j-1)*(th(ip,j-1)-th(im,j-1))
     +           +th(i,j+1)*(psi(ip,j+1)-psi(im,j+1))
     +           -th(i,j-1)*(psi(ip,j-1)-psi(im,j-1))
     +           -th(ip,j)*(psi(ip,j+1)-psi(ip,j-1))
     +           +th(im,j)*(psi(im,j+1)-psi(im,j-1))
        arak(i,j)=arak(i,j)/(12.d0*hx*hy)
  10  continue
c
      return
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     subroutine HEAT                                                 c
c     ================                                                 c
c     Calculates the heat transfer rate in terms of the Nusselt number c
c     defined as the ratio of heat flux to that due to conduction      c
c----------------------------------------------------------------------c
c
      subroutine heat(trans)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=128,nzmax=128)
      real*8 theta(0:nrmax,0:nzmax),thetan(0:nrmax,0:nzmax)
      real*8 thetap(0:nrmax,0:nzmax)
      common/sol2/thetan,theta,thetap
      common/num/dr,dz,nr,nz,dt
c
      trans=dz*(-25.d0*theta(0,0)+48.d0*theta(1,0)-36.d0*theta(2,0)
     +     +16.d0*theta(3,0)-3.d0*theta(4,0))/(24.d0*dr)
      trans=trans+dz*(-25.d0*theta(0,nz)+48.d0*theta(1,nz)-36.d0
     +     *theta(2,nz)+16.d0*theta(3,nz)-3.d0*theta(4,nz))/(24.d0*dr)
      do 10 j=1,nz-1
        trans=trans+dz*(-25.d0*theta(0,j)+48.d0*theta(1,j)-36.d0
     +       *theta(2,j)+16.d0*theta(3,j)-3.d0*theta(4,j))/(12.d0*dr)
  10  continue
c
      return
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     subroutine UPDATE                                                c
c     =================                                                c
c     updates the values for vorticity, velocity and temperature for   c
c     the next iteration. Boundary values for the three parameters are c
c     also calculated                                                  c
c----------------------------------------------------------------------c
c
      subroutine update
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=128,nzmax=128)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax),thetan(0:nrmax,0:nzmax)
      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      common/val/cf1,cf2,cf3,cf4,Ra,Re,Ptl
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
      common/sol2/thetan,theta,thetap
      common/num/dr,dz,nr,nz,dt
      common/radius/r0
c
      do 10 i=0,nr
        omegn(i,0)=2.d0*phi(i,1)/dz**2
        veln(i,0)=0.d0
        omegn(i,nz)=2.d0*phi(i,nz-1)/dz**2
        veln(i,nz)=0.d0
  10  continue
c
      do 30 j=0,nz
        omegn(0,j)=2.d0*phi(1,j)/dr**2
        veln(0,j)=0.d0
        thetan(0,j)=0.d0
        omegn(nr,j)=2.d0*phi(nr-1,j)/dr**2
        veln(nr,j)=0.d0
        thetan(nr,j)=1.d0
  30  continue
c
      do 40 i=0,nr
      do 40 j=0,nz
        omegp(i,j)=omeg(i,j)
        omeg(i,j)=omegn(i,j)
        velp(i,j)=vel(i,j)
        vel(i,j)=veln(i,j)
        thetap(i,j)=theta(i,j)
        theta(i,j)=thetan(i,j)
  40  continue
c
      return
      end
c
c----------------------------------------------------------------------c
c=====================================================================
      subroutine poissoninit(nr,nz,dr,dz)
      
c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c     ------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      parameter(kr=128,kz=128,mgen=100000)
      real*8 hr(6),hz(6),hhr(6),hhz(6)
      real*8 rval(0:kr)
      integer ibeg(6),n1(6),n2(6)
      common /step/hr,hz,hhr,hhz,rval
      common /ints/ibeg,n1,n2,nlev
      common /radius/r0
c
c
c---------------Determining number of multigrid levels
c
      nlev=0
      itest=nr
      itest2=nr
      jtest=nz
      jtest2=nz
c     print*,itest,itest2,jtest,jtest2
      do 1 i=2,8
        itest=itest/2
        jtest=jtest/2
        if(itest*2-itest2.eq.0.and.jtest*2-jtest2.eq.0)then
          nlev=nlev+1
          itest2=itest
          jtest2=jtest
        else
          goto 2
        end if
    1 continue
    2 continue
      if(nlev.gt.6)nlev=6
        print*,'Number of levels=',nlev
c
      ii=1
      do 20 i=1,nlev
        n1(i)=nr*2/(2**i)
        n2(i)=nz*2/(2**i)
        hr(i)=dr*dble(2**i)/2.d0
        hz(i)=dz*dble(2**i)/2.d0
        hhr(i)=hr(i)**2
        hhz(i)=hz(i)**2
        ibeg(i)=ii
        ii=ii+(n1(i)+1)*(n2(i)+1)+1
        print*,i,n1(i),n2(i),hr(i),hz(i),ibeg(i)
   20 continue
c
      return
      end
c=====================================================================
      subroutine poisson(psi,rhs)

c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c     ------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      parameter(kr=128,kz=128,mgen=100000)
      real*8 psi(0:kr,0:kz),rhs(0:kr,0:kz)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz),psif(0:kr,0:kz),theta(0:kr,0:kz)
      real*8 rc(0:kr,0:kz),rf(0:kr,0:kz),tempa(0:kr,0:kz)
      real*8 hr(6),hz(6),hhr(6),hhz(6),vel(0:kr,0:kz)
      real*8 rval(0:kr)
      integer ibeg(6),n1(6),n2(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hhr,hhz,rval
      common /ints/ibeg,n1,n2,nlev
      common /radius/r0
c
      nvs=10000
      nrel1=2
      nrel2=2
c
      do 10 i=1,mgen
        v(i)=0.d0
        f(i)=0.d0
   10 continue
c
      do 30 i=0,n1(1)
        rval(i)=dble(i)*hr(1) + r0
        do 30 j=0,n2(1)
          psic(i,j)=0.d0
          psif(i,j)=psi(i,j)
          rc(i,j)=0.d0
          rf(i,j)=rhs(i,j)
   30 continue
c
      call trans1(rf,f,ibeg(1),n1(1),n2(1))
      call trans1(psif,v,ibeg(1),n1(1),n2(1))
c
      do 50 i=1,nvs
c        print*,'Sweep ',i
        if(i.eq.1.or.nlev.eq.1)then
          call relax(nrel1,1,nconv)
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser(ilev)
            call relax(nrel1,ilev+1,nconv)
   51     continue
          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax(nrel2,ilev-1,nconv)
   52     continue 
        end if
        if(nconv.eq.1)then
c          print*,'Converged'
          do 60 ii=0,n1(1)
          do 60 jj=0,n2(1)
            psi(ii,jj)=psif(ii,jj)
   60     continue
c         call rufcon2(theta,psi,nr,nz,amax,amin)
          return
        end if
   50 continue
      return
      end
c==================================================
      subroutine relax(nrel,klev,nconv)
      implicit real*8(a-h,o-z)
      parameter(kr=128,kz=128,mgen=100000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz),psif(0:kr,0:kz)
      real*8 rc(0:kr,0:kz),rf(0:kr,0:kz),tempa(0:kr,0:kz)
      real*8 hr(6),hz(6),hhr(6),hhz(6)
      real*8 rval(0:kr)
c
      real*8 tdmaa(kr),tdmab(kr),tdmac(kr),tdmar(kr),tdmas(kr)
c
      integer ibeg(6),n1(6),n2(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hhr,hhz,rval
      common /ints/ibeg,n1,n2,nlev
c
      tolcon=1.d-6
      nconv=0
c
      call trans2(f,rf,ibeg(klev),n1(klev),n2(klev))
      call trans2(v,psif,ibeg(klev),n1(klev),n2(klev))
c
c------------------------------relaxing with some G & S.
c
      drdz=2.d0*(1.d0/hhr(klev)+1.d0/hhz(klev))
      levfact=(2**klev)/2
c
      do 200 irel=1,nrel
c                                 line relax in r-dir
        do 10 j=1,n2(klev)-1
          do 11 i=1,n1(klev)-1
            ii=i*levfact
            tdmaa(i)=1.d0/hhr(klev)-0.5d0/hr(klev)/rval(ii)
            tdmac(i)=1.d0/hhr(klev)+0.5d0/hr(klev)/rval(ii)
            tdmab(i)=-2.d0/hhr(klev)-2.d0/hhz(klev)-1.d0/rval(ii)**2
            tdmar(i)=-(psif(i,j+1)+psif(i,j-1))/hhz(klev)+rf(i,j)
   11     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n1(klev)-1)
          do 12 i=1,n1(klev)-1
            psif(i,j)=tdmas(i)
   12     continue
   10   continue
c                                 line relax in z-dir
        do 13 i=1,n1(klev)-1
          ii=i*levfact
          do 14 j=1,n2(klev)-1
            tdmaa(j)=1.d0/hhz(klev)
            tdmac(j)=1.d0/hhz(klev)
            tdmab(j)=-2.d0/hhr(klev)-2.d0/hhz(klev)-1.d0/rval(ii)**2
            tdmar(j)=-(psif(i+1,j)+psif(i-1,j))/hhr(klev)+rf(i,j)
     +              -(psif(i+1,j)-psif(i-1,j))*0.5d0/rval(ii)/hr(klev)
   14     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n2(klev)-1)
          do 15 j=1,n2(klev)-1
            psif(i,j)=tdmas(j)
   15     continue
   13   continue
c
c                                  Calculating residual.
c
        resmax=0.d0
        do 20 i=1,n1(klev)-1
        do 20 j=1,n2(klev)-1
           ii=i*levfact
          temp=(psif(i+1,j)+psif(i-1,j))/hhr(klev)
     +        +(psif(i+1,j)-psif(i-1,j))*0.5d0/hr(klev)/rval(ii)
     +        -psif(i,j)/rval(ii)**2
     +        +(psif(i,j+1)+psif(i,j-1))/hhz(klev)
     +        -psif(i,j)*drdz
          temp=rf(i,j)-temp
          rc(i,j)=temp
          temp=dabs(temp)
          if(temp.gt.resmax)resmax=temp
   20   continue
c
c       if(klev.eq.1)write(6,100)klev,resmax
  100   format(1x,'Level',i1,'  resmax=',f14.10)
c
  200 continue
c
  210 continue
      call trans1(psif,v,ibeg(klev),n1(klev),n2(klev))
      if(klev.eq.1.and.resmax.lt.tolcon)nconv=1
      return
      end
c==================================================
      subroutine corser(klev)
      implicit real*8(a-h,o-z)
      parameter(kr=128,kz=128,mgen=100000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz),psif(0:kr,0:kz)
      real*8 rc(0:kr,0:kz),rf(0:kr,0:kz),tempa(0:kr,0:kz)
      real*8 hr(6),hz(6),hhr(6),hhz(6)
      real*8 rval(0:kr)
      integer ibeg(6),n1(6),n2(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hhr,hhz,rval
      common /ints/ibeg,n1,n2,nlev
c
c------------------Start off with restriction
c
      call trans3(rc,rf,n1(klev),n2(klev))
c
      do 10 i=1,n1(klev+1)-1
      do 10 j=1,n2(klev+1)-1
        ii=i+i
        jj=j+j
        temp=rf(ii,jj)
        temp=temp*2.d0+rf(ii+1,jj)+rf(ii-1,jj)
        temp=temp+rf(ii,jj+1)+rf(ii,jj-1)
        temp=temp*2.d0+rf(ii+1,jj+1)+rf(ii-1,jj-1)
        temp=temp+rf(ii+1,jj-1)+rf(ii-1,jj+1)
        temp=temp/16.d0
        rc(i,j)=temp
        psic(i,j)=0.d0
   10 continue
      do 20 i=0,n1(klev+1)
        rc(i,0)=0.d0
        rc(i,n1(klev+1))=0.d0
        psic(i,0)=0.d0
        psic(i,n1(klev+1))=0.d0
   20 continue
      do 30 j=0,n1(klev+1)
        rc(0,j)=0.d0
        rc(n1(klev+1),j)=0.d0
        psic(0,j)=0.d0
        psic(n1(klev+1),j)=0.d0
   30 continue
      call trans1(rc,f,ibeg(klev+1),n1(klev+1),n2(klev+1))
      call trans1(psic,v,ibeg(klev+1),n1(klev+1),n2(klev+1))
      return
      end
c==================================================
      subroutine finer(klev)
      implicit real*8(a-h,o-z)
      parameter(kr=128,kz=128,mgen=100000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz),psif(0:kr,0:kz)
      real*8 rc(0:kr,0:kz),rf(0:kr,0:kz),tempa(0:kr,0:kz)
      real*8 hr(6),hz(6),hhr(6),hhz(6)
      real*8 rval(0:kr)
      integer ibeg(6),n1(6),n2(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hhr,hhz,rval
      common /ints/ibeg,n1,n2,nlev
c
      call trans2(v,psic,ibeg(klev),n1(klev),n2(klev))
c
c--------------------prolongation
c
      do 10 i=0,n1(klev-1)
      do 10 j=0,n2(klev-1)
        i1=i/2
        j1=j/2
        i2=i1*2
        j2=j1*2
        if(i2.eq.i)then
          if(j2.eq.j)then
            psif(i,j)=psic(i1,j1)
          else
            psif(i,j)=(psic(i1,j1)+psic(i1,j1+1))*0.5d0
          end if
        else
          if(j2.eq.j)then
            psif(i,j)=(psic(i1,j1)+psic(i1+1,j1))*0.5d0
          else
            temp=psic(i1,j1)+psic(i1+1,j1+1)
            temp=temp+psic(i1,j1+1)+psic(i1+1,j1)
            psif(i,j)=temp*0.25d0
          end if
        end if
   10 continue
c
      call trans2(v,tempa,ibeg(klev-1),n1(klev-1),n2(klev-1))
c
      do 20 i=0,n1(klev-1)
      do 20 j=0,n2(klev-1)
        psif(i,j)=psif(i,j)+tempa(i,j)
        psic(i,j)=0.d0
   20 continue
      call trans1(psic,v,ibeg(klev),n1(klev),n2(klev))
      call trans1(psif,v,ibeg(klev-1),n1(klev-1),n2(klev-1))
      return
      end
c==================================================
      subroutine trans1(a2d,a1d,i1,n1,n2)
      implicit real*8(a-h,o-z)
      parameter(kr=128,kz=128,mgen=100000)
      real*8 a2d(0:kr,0:kz),a1d(mgen)
c
      do 10 j=0,n2
      do 10 i=0,n1
        ii=j*(n1+1) + i + i1
        a1d(ii)=a2d(i,j)
   10 continue
      return
      end
c==================================================
      subroutine trans2(a1d,a2d,i1,n1,n2)
      implicit real*8(a-h,o-z)
      parameter(kr=128,kz=128,mgen=100000)
      real*8 a2d(0:kr,0:kz),a1d(mgen)
c
      do 10 j=0,n2
      do 10 i=0,n1
        ii=j*(n1+1) + i + i1
        a2d(i,j)=a1d(ii)
   10 continue
      return
      end
c==================================================
      subroutine trans3(a2d1,a2d2,n1,n2)
      implicit real*8(a-h,o-z)
      parameter(kr=128,kz=128,mgen=100000)
      real*8 a2d1(0:kr,0:kz),a2d2(0:kr,0:kz)
c
      do 10 j=0,n2
      do 10 i=0,n1
        a2d2(i,j)=a2d1(i,j)
   10 continue
      return
      end
c==================================================
      subroutine rufcon2(a,b,nx,ny,amax,amin)
      parameter(mx=128,my=128,mgen=100000)
      implicit real*8(a-h,o-z)
      real*8 a(0:mx,0:my)
      real*8 b(0:mx,0:my)
      character*1 c(0:401)
      character*1 ch(11)
      ch(1)='0'
      ch(2)='1'
      ch(3)='2'
      ch(4)='3'
      ch(5)='4'
      ch(6)='5'
      ch(7)='6'
      ch(8)='7'
      ch(9)='8'
      ch(10)='9'
      ch(11)='t'
      c(nx+1)=' '
      nxy=(nx+1)*(ny+1)
      amax=woppa(a,nx,ny)
      amin=titch(a,nx,ny)
      bmax=woppa(b,nx,ny)
      bmin=titch(b,nx,ny)
      print*,'arraymax,min=',amax,amin
      print*,'arraymax,min=',bmax,bmin
      adiff=amax-amin
      bdiff=bmax-bmin
      if(adiff.ne.0.)adiff=1./adiff
      if(bdiff.ne.0.)bdiff=1./bdiff
      do 30 j=0,ny,2
      jj=ny-j
      do 10 i=0,nx,2
      xa=1.5+10.*(a(i,jj)-amin)*adiff
      xb=1.5+10.*(b(i,jj)-bmin)*bdiff
      na=xa
      nb=xb
      c(i)=ch(na)
      c(i+nx+2)=ch(nb)
  10  continue
      m2=nx*2+2
      m1=0
      if(nx.gt.130)m1=nx*2+2-130
      write(6,20)(c(k),k=m1,m2,2)
  20  format(1x,200a1)
  30  continue
      return
      end
c==================================================
      double precision function titch(a,nr,nz)
      parameter(kr=128,kz=128,mgen=100000)
      implicit real*8(a-h,o-z)
      real*8 a(0:kr,0:kz)
c
      z=a(0,0)
      do 10 i=0,nr
      do 10 j=0,nz
      if(z.gt.a(i,j))z=a(i,j)
  10  continue
      titch=z
      return
      end
c==================================================
      double precision function woppa(a,nr,nz)
      parameter(kr=128,kz=128,mgen=100000)
      implicit real*8(a-h,o-z)
      real*8 a(0:kr,0:kz)
c
      z=a(0,0)
      do 10 i=0,nr
      do 10 j=0,nz
      if(z.lt.a(i,j))z=a(i,j)
  10  continue
      woppa=z
      return
      end
c==================================================
      subroutine tdma(a,b,c,y,r,n)
      implicit real*8(a-h,o-z)
      real*8 a(n),b(n),c(n),y(n),r(n)
c
      do 10 j=2,n
        fact=a(j)/b(j-1)
        b(j)=b(j)-fact*c(j-1)
        r(j)=r(j)-fact*r(j-1)
   10 continue
c
c
      y(n)=r(n)/b(n)
c
c
      do 20 j=n-1,1,-1
        y(j)=(r(j)-c(j)*y(j+1))/b(j)
   20 continue
      return
      end

c================================================















 
