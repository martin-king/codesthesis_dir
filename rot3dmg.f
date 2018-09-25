c23456789012345678901234567890123456789012345678901234567890123456789012
c     Program rot3dmg.f                                                c
c     ====================                                             c
c!mpk located disk14/enptwl/proj14.dir/rotstag14.f                     c
c     coded originally by T.Lewis(1999).                               c
c     Mathematical modelling of flow in a rotating annulus             c
c     with heat flux applied in the radial direction.                  c
c     Three dimensional velocity/vorticity model on staggered grid.    c
c     mpk.retains original non-d scalings by Lewis.                    c
c        .introduces gravitational terms(present only in Avort and     c
c         Bvort. A Grashof number(Gr) entails.                         c
c        .hybrid scheme introduced in ABCVort and Temp(commented out)  c
c        .corrections                                                  c
c        .u, v, w located at walls on bound                            c
c         kr=65,kz=65,kp=128,mgen=626953                               c
c        extrapolation and SOR added in subroutine POISSON?            c     
c         LAST UPDATED: 19-3-01.                                       c
c         emergency measure. tolcon in poisson? changed fr 1d-6 to 1d-5c
c         10-1-01.    compiled to rot3dmg.rub                          c                                       
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      real*8 urh(0:kr,0:kz,0:kp),vrh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
      real*8 httrans(0:mgen)
!mpk      real*8 ssumc(6,mgen),ssums(6,mgen)
      character*20 fname1,fname2,fname3,fname4,fname5
      character*20 fname6,fname7,fname8,fname9,ans
!mpk      common/grals/ssumc,ssums
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
      common/iteration/kt
c
      print*,'filename for A vorticity: '
      read*,fname1
      print*,'filename for B vorticity: '
      read*,fname2
      print*,'filename for C vorticity: '
      read*,fname3
      print*,'filename for u velocity: '
      read*,fname4                 
      print*,'filename for v velocity: '
      read*,fname5
      print*,'filename for w velocity: '
      read*,fname6
      print*,'filename for temperature: '
      read*,fname7
      print*,'filename for heat transfer data (filename or "no"): '
      read*,fname8 
      print*,'time step: '
      read*,dt
      print*,'number of intervals in r-direction: '
      read*,nr
      print*,'number of intervals in z-direction: '
      read*,nz
      print*,'number of intervals in phi-direction: '
      read*,np
      print*,'maximum number of time steps: '
      read*,nts
      print*,'r0,rmax,zmax'
      read*,r0,rmax,zmax
      print*,'rotational Rayleigh number: '
      read*,Ra
      print*,'rotational Reynolds number: '
      read*,Re
      print*,'Grashof number(=0.d0 if g terms are excluded: '
      read*,Gr
      print*,'Prandtl number: '
      read*,Ptl
      print*,'axial temperature perturbation factor'
      read*,pert
      print*,'circumferential amplitude factor'
      read*,amp
      print*,'nval'
      read*,nval 
c
c----------------------------------------------------------------------c
c     Calculate constants                                              c
c----------------------------------------------------------------------c
c
      zmin=0.d0
      pi=dacos(-1.d0)
      dr=(rmax-r0)/dble(nr)
      dz=zmax/dble(nz)
      dp=2.d0*pi/dble(np)
c
c----------------------------------------------------------------------c
c     Set initial conditions                                           c
c----------------------------------------------------------------------c
c
      do 10 i=0,nr+1
      do 10 j=0,nz+1
      do 10 k=0,np
        if(i.eq.0)then
          r=r0
        elseif(i.eq.nr+1)then
          r=rmax
        else
          r=r0+(dble(i)-0.5d0)*dr
        endif
!mpk        z=dble(j)*dz
        An(i,j,k)=0.d0
        Bn(i,j,k)=0.d0
        Cn(i,j,k)=0.d0
        A(i,j,k)=0.d0
        B(i,j,k)=0.d0
        C(i,j,k)=0.d0
        uel(i,j,k)=0.d0
c       vel(i,j,k)=-Ptl*Re*(r-r0)*(rmax-r)
c       vel(i,j,k)=-Ptl*Re*r
        vel(i,j,k)=0.d0
!mpk        vel(0,j,k)=0.d0
!mpk        vel(i,0,k)=0.d0
!mpk        vel(i,nz+1,k)=0.d0
!mpk        vel(nr+1,j,k)=0.d0
        wel(i,j,k)=0.d0
        phi=8.d0*datan(1.d0)*dble(k)/dble(np)
        th(i,j,k)=1.d0-dlog(r)/dlog(r0)
c    +      +pert*dsin(2.d0*pi*(r/rmax-0.5d0))*dcos(2.d0*pi*z)
     +      +amp*(rmax-r)*(r-r0)*dcos(dble(nval)*phi)
c    +      +amp*dcos(pi*z/zmax)*(rmax-r)*(r-r0)*dcos(dble(nval)*phi)
  10  continue
!mpk      do 11 j=1,nz-1
!mpk      th(nr/4,j,np/2)=th(nr/4,j,np/2)-0.2d0
!mpk  11  continue
          th(nr/4,nz/4,np/2)=th(nr/4,nz/4,np/2)-0.5d0*th(nr/4,nz/4,np/2)
c
!mpk      do 12 i=1,nr
!mpk      do 12 j=1,nz
!mpk      do 12 k=0,np
!mpk        r=r0+dble(i)*dr
!mpk        A(i,j,k)=-(vel(i,j+1,k)-vel(i,j,k))/dz
!mpk        A(i,0,k)=-dum10(vel,dz,i,0,k)
!mpk        A(i,nz,k)=-dum12(vel,dz,i,nz+1,k)
!mpk        C(i,j,k)=(vel(i+1,j,k)-vel(i,j,k))/dr
!mpk     +           +0.5d0*(vel(i+1,j,k)+vel(i,j,k))/r
!mpk        C(0,j,k)=dum9(vel,dr,0,j,k)
!mpk        C(nr,j,k)=dum11(vel,dr,nr+1,j,k)+vel(nr+1,j,k)/rmax
!mpk        A(i,0,k)=-2.d0*vel(i,1,k)/dz
!mpk        A(i,nz,k)=2.d0*vel(i,nz,k)/dz
!mpk        C(0,j,k)=2.d0*vel(1,j,k)/dr
!mpk        C(nr,j,k)=-2.d0*vel(nr,j,k)/dr
!mpk  12  continue
c
!mpk      call rufcon(th,nr+1,nz+1,np)
      call poissoninit(nr,nz,np,dr,dz,dp,nlev)
c
      print*,'initial conditions from file (yes/no):'
      read*,ans
      if(ans .ne. 'no')then
        print*,'filename for vorticity A: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)nr,nz,np
        read(11)rmax,r0,zmax,zmin
        read(11)(((A(i,j,k),i=0,nr+1),j=0,nz),k=0,np)
        close(11)
c
        print*,'filename for vorticity B: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)nr,nz,np
        read(11)rmax,r0,zmax,zmin
        read(11)(((B(i,j,k),i=0,nr),j=0,nz),k=0,np)
        close(11)
c
        print*,'filename for vorticity C: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)nr,nz,np
        read(11)rmax,r0,zmax,zmin
        read(11)(((C(i,j,k),i=0,nr),j=0,nz+1),k=0,np)
        close(11)
c
        print*,'filename for velocity u: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)nr,nz,np
        read(11)rmax,r0,zmax,zmin
        read(11)(((uel(i,j,k),i=0,nr),j=0,nz+1),k=0,np)
        close(11)
c
        print*,'filename for velocity v: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)nr,nz,np
        read(11)rmax,r0,zmax,zmin
        read(11)(((vel(i,j,k),i=0,nr+1),j=0,nz+1),k=0,np)
        close(11)
c
        print*,'filename for velocity w: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)nr,nz,np
        read(11)rmax,r0,zmax,zmin
        read(11)(((wel(i,j,k),i=0,nr+1),j=0,nz),k=0,np)
        close(11)
c
        print*,'filename for temperature: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)nr,nz,np
        read(11)rmax,r0,zmax,zmin
        read(11)(((th(i,j,k),i=0,nr+1),j=0,nz+1),k=0,np)
        close(11)
c
      endif
!mpk      call rufcon(th,nr+1,nz+1,np)
c
      do 50 i=0,nr+1
      do 50 j=0,nz+1
      do 50 k=0,np
        Ap(i,j,k)=A(i,j,k)
        Bp(i,j,k)=B(i,j,k)
        Cp(i,j,k)=C(i,j,k)
        thp(i,j,k)=th(i,j,k)
        uelp(i,j,k)=uel(i,j,k)
        velp(i,j,k)=vel(i,j,k)
        welp(i,j,k)=wel(i,j,k)
  50  continue
c
c----------------------------------------------------------------------c
c     Solve for vorticity, stream function, axial velocity and         c
c     temperature on internal points                                   c
c----------------------------------------------------------------------c
c
      ntot=0
c
      do 15 ku=1,9999
      do 20 kt=1,nts
        call Avorticity
        call Bvorticity
        call Cvorticity
        call temperature
        call heat(trans1)
!mpk        call aaaaaaaa(ntot+kt)
        call rhsw(wrh)
        call poisson1(wel,welp,wrh,nr,nz,np,dr,dz,dp,nlev)
        call rhsv(vrh)
        call poisson2(vel,velp,vrh,nr,nz,np,dr,dz,dp,nlev)
        call rhsu(urh)
        call poisson3(uel,uelp,urh,nr,nz,np,dr,dz,dp,nlev)
        call update
c
        httrans(ntot+kt)=-r0*dlog(r0/rmax)*trans1/(zmax*2.d0*pi)
        print*,httrans(ntot+kt),ntot+kt
c
!mpk        call rufcon(th,nr+1,nz+1,np)
  20  continue
c
      ntot=ntot+nts
c
      print*,'more iterations? (y/n): '
      read*,ans
      if(ans.eq.'n'.or.ans.eq.'no')goto 30
      print*,'no. of iterations: '
      read*,nts
c
  15  continue
c
  30  continue
c
c---------------------------------------------------------------------c
c     Write vorticity, stream function, velocity and temperature data c
c     to file                                                         c
c---------------------------------------------------------------------c
c
      open(11,file=fname1,form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((A(i,j,k),i=0,nr+1),j=0,nz),k=0,np)
      close(11)
c
      open(11,file=fname2,form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((B(i,j,k),i=0,nr),j=0,nz),k=0,np)
      close(11)
c
      open(11,file=fname3,form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((C(i,j,k),i=0,nr),j=0,nz+1),k=0,np)
      close(11)
c
      open(11,file=fname4,form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((uel(i,j,k),i=0,nr),j=0,nz+1),k=0,np)
      close(11)
c
      open(11,file=fname5,form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((vel(i,j,k),i=0,nr+1),j=0,nz+1),k=0,np)
      close(11)
c
      open(11,file=fname6,form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((wel(i,j,k),i=0,nr+1),j=0,nz),k=0,np)
      close(11)
c
      open(11,file=fname7,form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((th(i,j,k),i=0,nr+1),j=0,nz+1),k=0,np)
      close(11)
c
      if(fname8 .ne. 'no')then
        open(11,file=fname8,status='unknown')
        write(11,200)((i-1)*dt,httrans(i),i=1,ntot)
!mpk                                      ssumc(1,i),ssums(1,i),
!mpk     +  ssumc(2,i),ssums(2,i),ssumc(3,i),ssums(3,i),ssumc(4,i),
!mpk     +  ssums(4,i),ssumc(5,i),ssums(5,i),ssumc(6,i),ssums(6,i),i=1,ntot)
        close(11)
      endif
c
 200  format(2(1x,e16.8))
c
      stop
      end
c
c----------------------------------------------------------------------c
c     Subroutine AVORTICITY                                            c
c     =====================                                            c
c     Calculates values of vorticity using the duFort-Frankel method   c
c----------------------------------------------------------------------c
c
      subroutine Avorticity
c    
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      do 10 i=1,nr
      do 10 j=1,nz-1
      do 10 k=0,np
c
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        km=k-1
        ka=k+1
        if(k.eq.0)km=np-1
        if(k.eq.np)ka=1
c
        r=r0+(dble(i)-0.5d0)*dr
        u=(uel(i,j,k)+uel(i,jp,k)+uel(i,j,ka)+uel(i,jp,ka)
     +   +uel(im,j,k)+uel(im,jp,k)+uel(im,j,ka)+uel(im,jp,ka))*0.125d0
        v=(vel(i,jp,k)+vel(i,j,k))*0.5d0
        w=(wel(i,j,ka)+wel(i,j,k))*0.5d0
        tem=(th(i,j,k)+th(i,j,ka)+th(i,jp,k)+th(i,jp,ka))*0.25d0
        BB=(B(i,j,k)+B(im,j,k)+B(i,j,ka)+B(im,j,ka))*0.25d0
        CC=(C(i,j,k)+C(i,jp,k)+C(im,j,k)+C(im,jp,k))*0.25d0
c
        if(i.eq.1)then
         Arra=(3.2d0*A(i-1,j,k)+2.d0*A(i+1,j,k)-0.2d0*A(i+2,j,k))/dr**2
         Ara=(-4.d0*A(i-1,j,k)/3.d0+A(i+1,j,k)/3.d0)/dr
         Ar=(-4.d0*A(i-1,j,k)/3.d0+A(i,j,k)+A(i+1,j,k)/3.d0)/dr
         cnr=-2.5d0/dr**2+0.5d0/dr/r
        elseif(i.eq.nr)then
         Arra=(3.2d0*A(i+1,j,k)+2.d0*A(i-1,j,k)-0.2d0*A(i-2,j,k))/dr**2
         Ara=(-A(i-1,j,k)/3.d0+4.d0*A(i+1,j,k)/3.d0)/dr
         Ar=(-4.d0*A(i-1,j,k)/3.d0+A(i,j,k)+A(i+1,j,k)/3.d0)/dr
         cnr=-2.5d0/dr**2-0.5d0/dr/r
        else
         Arra=(A(i-1,j,k)+A(i+1,j,k))/dr**2
         Ara=(A(i+1,j,k)-A(i-1,j,k))*0.5d0/dr
         Ar=Ara
         cnr=-1.d0/dr**2
c
c-------------------hybrid scheme introduced by mpk--------------------c
         Pe=dabs(u*dr/Ptl)
         if(Pe.gt.2)then        
!mpk          print*,'Pe_uA=',Pe
          if(u.gt.0)then
             Ar=(A(i,j,k)-A(i-1,j,k))/dr
          elseif(u.lt.0)then
             Ar=(A(i+1,j,k)-A(i,j,k))/dr
          endif
         endif
c-----------------hybrid scheme ends for Ar in advectice term----------c
c
        endif
c
        uelr=((uel(i,j,k)+uel(i,jp,k)+uel(i,j,ka)+uel(i,jp,ka))*0.25d0
     +      -(uel(im,j,k)+uel(im,jp,k)+uel(im,j,ka)+uel(im,jp,ka))
     +      *0.25d0)/dr
       Azza=(A(i,j-1,k)+A(i,j+1,k))/dz**2
       Az=(A(i,j+1,k)-A(i,j-1,k))*0.5d0/dz
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(w*dz/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_wA=',Pe
         if(w.gt.0)then
           Az=(A(i,j,k)-A(i,j-1,k))/dz
         elseif(w.lt.0)then
           Az=(A(i,j+1,k)-A(i,j,k))/dz
         endif
        endif
c-----------------hybrid scheme ends for Az in advectice term----------c
c
        uelz=((uel(i,jp,k)+uel(i,jp,ka)+uel(im,jp,k)+uel(im,jp,ka))*
     +       0.25d0-(uel(i,j,k)+uel(i,j,ka)+uel(im,j,k)+uel(im,j,ka))*
     +       0.25d0)/dz
        thz=(th(i,jp,k)+th(i,jp,ka)-th(i,j,k)-th(i,j,ka))*0.5d0/dz
        thh=(th(i,j,ka)+th(i,jp,ka)-th(i,j,k)-th(i,jp,k))*0.5d0/dp
c
        Ahha=(A(i,j,ka)+A(i,j,km))/dp**2
        Ah=(A(i,j,ka)-A(i,j,km))*0.5d0/dp
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(v*dp*r/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_vA=',Pe
         if(v.gt.0)then
           Ah=(A(i,j,k)-A(i,j,km))/dp
         elseif(v.lt.0)then
           Ah=(A(i,j,ka)-A(i,j,k))/dp
         endif
        endif
c-----------------hybrid scheme ends for Ah in advectice term----------c
c
        Bh=(B(i,j,ka)+B(im,j,ka)-B(i,j,k)-B(im,j,k))*0.5d0/dp
        uelh=((uel(i,j,ka)+uel(i,jp,ka)+uel(im,j,ka)+uel(im,jp,ka))*
     +       0.25d0-(uel(i,j,k)+uel(i,jp,k)+uel(im,j,k)+uel(im,jp,k))*
     +       0.25d0)/dp
c
        cnz=-1.d0/dz**2
        cnh=-1.d0/dp**2/r**2        
c
        RHS1=Arra+Ara/r+Ahha/r**2+Azza-2.d0*Bh/r**2+A(i,j,k)/(Ptl*dt)
        ULT=(u*Ar+v*Ah/r+w*Az-uelh*BB/r-uelz*CC)/Ptl
        RHS2=-2.d0*Re*uelz+2.d0*Ra*(tem*uelz+thz*u)/Re/Ptl
     +       -Gr*Ptl/r*thh
        cf1=1.d0/(Ptl*dt)-uelr/(2.d0*Ptl)+1.d0/(2.d0*r**2)-cnr-cnz-cnh
        cf2=uelr/(2.d0*Ptl)-1.d0/(2.d0*r**2)+cnr+cnz+cnh
        An(i,j,k)=(cf2*Ap(i,j,k)-ULT-RHS2+RHS1)/cf1
c
  10  continue
c
      return
      end
c
c----------------------------------------------------------------------c
c     Subroutine BVORTICITY                                            c
c     =====================                                            c
c     Calculates values of vorticity using the duFort-Frankel method   c
c----------------------------------------------------------------------c
c
      subroutine Bvorticity
c   
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      do 10 i=1,nr-1
      do 10 j=1,nz-1
      do 10 k=0,np
c
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        km=k-1
        ka=k+1
        if(k.eq.0)km=np-1
        if(k.eq.np)ka=1
c
        r=r0+dble(i)*dr
        u=(uel(i,j,k)+uel(i,jp,k))*0.5d0
        v=(vel(i,j,km)+vel(i,jp,km)+vel(i,j,k)+vel(i,jp,k)
     +   +vel(ip,j,km)+vel(ip,jp,km)+vel(ip,j,k)+vel(ip,jp,k))*0.125d0
        w=(wel(i,j,k)+wel(ip,j,k))*0.5d0
        tem=(th(i,jp,k)+th(ip,jp,k)+th(i,j,k)+th(ip,j,k))*0.25d0
        AA=(A(i,j,k)+A(ip,j,k)+A(i,j,km)+A(ip,j,km))*0.25d0
        CC=(C(i,j,k)+C(i,jp,k)+C(i,j,km)+C(i,jp,km))*0.25d0
c
        Brra=(B(i-1,j,k)+B(i+1,j,k))/dr**2
        Bra=(B(i+1,j,k)-B(i-1,j,k))*0.5d0/dr
        Br=Bra
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(u*dr/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_uB=',Pe
         if(u.gt.0)then
           Br=(B(i,j,k)-B(i-1,j,k))/dr
         elseif(u.lt.0)then
           Br=(B(i+1,j,k)-B(i,j,k))/dr
         endif
        endif
c-----------------hybrid scheme ends for Br in advectice term----------c
c
        velr=((vel(ip,j,km)+vel(ip,jp,km)+vel(ip,j,k)+vel(ip,jp,k))*
     +      0.25d0-(vel(i,j,km)+vel(i,jp,km)+vel(i,j,k)+vel(i,jp,k))*
     +      0.25d0)/dr
        Bzza=(B(i,j-1,k)+B(i,j+1,k))/dz**2
        Bz=(B(i,j+1,k)-B(i,j-1,k))*0.5d0/dz
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(w*dz/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_wB=',Pe
         if(w.gt.0)then
           Bz=(B(i,j,k)-B(i,j-1,k))/dz
         elseif(w.lt.0)then
           Bz=(B(i,j+1,k)-B(i,j,k))/dz
         endif
        endif
c-----------------hybrid scheme ends for Bz in advectice term----------c
c
        velz=((vel(i,jp,km)+vel(ip,jp,km)+vel(i,jp,k)+vel(ip,jp,k))*
     +       0.25d0-(vel(i,j,km)+vel(ip,j,km)+vel(i,j,k)+vel(ip,j,k))*
     +       0.25d0)/dz
        thz=(th(i,jp,k)+th(ip,jp,k)-th(i,j,k)-th(ip,j,k))*0.5d0/dz
        thr=(th(ip,jp,k)+th(ip,j,k)-th(i,jp,k)-th(i,j,k))*0.5d0/dr
        Bhha=(B(i,j,ka)+B(i,j,km))/dp**2
        Bh=(B(i,j,ka)-B(i,j,km))*0.5d0/dp
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(v*dp*r/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_vB=',Pe
         if(v.gt.0)then
           Bh=(B(i,j,k)-B(i,j,km))/dp
         elseif(v.lt.0)then
           Bh=(B(i,j,ka)-B(i,j,k))/dp
         endif
        endif
c-----------------hybrid scheme ends for Bh in advectice term----------c
c
        Ah=(A(i,j,k)+A(ip,j,k)-A(i,j,km)-A(ip,j,km))*0.5d0/dp
        velh=((vel(i,j,k)+vel(ip,j,k)+vel(i,jp,k)+vel(ip,jp,k))*
     +       0.25d0-(vel(i,j,km)+vel(ip,j,km)+vel(i,jp,km)+vel(ip,jp,km)
     +       )*0.25d0)/dp
c
        cnr=-1.d0/dr**2
        cnz=-1.d0/dz**2
        cnh=-1.d0/dp**2/r**2
c
        RHS1=Brra+Bra/r+Bhha/r**2+Bzza+2.d0*Ah/r**2+B(i,j,k)/(Ptl*dt)
        ULT=(u*Br+v*Bh/r+w*Bz-(velr-v/r)*AA-velz*CC)/Ptl
        RHS2=-2.d0*Re*velz+r*Ra*thz+2.d0*Ra*(tem*velz+thz*v)/Re/Ptl
     +       +Gr*Ptl*thr
        cf1=1.d0/(Ptl*dt)-(u/r+velh/r)/(2.d0*Ptl)
     +     +1.d0/(2.d0*r**2)-cnr-cnz-cnh
        cf2=(u/r+velh/r)/(2.d0*Ptl)
     +     -1.d0/(2.d0*r**2)+cnr+cnz+cnh
        Bn(i,j,k)=(cf2*Bp(i,j,k)-ULT-RHS2+RHS1)/cf1
c
  10  continue
c
      return
      end
c
c----------------------------------------------------------------------c
c     Subroutine CVORTICITY                                            c
c     =====================                                            c
c     Calculates values of vorticity using the duFort-Frankel method   c
c----------------------------------------------------------------------c
c
      subroutine Cvorticity
c   
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      do 10 i=1,nr-1
      do 10 j=1,nz
      do 10 k=0,np
c
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        km=k-1
        ka=k+1
        if(k.eq.0)km=np-1
        if(k.eq.np)ka=1
c
        r=r0+dble(i)*dr
        u=(uel(i,j,ka)+uel(i,j,k))*0.5d0
        v=(vel(i,j,k)+vel(ip,j,k))*0.5d0
        w=(wel(i,j,k)+wel(ip,j,k)+wel(i,j,ka)+wel(ip,j,ka)+
     +   wel(i,jm,k)+wel(ip,jm,k)+wel(i,jm,ka)+wel(ip,jm,ka))*0.125d0
        AA=(A(i,j,k)+A(i,jm,k)+A(ip,j,k)+A(ip,jm,k))*0.25d0
        BB=(B(i,j,k)+B(i,jm,k)+B(i,j,ka)+B(i,jm,ka))*0.25d0
        tem=(th(i,j,k)+th(ip,j,k)+th(ip,j,ka)+th(i,j,ka))*0.25d0
c
        Crra=(C(i-1,j,k)+C(i+1,j,k))/dr**2
        Cra=(C(i+1,j,k)-C(i-1,j,k))*0.5d0/dr
        Cr=Cra
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(u*dr/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_uC=',Pe
         if(u.gt.0)then
           Cr=(C(i,j,k)-C(i-1,j,k))/dr
         elseif(u.lt.0)then
           Cr=(C(i+1,j,k)-C(i,j,k))/dr
         endif
        endif
c-----------------hybrid scheme ends for Cr in advectice term----------c
c
        welr=((wel(ip,j,k)+wel(ip,j,ka)+wel(ip,jm,k)+wel(ip,jm,ka))*
     +       0.25d0-(wel(i,j,k)+wel(i,j,ka)+wel(i,jm,k)+wel(i,jm,ka))*
     +       0.25d0)/dr
        uelr=(uel(ip,j,k)+uel(ip,j,ka)-uel(im,j,k)-uel(im,j,ka))*
     +       0.25d0/dr
        thr=(th(ip,j,k)+th(ip,j,ka)-th(i,j,k)-th(i,j,ka))*0.5d0/dr
        if(j.eq.1)then
         Czza=(3.2d0*C(i,j-1,k)+2.d0*C(i,j+1,k)-0.2d0*C(i,j+2,k))/dz**2
         Cz=(-4.d0*C(i,j-1,k)/3.d0+C(i,j,k)+C(i,j+1,k)/3.d0)/dz
          cnz=-2.5d0/dz**2
        elseif(j.eq.nz)then
         Czza=(3.2d0*C(i,j+1,k)+2.d0*C(i,j-1,k)-0.2d0*C(i,j-2,k))/dz**2
         Cz=(-C(i,j-1,k)/3.d0-C(i,j,k)+4.d0*C(i,j+1,k)/3.d0)/dz  
          cnz=-2.5d0/dz**2
        else
         Czza=(C(i,j-1,k)+C(i,j+1,k))/dz**2
         Cz=(C(i,j+1,k)-C(i,j-1,k))*0.5d0/dz
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(w*dz/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_wC=',Pe
         if(w.gt.0)then
           Cz=(C(i,j,k)-C(i,j-1,k))/dz
         elseif(w.lt.0)then
           Cz=(C(i,j+1,k)-C(i,j,k))/dz
         endif
        endif
c-----------------hybrid scheme ends for Cz in advectice term----------c
c
          cnz=-1.d0/dz**2
        endif
        welz=((wel(i,j,k)+wel(ip,j,k)+wel(i,j,ka)+wel(ip,j,ka))*0.25d0-
     +       (wel(i,jm,k)+wel(ip,jm,k)+wel(i,jm,ka)+wel(ip,jm,ka))*
     +       0.25d0)/dz
        Chha=(C(i,j,ka)+C(i,j,km))/dp**2
        Ch=(C(i,j,ka)-C(i,j,km))*0.5d0/dp
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(v*dp*r/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_vC=',Pe
         if(v.gt.0)then
           Ch=(C(i,j,k)-C(i,j,km))/dp
         elseif(v.lt.0)then
           Ch=(C(i,j,ka)-C(i,j,k))/dp
         endif
        endif
c-----------------hybrid scheme ends for Ch in advectice term----------c
c
        welh=((wel(i,j,ka)+wel(ip,j,ka)+wel(i,jm,ka)+wel(ip,jm,ka))*
     +       0.25d0-(wel(i,j,k)+wel(ip,j,k)+wel(i,jm,k)+wel(ip,jm,k))*
     +       0.25d0)/dp
        velh=(vel(i,j,ka)+vel(ip,j,ka)-vel(i,j,km)-vel(ip,j,km))*
     +       0.25d0/dp
        thh=(th(i,j,ka)+th(ip,j,ka)-th(i,j,k)-th(ip,j,k))*0.5d0/dp
c
        cnr=-1.d0/dr**2
        cnh=-1.d0/dp**2/r**2
c
        RHS1=Crra+Cra/r+Chha/r**2+Czza+C(i,j,k)/(Ptl*dt)
        ULT=(u*Cr+v*Ch/r+w*Cz-welr*AA-welh*BB/r)/Ptl
        RHS2=-2.d0*Re*welz-Ra*thh+2.d0*Ra*(tem*welz-thr*u-thh*v/r)
     +       /Re/Ptl
        cf1=1.d0/(Ptl*dt)-welz/(2.d0*Ptl)-cnr-cnz-cnh
        cf2=welz/(2.d0*Ptl)+cnr+cnz+cnh
        Cn(i,j,k)=(cf2*Cp(i,j,k)-ULT-RHS2+RHS1)/cf1
c
  10  continue
c
      return
      end
c
c----------------------------------------------------------------------c
c     subroutine TEMPERATURE                                           c
c     ======================                                           c
c     Calculates temperature using the duFort-Frankel mehod            c
c----------------------------------------------------------------------c
c
      subroutine temperature
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      do 10 i=1,nr
      do 10 j=0,nz+1
      do 10 k=0,np
c
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        km=k-1
        ka=k+1
        if(k.eq.0)km=np-1
        if(k.eq.np)ka=1
c
        r=r0+(dble(i)-0.5d0)*dr
        if((j.eq.0).or.(j.eq.nz+1))then
          u=0.d0
          v=0.d0
          w=0.d0
        else
          u=(uel(i,j,k)+uel(im,j,k))*0.5d0
          v=(vel(i,j,k)+vel(i,j,km))*0.5d0
          w=(wel(i,j,k)+wel(i,jm,k))*0.5d0
        endif
c
        if(i.eq.1)then
         thrra=(3.2d0*th(i-1,j,k)+2.d0*th(i+1,j,k)-0.2d0*th(i+2,j,k))/
     +         dr**2
         thra=(-4.d0*th(i-1,j,k)/3.d0+th(i+1,j,k)/3.d0)/dr
         thr=(-4.d0*th(i-1,j,k)/3.d0+th(i,j,k)+th(i+1,j,k)/3.d0)/dr
          cnr=-2.5d0/dr**2+0.5d0/dr/r
        elseif(i.eq.nr)then
          thrra=(3.2d0*th(i+1,j,k)+2.d0*th(i-1,j,k)-0.2d0*th(i-2,j,k))/
     +          dr**2
          thra=(-th(i-1,j,k)/3.d0+4.d0*th(i+1,j,k)/3.d0)/dr
          thr=(-th(i-1,j,k)/3.d0-th(i,j,k)+4.d0*th(i+1,j,k)/3.d0)/dr
          cnr=-2.5d0/dr**2-0.5d0/dr/r
        else
          thrra=(th(i-1,j,k)+th(i+1,j,k))/dr**2
          thra=(th(i+1,j,k)-th(i-1,j,k))*0.5d0/dr
          thr=thra
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(u*dr/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_uth=',Pe
         if(u.gt.0)then
           thr=(th(i,j,k)-th(i-1,j,k))/dr
         elseif(u.lt.0)then
           thr=(th(i+1,j,k)-th(i,j,k))/dr
         endif
        endif
c-----------------hybrid scheme ends for thr in advectice term---------c
c
          cnr=-1.d0/dr**2
        endif
c
        if(j.eq.1)then
          thzza=(3.2d0*th(i,j-1,k)+2.d0*th(i,j+1,k)-0.2d0*th(i,j+2,k))/
     +          dz**2
          thz=(-4.d0*th(i,j-1,k)/3.d0+th(i,j,k)+th(i,j+1,k)/3.d0)/dz
          cnz=-2.5d0/dz**2
        elseif(j.eq.nz)then
          thzza=(3.2d0*th(i,j+1,k)+2.d0*th(i,j-1,k)-0.2d0*th(i,j-2,k))/
     +          dz**2
          thz=(-th(i,j-1,k)/3.d0-th(i,j,k)+4.d0*th(i,j+1,k)/3.d0)/dz
          cnz=-2.5d0/dz**2
        elseif(j.eq.0)then
          thzza=8.d0*th(i,j+1,k)/(dz**2)
          thz=0.d0
          cnz=-4.d0/dz**2
        elseif(j.eq.nz+1)then
          thzza=8.d0*th(i,j-1,k)/(dz**2)
          thz=0.d0
          cnz=-4.d0/dz**2
        else
          thzza=(th(i,j-1,k)+th(i,j+1,k))/dz**2
          thz=(th(i,j+1,k)-th(i,j-1,k))*0.5d0/dz
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(w*dz/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_uth=',Pe
         if(w.gt.0)then
           thz=(th(i,j,k)-th(i,j-1,k))/dz
         elseif(w.lt.0)then
           thz=(th(i,j+1,k)-th(i,j,k))/dz
         endif
        endif
c-----------------hybrid scheme ends for thz in advectice term---------c
c
          cnz=-1.d0/dz**2
        endif
c
        thhha=(th(i,j,ka)+th(i,j,km))/dp**2
        thh=(th(i,j,ka)-th(i,j,km))*0.5d0/dp
c
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(v*dp*r/Ptl)
        if(Pe.gt.2)then
!mpk         print*,'Pe_uth=',Pe
         if(v.gt.0)then
           thh=(th(i,j,k)-th(i,j,km))/dp
         elseif(v.lt.0)then
           thh=(th(i,j,ka)-th(i,j,k))/dp
         endif
        endif
c-----------------hybrid scheme ends for thh in advectice term---------c
c
        cnp=-1.d0/dp**2/r**2
c
        RHS1=thrra+thra/r+thhha/r**2+thzza+th(i,j,k)/dt  
        ULT=u*thr+v*thh/r+w*thz
        cf1=1.d0/(dt)-cnr-cnz-cnp
        cf2=cnr+cnz+cnp
        thn(i,j,k)=(cf2*thp(i,j,k)-ULT+RHS1)/cf1
c
  10  continue
c
      return
      end
c
c----------------------------------------------------------------------c
c     subroutine HEAT                                                  c
c     ================                                                 c
c     Calculates the heat transfer rate in terms of the Nusselt number c
c     defined as the ratio of heat flux to that due to conduction      c
c----------------------------------------------------------------------c
c
      subroutine heat(trans1)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      common/sol2/thn,th,thp
      common/num/dt,dr,dz,dp,nr,nz,np
c
      trans1=0.d0
      do 10 k=1,np
      trans=dz*(-184.d0*th(0,0,k)+225.d0*th(1,0,k)
     +     -50.d0*th(2,0,k)+9.d0*th(3,0,k))/(240.d0*dr)
      trans=trans+dz*(-184.d0*th(0,1,k)+225.d0*th(1,1,k)
     +     -50.d0*th(2,1,k)+9.d0*th(3,1,k))/(80.d0*dr)
      trans=trans+dz*(-184.d0*th(0,nz,k)+225.d0*th(1,nz,k)
     +     -50.d0*th(2,nz,k)+9.d0*th(3,nz,k))/(80.d0*dr)
      trans=trans+dz*(-184.d0*th(0,nz+1,k)+225.d0*th(1,nz+1,k)
     +     -50.d0*th(2,nz+1,k)+9.d0*th(3,nz+1,k))/(240.d0*dr)
c
      do 20 j=2,nz-1
        trans=trans+dz*(-184.d0*th(0,j,k)+225.d0*th(1,j,k)-50.d0
     +       *th(2,j,k)+9.d0*th(3,j,k))/(60.d0*dr)
  20  continue
        trans1=trans1+trans*dp
  10  continue
c
      return
      end
c
!mpkc----------------------------------------------------------------------c
!mpkc     subroutine INT                                                   c
!mpkc     ==============                                                   c
!mpkc     Calculates the volume integrals                                  c
!mpkc----------------------------------------------------------------------c
!mpkc
!mpk      subroutine aaaaaaaa(kt)
!mpkc
!mpk      implicit real*8(a-h,o-z)
!mpk      parameter(kr=33,kz=33,kp=160,mgen=100000)
!mpk      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
!mpk      real*8 ssumc(6,mgen),ssums(6,mgen)
!mpk      common/grals/ssumc,ssums
!mpk      common/sol2/thn,th,thp
!mpk      common/num/dt,dr,dz,dp,nr,nz,np
!mpkc
!mpk      do 10 n=1,6
!mpk        ssumc(n,kt)=0.d0
!mpk        ssums(n,kt)=0.d0 
!mpk        sumc=0.d0
!mpk        sums=0.d0
!mpk        do 20 k=1,np
!mpk          ang=dble(k)*dp
!mpk          do 30 i=0,nr+1
!mpk            r=r0+dble(i)*dr
!mpk          do 30 j=0,nz+1
!mpk            fact=1.d0
!mpk            if(i.eq.0.or.i.eq.nr+1)fact=fact/2.d0
!mpk            if(j.eq.0.or.j.eq.nz+1)fact=fact/2.d0
!mpk            fc=r*th(i,j,k)*dcos(dble(n)*ang)
!mpk            fs=r*th(i,j,k)*dsin(dble(n)*ang)
!mpk            sumc=sumc+fc*fact
!mpk            sums=sums+fs*fact
!mpk   30     continue
!mpk          sumc=sumc*dr*dz
!mpk          sums=sums*dr*dz
!mpkc
!mpk          ssumc(n,kt)=ssumc(n,kt)+sumc*dp
!mpk          ssums(n,kt)=ssums(n,kt)+sums*dp
!mpk   20   continue
!mpk   10 continue
!mpkc
!mpk      return
!mpk      end
c
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
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c------left and right side walls (omitting edges)
      do 10 i=1,nr
      do 10 k=0,np
         An(i,0,k)=-(3.d0*vel(i,1,k)-vel(i,2,k)/3.d0)/dz
         An(i,nz,k)=-(vel(i,nz-1,k)/3.d0-3.d0*vel(i,nz,k))/dz
  10  continue
      do 11 i=1,nr-1
      do 11 k=0,np
        Bn(i,0,k)=(3.d0*uel(i,1,k)-uel(i,2,k)/3.d0)/dz
        Bn(i,nz,k)=(-3.d0*uel(i,nz,k)+uel(i,nz-1,k)/3.d0)/dz
  11  continue
c------inner and outer cylindrical walls (omitting edges)
      do 12 j=1,nz-1
      do 12 k=0,np
        Bn(0,j,k)=-(3.d0*wel(1,j,k)-wel(2,j,k)/3.d0)/dr
        Bn(nr,j,k)=-(wel(nr-1,j,k)/3.d0-3.d0*wel(nr,j,k))/dr
  12  continue 
      do 13 j=1,nz
      do 13 k=0,np
        Cn(0,j,k)=(3.d0*vel(1,j,k)-vel(2,j,k)/3.d0)/dr
        Cn(nr,j,k)=(vel(nr-1,j,k)/3.d0-3.d0*vel(nr,j,k))/dr
  13  continue
c
      do 30 j=0,nz+1
      do 30 k=0,np
        thn(0,j,k)=0.d0
        thn(nr+1,j,k)=1.d0
  30  continue
c
      do 40 i=0,nr+1
      do 40 j=0,nz+1
      do 40 k=0,np
        Ap(i,j,k)=A(i,j,k)
        A(i,j,k)=An(i,j,k)
        Bp(i,j,k)=B(i,j,k)
        B(i,j,k)=Bn(i,j,k)
        Cp(i,j,k)=C(i,j,k)
        C(i,j,k)=Cn(i,j,k)
        thp(i,j,k)=th(i,j,k)
        th(i,j,k)=thn(i,j,k)
  40  continue
c
      return
      end
c
c----------------------------------------------------------------------c
      subroutine rhsw(wrh)
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      real*8 wrh(0:kr,0:kz,0:kp)
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      do 10 i=1,nr
      do 10 j=1,nz-1
      do 10 k=0,np
        km=k-1
        if(k.eq.0)km=np-1
        r=r0+(dble(i)-0.5d0)*dr
        Ah=(An(i,j,k)-An(i,j,km))/dp
        Br=(Bn(i,j,k)-Bn(i-1,j,k))/dr
        if(i.eq.1)then
          wrh(i,j,k)=Ah/r-Bn(1,j,k)/dr-0.5d0*Bn(1,j,k)/r
        elseif(i.eq.nr)then
          wrh(i,j,k)=Ah/r+Bn(nr-1,j,k)/dr-0.5d0*Bn(nr-1,j,k)/r
        else
          wrh(i,j,k)=Ah/r-Br-0.5d0*(Bn(i,j,k)+Bn(i-1,j,k))/r
        endif
   10 continue
c
      return 
      end
c
c----------------------------------------------------------------------c
      subroutine rhsv(vrh)
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      real*8 vrh(0:kr,0:kz,0:kp)
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      do 10 i=1,nr
      do 10 j=1,nz
      do 10 k=0,np
        r=r0+(dble(i)-0.5d0)*dr
        Az=(An(i,j,k)-An(i,j-1,k))/dz
        Cr=(Cn(i,j,k)-Cn(i-1,j,k))/dr
        if(i.eq.1)then
          vrh(i,j,k)=Cn(1,j,k)/dr+Cn(1,j,k)/r
        elseif(i.eq.nr)then
          vrh(i,j,k)=-Cn(nr-1,j,k)/dr+Cn(nr-1,j,k)/r
        else
          vrh(i,j,k)=Cr+(Cn(i,j,k)+Cn(i-1,j,k))/r
        endif
        if(j.eq.1)then
          vrh(i,j,k)=vrh(i,j,k)-An(i,1,k)/dz
        elseif(j.eq.nz)then
          vrh(i,j,k)=vrh(i,j,k)+An(i,nz-1,k)/dz
        else
          vrh(i,j,k)=vrh(i,j,k)-Az
        endif
   10 continue
c
      return
      end
c
c----------------------------------------------------------------------c
      subroutine rhsu(urh)
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=100000)
      real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
      real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
      real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
      real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
      real*8 uelp(0:kr,0:kz,0:kp),velp(0:kr,0:kz,0:kp)
      real*8 welp(0:kr,0:kz,0:kp)
      real*8 urh(0:kr,0:kz,0:kp)
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      do 10 i=1,nr-1
      do 10 j=1,nz
      do 10 k=0,np
        km=k-1
        if(k.eq.0)km=np-1
        r=r0+dble(i)*dr
        Bz=(Bn(i,j,k)-Bn(i,j-1,k))/dz
        welz=0.5d0*(wel(i,j,k)+wel(i+1,j,k)-wel(i,j-1,k)
     +      -wel(i+1,j-1,k))/dz
        Ch=(Cn(i,j,k)-Cn(i,j,km))/dp
        if(j.eq.1)then
          urh(i,j,k)=Bn(i,1,k)/dz-Ch/r-2.d0*welz/r
        elseif(j.eq.nz)then
          urh(i,j,k)=-Bn(i,nz-1,k)/dz-Ch/r-2.d0*welz/r
        else
          urh(i,j,k)=Bz-Ch/r-2.d0*welz/r
        endif
   10 continue
c
      return
      end
c
c----------------------------------------------------------------------c
      subroutine poissoninit(nr,nz,np,dr,dz,dp,nlev)
c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
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
          if(itest2.ge.4.and.jtest2.ge.4)nlev=nlev+1
          itest2=itest
          jtest2=jtest
        else
          goto 2
        end if
    1 continue
    2 continue
      print*,'number of MG levels=',nlev
      if(nlev.gt.4)nlev=4
c
      ii=1
      do 20 i=1,nlev
        n1(i)=nr*2/(2**i)
        n2(i)=nz*2/(2**i)
        n3(i)=np*2/(2**i)
        hr(i)=dr*dble(2**i)*0.5d0
        hz(i)=dz*dble(2**i)*0.5d0
        hp(i)=dp*dble(2**i)*0.5d0
        hhr(i)=hr(i)**2
        hhz(i)=hz(i)**2
        hhp(i)=hp(i)**2
        ibeg(i)=ii
        ii=ii+(n1(i)+1)*(n2(i)+1)*(n3(i)+1)+1
c       print*,i
c       print*,n1(i),n2(i),n3(i),ibeg(i)
c       print*,hr(i),hz(i),hp(i)
   20 continue
c
      return
      end
c
c----------------------------------------------------------------------c
      subroutine poisson1(psi,psip,rhs,nr,nz,np,dr,dz,dp,nlev)
c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 psi(0:kr,0:kz,0:kp),psip(0:kr,0:kz,0:kp)
      real*8 rhs(0:kr,0:kz,0:kp)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
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
      do 30 i=0,nr+1
      do 30 j=0,nz
      do 30 k=0,np
        psic(i,j,k)=0.d0
        psif(i,j,k)=2.d0*psi(i,j,k)-psip(i,j,k)   !to improve initial guess
!mpk        psif(i,j,k)=psi(i,j,k)
        rc(i,j,k)=0.d0
        rf(i,j,k)=rhs(i,j,k)
!mpk        psip(i,j,k)=psi(i,j,k)
   30 continue
c
      call trans1(rf,f,ibeg(1),n1(1),n2(1),n3(1))
      call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
c
      do 50 i=1,nvs
        if(i.gt.500)stop
c       print*,'Sweep ',i
        if(i.ne.1)then
        do 31 is=0,nr+1
        do 31 js=0,nz
        do 31 ks=0,np
        psif(is,js,ks)=psip(is,js,ks)+1.2d0*(psif(is,js,ks)-
     +    psip(is,js,ks))        
        psip(is,js,ks)=psif(is,js,ks)
   31   continue
        call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
        endif
        if(i.eq.1.or.nlev.eq.1)then
          call relax1(nrel1,1,nconv,nlev)
            do 54 ii=0,nr+1
            do 54 jj=0,nz
            do 54 kk=0,np
             psip(ii,jj,kk)=psif(ii,jj,kk)
   54       continue
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser(ilev)
            call relax1(nrel1,ilev+1,nconv,nlev)
   51     continue
          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax1(nrel2,ilev-1,nconv,nlev)
   52     continue
        end if
        if(nconv.eq.1)then
c         print*,'Converged'
          do 63 ii=0,nr+1
          do 63 jj=0,nz
          do 63 kk=0,np
            psip(ii,jj,kk)=psi(ii,jj,kk)
   63     continue
c         UPDATING
c         internal points 
          do 60 ii=1,nr
          do 60 jj=1,nz-1
          do 60 kk=0,np
            psi(ii,jj,kk)=psif(ii,jj,kk)
   60     continue
c     endwalls
          do 61 ii=0,nr+1
          do 61 kk=0,np
            psi(ii,0,kk)=0.d0
            psi(ii,nz,kk)=0.d0
   61     continue
c     cylinders
          do 62 jj=0,nz
          do 62 kk=0,np
              psi(0,jj,kk)=0.d0
              psi(nr+1,jj,kk)=0.d0
   62     continue
          print*,'w vel; sweeps=',i
c         call rufcon(psi,nr,nz,np)
          return
        end if
   50 continue
      return
      end
c
c----------------------------------------------------------------------c
      subroutine poisson2(psi,psip,rhs,nr,nz,np,dr,dz,dp,nlev)
c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 psi(0:kr,0:kz,0:kp),psip(0:kr,0:kz,0:kp)
      real*8 rhs(0:kr,0:kz,0:kp)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
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
      do 30 i=0,nr+1
      do 30 j=0,nz+1
      do 30 k=0,np
        psic(i,j,k)=0.d0
        psif(i,j,k)=2.d0*psi(i,j,k)-psip(i,j,k)
!mpk        psif(i,j,k)=psi(i,j,k)
        rc(i,j,k)=0.d0
        rf(i,j,k)=rhs(i,j,k)
!mpk        psip(i,j,k)=psi(i,j,k)
   30 continue
c
      call trans1(rf,f,ibeg(1),n1(1),n2(1),n3(1))
      call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
c
      do 50 i=1,nvs
        if(i.gt.500)stop
c       print*,'Sweep ',i
      if(i.ne.1)then
        do 31 is=0,nr+1
        do 31 js=0,nz+1
        do 31 ks=0,np
        psif(is,js,ks)=psip(is,js,ks)+1.2d0*(psif(is,js,ks)-
     +   psip(is,js,ks))
        psip(is,js,ks)=psif(is,js,ks)
   31   continue
        call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
      endif
        if(i.eq.1.or.nlev.eq.1)then
          call relax2(nrel1,1,nconv,nlev)
          do 55 ii=0,nr+1
            do 55 jj=0,nz+1
            do 55 kk=0,np
              psip(ii,jj,kk)=psif(ii,jj,kk)
   55     continue
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser(ilev)
            call relax2(nrel1,ilev+1,nconv,nlev)
   51     continue
          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax2(nrel2,ilev-1,nconv,nlev)
   52     continue
        end if
        if(nconv.eq.1)then
c         print*,'Converged'
          do 64 ii=0,nr+1
          do 64 jj=0,nz+1
          do 64 kk=0,np
            psip(ii,jj,kk)=psi(ii,jj,kk)
   64     continue
c       UPDATING
c       internal points
          do 60 ii=1,nr
          do 60 jj=1,nz
          do 60 kk=0,np
            psi(ii,jj,kk)=psif(ii,jj,kk)
   60     continue
c       cylinders
          do 62 jj=0,nz+1
          do 62 kk=0,np
            psi(0,jj,kk)=0.d0        
            psi(nr+1,jj,kk)=0.d0     
   62     continue
c       endwalls
          do 63 ii=1,nr
          do 63 kk=0,np
            psi(ii,0,kk)=0.d0        
            psi(ii,nz+1,kk)=0.d0
   63     continue    
          print*,'v vel; sweeps=',i
c         call rufcon(psi,nr,nz,np)
          return
        end if
   50 continue
      return
      end
c
c----------------------------------------------------------------------c
      subroutine poisson3(psi,psip,rhs,nr,nz,np,dr,dz,dp,nlev)
c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 psi(0:kr,0:kz,0:kp),psip(0:kr,0:kz,0:kp)
      real*8 rhs(0:kr,0:kz,0:kp)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
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
      do 30 i=0,nr
      do 30 j=0,nz+1
      do 30 k=0,np
        psic(i,j,k)=0.d0
        psif(i,j,k)=2.d0*psi(i,j,k)-psip(i,j,k)
!mpk        psif(i,j,k)=psi(i,j,k)
        rc(i,j,k)=0.d0
        rf(i,j,k)=rhs(i,j,k)
!mpk        psip(i,j,k)=psi(i,j,k)
   30 continue
c
      call trans1(rf,f,ibeg(1),n1(1),n2(1),n3(1))
      call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
c
      do 50 i=1,nvs
        if(i.gt.500)stop
c       print*,'Sweep ',i
       if(i.ne.1)then
        do 31 is=0,nr
        do 31 js=0,nz+1
        do 31 ks=0,np
        psif(is,js,ks)=psip(is,js,ks)+1.2d0*(psif(is,js,ks)-
     +    psip(is,js,ks))
        psip(is,js,ks)=psif(is,js,ks)
   31   continue
        call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
       endif
        if(i.eq.1.or.nlev.eq.1)then
          call relax3(nrel1,1,nconv,nlev)
            do 54 ii=0,nr
            do 54 jj=0,nz+1
            do 54 kk=0,np
             psip(ii,jj,kk)=psif(ii,jj,kk)
   54       continue
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser(ilev)
            call relax3(nrel1,ilev+1,nconv,nlev)
   51     continue
          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax3(nrel2,ilev-1,nconv,nlev)
   52     continue
        end if
        if(nconv.eq.1)then
c         print*,'Converged'
          do 63 ii=0,nr
          do 63 jj=0,nz+1
          do 63 kk=0,np
            psip(ii,jj,kk)=psi(ii,jj,kk)
   63     continue
c       UPDATING
c       internal points
          do 60 ii=1,nr-1
          do 60 jj=1,nz
          do 60 kk=0,np
            psi(ii,jj,kk)=psif(ii,jj,kk)
   60     continue
c       cylinders
          do 61 jj=0,nz+1
          do 61 kk=0,np
            psi(0,jj,kk)=0.d0
            psi(nr,jj,kk)=0.d0
   61     continue
c       endwalls
          do 62 ii=0,nr
          do 62 kk=0,np
             psi(ii,0,kk)=0.d0
             psi(ii,nz+1,kk)=0.d0
   62     continue
          print*,'u vel; sweeps=',i
c         call rufcon(psi,nr,nz,np)
          return
        end if
   50 continue
      return
      end
c
c----------------------------------------------------------------------c
      subroutine relax1(nrel,klev,nconv,nlev)
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      real*8 pdmaa(kp),pdmab(kp),pdmac(kp),pdmad(kp),pdmae(kp),pdmas(kp)
      real*8 pdmar(kp)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
c
      tolcon=1.d-5
      nconv=0
c
      call trans2(f,rf,ibeg(klev),n1(klev),n2(klev),n3(klev))
      call trans2(v,psif,ibeg(klev),n1(klev),n2(klev),n3(klev))
c
c------------------------------relaxing with some G & S.
c
      do 200 irel=1,nrel
c                                 line relax in r-dir
        do 10 k=1,n3(klev)
        do 10 j=1,n2(klev)-1
          ka=k+1
          km=k-1
          if(k.eq.n3(klev))ka=1
          do 11 i=1,n1(klev)
            rval=r0+(dble(i)-0.5d0)*hr(klev)
            if(i.eq.1)then
              pdmaa(i)=0.d0
              pdmab(i)=0.d0
              pdmac(i)=-2.d0/hhr(klev)-1.d0/hr(klev)/rval*0.5d0
     +                -2.d0/hhz(klev)-2.d0/hhp(klev)/rval**2
              pdmad(i)=5.d0/hhr(klev)/3.d0+1.d0/hr(klev)/rval*0.5d0
              pdmae(i)=-0.2d0/hhr(klev)
            elseif(i.eq.n1(klev))then
              pdmaa(i)=-0.2d0/hhr(klev)
              pdmab(i)=5.d0/hhr(klev)/3.d0-1.d0/hr(klev)/rval*0.5d0
              pdmac(i)=-2.d0/hhr(klev)+1.d0/hr(klev)/rval*0.5d0
     +                -2.d0/hhz(klev)-2.d0/hhp(klev)/rval**2
              pdmad(i)=0.d0
              pdmae(i)=0.d0
            else
              pdmaa(i)=0.d0
              pdmab(i)=1.d0/hhr(klev)-0.5d0/hr(klev)/rval
              pdmac(i)=-(2.d0/hhr(klev)+2.d0/hhz(klev)+2.d0/hhp(klev)
     +              /rval**2)
              pdmad(i)=1.d0/hhr(klev)+0.5d0/hr(klev)/rval
              pdmae(i)=0.d0
            endif
            pdmar(i)=-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)
     +              -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +              +rf(i,j,k)
   11     continue
          call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n1(klev))
          do 12 i=1,n1(klev)
            psif(i,j,k)=pdmas(i)
            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
   12     continue
   10   continue
c                               line relax in z-dir
        do 13 k=1,n3(klev)
        do 13 i=1,n1(klev)
          rval=r0+(dble(i)-0.5d0)*hr(klev)
          ka=k+1
          km=k-1
          if(k.eq.n3(klev))ka=1
          do 14 j=1,n2(klev)-1
            tdmaa(j)=1.d0/hhz(klev)
            tdmab(j)=-(2.d0/hhz(klev)+2.d0/hhp(klev)/rval**2)
            tdmac(j)=1.d0/hhz(klev)
            if(i.eq.1)then
              tdmab(j)=tdmab(j)-2.d0/hhr(klev)-1.d0/hr(klev)
     +                /rval*0.5d0
              tdmar(j)=-(5.d0*psif(i+1,j,k)/3.d0
     +                -0.2*psif(i+2,j,k))/hhr(klev)
     +                -1.d0*psif(i+1,j,k)/hr(klev)/rval*0.5d0
     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            elseif(i.eq.n1(klev))then
              tdmab(j)=tdmab(j)-2.d0/hhr(klev)+1.d0/hr(klev)
     +                /rval*0.5d0
              tdmar(j)=-(5.d0*psif(i-1,j,k)/3.d0
     +                -0.2*psif(i-2,j,k))/hhr(klev)
     +                +1.d0*psif(i-1,j,k)/hr(klev)/rval*0.5d0
     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            else
              tdmab(j)=tdmab(j)-2.d0/hhr(klev)
              tdmar(j)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +            -0.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)/rval
     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            endif
   14     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n2(klev)-1)
          do 15 j=1,n2(klev)-1
            psif(i,j,k)=tdmas(j)
            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
   15     continue
   13   continue
c                             line relax in phi-dir
        do 16 j=1,n2(klev)-1
        do 16 i=1,n1(klev)
          rval=r0+(dble(i)-0.5d0)*hr(klev)
          do 17 k=1,n3(klev)
            tdmaa(k)=1.d0/hhp(klev)/rval**2
            tdmab(k)=-(2.d0/hhz(klev)+2.d0/hhp(klev)/rval**2)
            tdmac(k)=1.d0/hhp(klev)/rval**2
            if(i.eq.1)then
              tdmab(k)=tdmab(k)-2.d0/hhr(klev)-1.d0/hr(klev)
     +                /rval*0.5d0
              tdmar(k)=-(5.d0*psif(i+1,j,k)/3.d0
     +                -0.2*psif(i+2,j,k))/hhr(klev)
     +                -psif(i+1,j,k)/hr(klev)/rval*0.5d0
     +                -(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)
     +                +rf(i,j,k)
            elseif(i.eq.n1(klev))then
              tdmab(k)=tdmab(k)-2.d0/hhr(klev)+1.d0/hr(klev)
     +                /rval*0.5d0
              tdmar(k)=-(5.d0*psif(i-1,j,k)/3.d0
     +                -0.2*psif(i-2,j,k))/hhr(klev)
     +                +psif(i-1,j,k)/hr(klev)/rval*0.5d0
     +                -(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)
     +                +rf(i,j,k)
            else
              tdmab(k)=tdmab(k)-2.d0/hhr(klev)
              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +          -0.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)/rval
     +              -(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)+rf(i,j,k)
            endif
   17     continue
          call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n3(klev))
          do 18 k=0,n3(klev)
            if(k.eq.0)then
              psif(i,j,k)=tdmas(n3(klev))
            else
              psif(i,j,k)=tdmas(k)
            endif
   18     continue
   16   continue
c
c                                  Calculating residual.
c
          resmax=0.d0
          do 20 i=1,n1(klev)
            rval=r0+(dble(i)-0.5d0)*hr(klev)
          do 20 j=1,n2(klev)-1
          do 20 k=1,n3(klev)
          ka=k+1
          km=k-1
          if(k.eq.n3(klev))ka=1
          temp=(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)
     +        +(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +        -psif(i,j,k)*(2.d0/hhz(klev)+2.d0/hhp(klev)/rval**2)
          if(i.eq.1)then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0*psif(i+1,j,k)/3.d0
     +          -0.2d0*psif(i+2,j,k))/hhr(klev)
     +         +0.5d0*(-psif(i,j,k)+psif(i+1,j,k))/hr(klev)/rval
          elseif(i.eq.n1(klev))then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0*psif(i-1,j,k)/3.d0
     +          -0.2d0*psif(i-2,j,k))/hhr(klev)+0.5d0*(psif(i,j,k)
     +          -psif(i-1,j,k))/hr(klev)/rval
          else
            temp=temp+(psif(i-1,j,k)-2.d0*psif(i,j,k)+psif(i+1,j,k))
     +          /hhr(klev)+0.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +          /rval
          endif
            temp=rf(i,j,k)-temp
            rc(i,j,k)=temp
            temp=dabs(temp)
            if(temp.gt.resmax)resmax=temp
   20     continue
c
c         if(klev.eq.1)write(6,100)klev,resmax
c         write(6,100)klev,resmax
c       end if
  100   format(1x,'Level',i1,'  resmax=',f14.10)
c
  200 continue
c
      call trans1(psif,v,ibeg(klev),n1(klev),n2(klev),n3(klev))
      if(klev.eq.1.and.resmax.lt.tolcon)nconv=1
c
      return
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine relax2(nrel,klev,nconv,nlev)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      real*8 pdmaa(kp),pdmab(kp),pdmac(kp),pdmad(kp),pdmae(kp),pdmas(kp)
      real*8 pdmar(kp)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
c
      tolcon=1.d-5
      nconv=0
c
      call trans2(f,rf,ibeg(klev),n1(klev),n2(klev),n3(klev))
      call trans2(v,psif,ibeg(klev),n1(klev),n2(klev),n3(klev))
c
c------------------------------relaxing with some G & S.
c
      do 200 irel=1,nrel
c                                 line relax in r-dir
        do 10 k=1,n3(klev)
        do 10 j=1,n2(klev)
          ka=k+1
          km=k-1
          if(k.eq.0)km=n3(klev)-1
          if(k.eq.n3(klev))ka=1
            do 11 i=1,n1(klev)
             rval=r0+(dble(i)-0.5d0)*hr(klev)            
            if(i.eq.1)then
              pdmaa(i)=0.d0
              pdmab(i)=0.d0
              pdmac(i)=-2.d0/hhr(klev)-2.d0/hhp(klev)/rval**2
     +                 +1.d0/rval**2       !mpk +1.d0/rval/hr(klev)
              pdmad(i)=5.d0/3.d0/hhr(klev)+4.d0/3.d0/hr(klev)/rval
              pdmae(i)=-0.2d0/hhr(klev)
            elseif(i.eq.n1(klev))then
              pdmaa(i)=-0.2d0/hhr(klev)
              pdmab(i)=5.d0/3.d0/hhr(klev)-4.d0/3.d0/hr(klev)/rval
              pdmac(i)=-2.d0/hhr(klev)-2.d0/hhp(klev)/rval**2
     +                +1.d0/rval**2        !mpk -1.d0/rval/hr(klev)
              pdmad(i)=0.d0
              pdmae(i)=0.d0
            else
              pdmaa(i)=0.d0
              pdmab(i)=1.d0/hhr(klev)-1.5d0/hr(klev)/rval
              pdmac(i)=-(2.d0/hhr(klev)+2.d0/hhp(klev)/rval**2)
     +                +1.d0/rval**2
              pdmad(i)=1.d0/hhr(klev)+1.5d0/hr(klev)/rval
              pdmae(i)=0.d0
            endif
c
            if(j.eq.1)then
             pdmac(i)=pdmac(i)-2.d0/hhz(klev)
             pdmar(i)=-(5.d0/3.d0*psif(i,j+1,k)-0.2*psif(i,j+2,k))
     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            elseif(j.eq.n2(klev))then
             pdmac(i)=pdmac(i)-2.d0/hhz(klev)
             pdmar(i)=-(5.d0/3.d0*psif(i,j-1,k)-0.2*psif(i,j-2,k))
     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            else
             pdmac(i)=pdmac(i)-2.d0/hhz(klev)
             pdmar(i)=-(psif(i,j-1,k)+psif(i,j+1,k))
     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            endif
c
   11     continue
          call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n1(klev))
          do 12 i=1,n1(klev)
            psif(i,j,k)=pdmas(i)
            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
   12     continue
   10   continue
c                                 line relax in z-dir
        do 13 k=1,n3(klev)
        do 13 i=1,n1(klev)
            rval=r0+(dble(i)-0.5d0)*hr(klev)
          ka=k+1
          km=k-1
          if(k.eq.0)km=n3(klev)-1
          if(k.eq.n3(klev))ka=1
          do 14 j=1,n2(klev)
            if(j.eq.1)then
              pdmaa(j)=0.d0
              pdmab(j)=0.d0
              pdmac(j)=-2.d0/hhz(klev)-2.d0/hhp(klev)/rval**2
     +                +1.d0/rval**2
              pdmad(j)=5.d0/3.d0/hhz(klev)
              pdmae(j)=-0.2d0/hhz(klev)
            elseif(j.eq.n2(klev))then
              pdmaa(j)=-0.2d0/hhz(klev)
              pdmab(j)=5.d0/3.d0/hhz(klev)
              pdmac(j)=-2.d0/hhz(klev)-2.d0/hhp(klev)/rval**2
     +                +1.d0/rval**2
              pdmad(j)=0.d0
              pdmae(j)=0.d0
            else
              pdmaa(j)=0.d0
              pdmab(j)=1.d0/hhz(klev)
              pdmac(j)=-(2.d0/hhz(klev)+2.d0/hhp(klev)/rval**2)
     +                +1.d0/rval**2
              pdmad(j)=1.d0/hhz(klev)
              pdmae(j)=0.d0
            endif
c
            if(i.eq.1)then
            pdmac(j)=pdmac(j)-2.d0/hhr(klev)     !mpk +1.d0/rval/hr(klev)
            pdmar(j)=-(5.d0/3.d0*psif(i+1,j,k)-0.2d0*psif(i+2,j,k))/hhr
     +                (klev)-4.d0/3.d0*psif(i+1,j,k)/hr(klev)/rval
     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            elseif(i.eq.n1(klev))then
            pdmac(j)=pdmac(j)-2.d0/hhr(klev)      !mpk -1.d0/rval/hr(klev)
            pdmar(j)=-(5.d0/3.d0*psif(i-1,j,k)-0.2d0*psif(i-2,j,k))/
     +                hhr(klev)+4.d0/3.d0*psif(i-1,j,k)/hr(klev)/rval
     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            else
              pdmac(j)=pdmac(j)-2.d0/hhr(klev)
              pdmar(j)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)-1.5d0*
     +                (psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)/rval
     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            endif
   14     continue
          call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n2(klev))
          do 15 j=1,n2(klev)
            psif(i,j,k)=pdmas(j)
            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
   15     continue
   13   continue
c                             line relax in phi-dir
        do 16 j=1,n2(klev)
        do 16 i=1,n1(klev)
            rval=r0+(dble(i)-0.5d0)*hr(klev)
          do 17 k=1,n3(klev)
            tdmaa(k)=1.d0/hhp(klev)/rval**2
            tdmac(k)=1.d0/hhp(klev)/rval**2
            tdmab(k)=-(2.d0/hhp(klev)/rval**2)+1.d0/rval**2
            if(i.eq.1)then
            tdmab(k)=tdmab(k)-2.d0/hhr(klev)      !mpk +1.d0/rval/hr(klev)
            tdmar(k)=-(5.d0/3.d0*psif(i+1,j,k)-0.2d0*psif(i+2,j,k))/hhr
     +                (klev)-4.d0/3.d0*psif(i+1,j,k)/hr(klev)/rval
            elseif(i.eq.n1(klev))then
            tdmab(k)=tdmab(k)-2.d0/hhr(klev)       !mpk -1.d0/rval/hr(klev)
            tdmar(k)=-(5.d0/3.d0*psif(i-1,j,k)-0.2d0*psif(i-2,j,k))
     +              /hhr(klev)+4.d0/3.d0*psif(i-1,j,k)/hr(klev)/rval
            else
             tdmab(k)=tdmab(k)-2.d0/hhr(klev)
             tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)-1.5d0*
     +                (psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)/rval
            endif
c
            if(j.eq.1)then
             tdmab(k)=tdmab(k)-2.d0/hhz(klev)
             tdmar(k)=tdmar(k)-(5.d0/3.d0*psif(i,j+1,k)
     +               -0.2d0*psif(i,j+2,k))/hhz(klev)+rf(i,j,k)
            elseif(j.eq.n2(klev))then
             tdmab(k)=tdmab(k)-2.d0/hhz(klev)
             tdmar(k)=tdmar(k)-(5.d0/3.d0*psif(i,j-1,k)
     +               -0.2d0*psif(i,j-2,k))/hhz(klev)+rf(i,j,k)
            else
             tdmab(k)=tdmab(k)-2.d0/hhz(klev)
             tdmar(k)=tdmar(k)-(psif(i,j-1,k)+psif(i,j+1,k))
     +               /hhz(klev)+rf(i,j,k)
            endif
c
   17     continue
          call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n3(klev))
          do 18 k=0,n3(klev)
            if(k.eq.0)then
              psif(i,j,k)=tdmas(n3(klev))
            else
              psif(i,j,k)=tdmas(k)
            endif
   18     continue
   16   continue
c
c                                  Calculating residual.
c
          resmax=0.d0
          do 20 i=1,n1(klev)
             rval=r0+(dble(i)-0.5d0)*hr(klev)
          do 20 j=1,n2(klev)
          do 20 k=1,n3(klev)
          ka=k+1
          km=k-1
          if(k.eq.0)km=n3(klev)-1
          if(k.eq.n3(klev))ka=1
          temp=(psif(i,j,ka)-2.d0*psif(i,j,k)+psif(i,j,km))
     +        /hhp(klev)/rval**2+psif(i,j,k)/rval**2
          if(i.eq.1)then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0/3.d0*psif(i+1,j,k)-
     +           0.2d0*psif(i+2,j,k))/hhr(klev)+
     +           4.d0/3.d0*psif(i+1,j,k)/hr(klev)/rval
          elseif(i.eq.n1(klev))then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0/3.d0*psif(i-1,j,k)
     +          -0.2d0*psif(i-2,j,k))/hhr(klev)-
     +          4.d0/3.d0*psif(i-1,j,k)/hr(klev)/rval
          else
            temp=temp+(psif(i-1,j,k)-2.d0*psif(i,j,k)+psif(i+1,j,k))
     +          /hhr(klev)+1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +          /rval
          endif
c
          if(j.eq.1)then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0/3.d0*psif(i,j+1,k)
     +          -0.2d0*psif(i,j+2,k))/hhz(klev)
          elseif(j.eq.n2(klev))then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0/3.d0*psif(i,j-1,k)
     +          -0.2d0*psif(i,j-2,k))/hhz(klev)
          else
            temp=temp+(psif(i,j+1,k)-2.d0*psif(i,j,k)+psif(i,j-1,k))
     +          /hhz(klev)
          endif
c
            temp=rf(i,j,k)-temp
            rc(i,j,k)=temp
            temp=dabs(temp)
            if(temp.gt.resmax)resmax=temp
   20     continue
c
c         if(klev.eq.1)write(6,*)klev,resmax
c         if(klev.eq.1)call rufcon(vel,n1(klev),n2(klev),n3(klev))
c         write(6,100)klev,resmax
c       end if
  100   format(1x,'Level',i1,'  resmax=',f14.10)
c
  200 continue
c
      call trans1(psif,v,ibeg(klev),n1(klev),n2(klev),n3(klev))
      if(klev.eq.1.and.resmax.lt.tolcon)nconv=1

      return
      end
c 
c----------------------------------------------------------------------c
      subroutine relax3(nrel,klev,nconv,nlev)
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      real*8 pdmaa(kp),pdmab(kp),pdmac(kp),pdmad(kp),pdmae(kp),pdmas(kp)
      real*8 pdmar(kp)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
      common/iteration/kt
c
      tolcon=1.d-5
      nconv=0
c
      call trans2(f,rf,ibeg(klev),n1(klev),n2(klev),n3(klev))
      call trans2(v,psif,ibeg(klev),n1(klev),n2(klev),n3(klev))
c
c------------------------------relaxing with some G & S.
c
      do 200 irel=1,nrel
c                                 line relax in r-dir        
        do 10 k=1,n3(klev)
        do 10 j=1,n2(klev)
          ka=k+1
          km=k-1
          if(k.eq.n3(klev))ka=1
          do 11 i=1,n1(klev)-1
            rval=dble(i)*hr(klev)+r0
            tdmaa(i)=1.d0/hhr(klev)-1.5d0/hr(klev)/rval
            tdmab(i)=-(2.d0/hhr(klev)-1.d0/rval**2
     +              +2.d0/hhp(klev)/rval**2)
            tdmac(i)=1.d0/hhr(klev)+1.5d0/hr(klev)/rval
            if(j.eq.1)then
              tdmab(i)=tdmab(i)-2.d0/hhz(klev)
              tdmar(i)=-(5.d0*psif(i,j+1,k)/3.d0-0.2*psif(i,j+2,k))
     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            elseif(j.eq.n2(klev))then
              tdmab(i)=tdmab(i)-2.d0/hhz(klev)
              tdmar(i)=-(5.d0*psif(i,j-1,k)/3.d0-0.2*psif(i,j-2,k))
     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            else
              tdmab(i)=tdmab(i)-2.d0/hhz(klev)
              tdmar(i)=-(psif(i,j-1,k)+psif(i,j+1,k))/hhz(klev)
     +               -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            endif
   11     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n1(klev)-1)
          do 12 i=1,n1(klev)-1
            psif(i,j,k)=tdmas(i)
            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
   12     continue
   10   continue
c                                 line relax in z-dir
        do 13 k=1,n3(klev)
        do 13 i=1,n1(klev)-1
            rval=dble(i)*hr(klev)+r0
          ka=k+1
          km=k-1
          if(k.eq.0)km=n3(klev)-1
          if(k.eq.n3(klev))ka=1
          do 14 j=1,n2(klev)
            if(j.eq.1)then
              pdmaa(j)=0.d0
              pdmab(j)=0.d0
              pdmac(j)=-2.d0/hhz(klev)-2.d0/hhr(klev)+1.d0/rval**2
     +                -2.d0/hhp(klev)/rval**2
              pdmad(j)=5.d0/hhz(klev)/3.d0
              pdmae(j)=-0.2d0/hhz(klev)
            elseif(j.eq.n2(klev))then
              pdmaa(j)=-0.2d0/hhz(klev)
              pdmab(j)=5.d0/hhz(klev)/3.d0
              pdmac(j)=-2.d0/hhz(klev)-2.d0/hhr(klev)+1.d0/rval**2
     +                -2.d0/hhp(klev)/rval**2
              pdmad(j)=0.d0
              pdmae(j)=0.d0
            else
              pdmaa(j)=0.d0
              pdmab(j)=1.d0/hhz(klev)
              pdmac(j)=-(2.d0/hhz(klev)+2.d0/hhr(klev)+2.d0/hhp(klev)
     +              /rval**2)+1.d0/rval**2
              pdmad(j)=1.d0/hhz(klev)
              pdmae(j)=0.d0
            endif
            pdmar(j)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +               /rval-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
     +               /rval**2+rf(i,j,k)
   14     continue
          call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n2(klev))
          do 15 j=1,n2(klev)
            psif(i,j,k)=pdmas(j)
            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
   15     continue
   13   continue
c                                line relax in phi-dir
        do 16 j=1,n2(klev)
        do 16 i=1,n1(klev)-1
            rval=dble(i)*hr(klev)+r0
          do 17 k=1,n3(klev)
            tdmaa(k)=1.d0/hhp(klev)/rval**2
            tdmab(k)=-(2.d0/hhr(klev)-1.d0/rval**2
     +              +2.d0/hhp(klev)/rval**2)
            tdmac(k)=1.d0/hhp(klev)/rval**2
            if(j.eq.1)then
              tdmab(k)=tdmab(k)-2.d0/hhz(klev)
              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +               /rval-(5.d0*psif(i,j+1,k)/3.d0
     +               -0.2*psif(i,j+2,k))/hhz(klev)+rf(i,j,k)
            elseif(j.eq.n2(klev))then
              tdmab(k)=tdmab(k)-2.d0/hhz(klev)
              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +               /rval-(5.d0*psif(i,j-1,k)/3.d0
     +               -0.2*psif(i,j-2,k))/hhz(klev)+rf(i,j,k)
            else
              tdmab(k)=tdmab(k)-2.d0/hhz(klev)
              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +               /rval-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)
     +               +rf(i,j,k)
            endif
   17     continue
          call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n3(klev))
          do 18 k=0,n3(klev)
            if(k.eq.0)then
              psif(i,j,k)=tdmas(n3(klev))
            else
              psif(i,j,k)=tdmas(k)
            endif
   18     continue
   16   continue
c
c                                 Calculating residual.
c
          resmax=0.d0
          do 20 i=1,n1(klev)-1
            rval=dble(i)*hr(klev)+r0
          do 20 j=1,n2(klev)
          do 20 k=1,n3(klev)
          ka=k+1
          km=k-1
          if(k.eq.n3(klev))ka=1
          temp=(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +        +1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +        /rval+psif(i,j,k)/rval**2
     +        +(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
     +        -psif(i,j,k)*(2.d0/hhr(klev)+2.d0/hhp(klev)/rval**2)
          if(j.eq.1)then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0*psif(i,j+1,k)/3.d0
     +          -0.2d0*psif(i,j+2,k))/hhz(klev)
          elseif(j.eq.n2(klev))then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0*psif(i,j-1,k)/3.d0
     +          -0.2d0*psif(i,j-2,k))/hhz(klev)
          else
            temp=temp+(psif(i,j+1,k)-2.d0*psif(i,j,k)+psif(i,j-1,k))
     +          /hhz(klev)
          endif
            temp=rf(i,j,k)-temp
            rc(i,j,k)=temp
            temp=dabs(temp)
            if(temp.gt.resmax)resmax=temp
   20     continue
c
c         if(klev.eq.1)write(6,100)klev,resmax
c         write(6,100)klev,resmax
c       end if
  100   format(1x,'Level',i1,'  resmax=',f14.10)
c
  200 continue
c
      call trans1(psif,v,ibeg(klev),n1(klev),n2(klev),n3(klev))
      if(klev.eq.1.and.resmax.lt.tolcon)nconv=1

      return
      end
c
c----------------------------------------------------------------------c
      subroutine corser(klev)
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
c
c------------------Start off with restriction
c
      call trans3(rc,rf,n1(klev),n2(klev),n3(klev))
c
      do 10 i=1,n1(klev+1)-1
      do 10 j=1,n2(klev+1)-1
      do 10 k=1,n3(klev+1)-1
        ii=i+i
        jj=j+j
        kk=k+k
        iip=ii+1
        iim=ii-1
        jjp=jj+1
        jjm=jj-1
        kkp=kk+1
        kkm=kk-1
        if(k.eq.0)jjm=n3(klev)-1
        if(k.eq.n3(klev+1))kkp=1
c
        psic(i,j,k)=0.d0
c
        temp=rf(ii,jj,kk)
        temp=temp+temp+rf(iip,jj,kk)+rf(iim,jj,kk)
     +                +rf(ii,jjp,kk)+rf(ii,jjm,kk)
     +                +rf(ii,jj,kkp)+rf(ii,jj,kkm)
        temp=temp+temp+rf(iip,jjp,kk)+rf(iip,jjm,kk)
     +                +rf(iim,jjp,kk)+rf(iim,jjm,kk)
     +                +rf(iip,jj,kkp)+rf(iip,jj,kkm)
     +                +rf(iim,jj,kkp)+rf(iim,jj,kkm)
     +                +rf(ii,jjp,kkp)+rf(ii,jjp,kkm)
     +                +rf(ii,jjm,kkp)+rf(ii,jjm,kkm)
        temp=temp+temp+rf(iip,jjp,kkp)+rf(iip,jjp,kkm)
     +                +rf(iip,jjm,kkp)+rf(iip,jjm,kkm)
     +                +rf(iim,jjm,kkp)+rf(iim,jjm,kkm)
     +                +rf(iim,jjp,kkp)+rf(iim,jjp,kkm)
c        temp=temp/64.d0
        temp=temp*0.015625d0
        rc(i,j,k)=temp
   10 continue
c
c-----------setting faces to zero
c
      do 20 j=0,n2(klev+1)
      do 20 k=0,n3(klev+1)
        rc(0,j,k)=0.d0
        rc(n1(klev+1),j,k)=0.d0
        psic(0,j,k)=0.d0
        psic(n1(klev+1),j,k)=0.d0
   20 continue
c
      do 30 i=0,n1(klev+1)
      do 30 k=0,n3(klev+1)
        rc(i,0,k)=0.d0
        rc(i,n2(klev+1),k)=0.d0
        psic(i,0,k)=0.d0
        psic(i,n2(klev+1),k)=0.d0
   30 continue
c
c     do 40 i=0,n1(klev+1)
c     do 40 j=0,n2(klev+1)
c       rc(i,j,0)=0.d0
c       rc(i,j,n3(klev+1))=0.d0
c       psic(i,j,0)=0.d0
c       psic(i,j,n3(klev+1))=0.d0
c  40 continue
c
      call trans1(rc,f,ibeg(klev+1),n1(klev+1),n2(klev+1),n3(klev+1))
      call trans1(psic,v,ibeg(klev+1),n1(klev+1),n2(klev+1),n3(klev+1))
c
c     print*,'In corser:    PSI then R'
c     call rufcon(psic,n1(klev+1),n2(klev+1),n3(klev+1))
c     call rufcon(rc,n1(klev+1),n2(klev+1),n3(klev+1))
c
      return
      end
c
c----------------------------------------------------------------------c
      subroutine finer(klev)
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz,0:kp),psif(0:kr,0:kz,0:kp)
      real*8 rc(0:kr,0:kz,0:kp),rf(0:kr,0:kz,0:kp),tempa(0:kr,0:kz,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
c
      call trans2(v,psic,ibeg(klev),n1(klev),n2(klev),n3(klev))
c
c--------------------prolongation
c
      do 10 i=0,n1(klev-1)
      do 10 j=0,n2(klev-1)
      do 10 k=0,n3(klev-1)
        i1=i/2
        j1=j/2
        k1=k/2
        i2=i1*2
        j2=j1*2
        k2=k1*2
        if(i2.eq.i)then
          if(j2.eq.j)then
            if(k2.eq.k)then
              psif(i,j,k)=psic(i1,j1,k1)
            else
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1,j1,k1+1))*0.5d0
            end if
          else
            if(k2.eq.k)then
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1,j1+1,k1))*0.5d0
            else
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1,j1+1,k1)
     +                    +psic(i1,j1,k1+1)+psic(i1,j1+1,k1+1))*0.25d0
            end if
          end if
        else
          if(j2.eq.j)then
            if(k2.eq.k)then
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1+1,j1,k1))*0.5d0
            else
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1+1,j1,k1)
     +                    +psic(i1,j1,k1+1)+psic(i1+1,j1,k1+1))*0.25d0
            end if
          else
            if(k2.eq.k)then
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1+1,j1,k1)
     +                    +psic(i1,j1+1,k1)+psic(i1+1,j1+1,k1))*0.25d0
            else
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1+1,j1,k1)
     +                     +psic(i1,j1+1,k1)+psic(i1+1,j1+1,k1)
     +                     +psic(i1,j1,k1+1)+psic(i1+1,j1,k1+1)
     +               +psic(i1,j1+1,k1+1)+psic(i1+1,j1+1,k1+1))*0.125d0
            end if
          end if
        end if
   10 continue
c
      call trans2(v,tempa,ibeg(klev-1),n1(klev-1),n2(klev-1),n3(klev-1))
c
      do 20 i=0,n1(klev-1)
      do 20 j=0,n2(klev-1)
      do 20 k=0,n3(klev-1)
        psif(i,j,k)=psif(i,j,k)+tempa(i,j,k)
        psic(i,j,k)=0.d0
   20 continue
      call trans1(psic,v,ibeg(klev),n1(klev),n2(klev),n3(klev))
      call trans1(psif,v,ibeg(klev-1),n1(klev-1),n2(klev-1),n3(klev-1))
      return
      end
c==================================================
      subroutine trans1(a3d,a1d,i1,n1,n2,n3)
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 a3d(0:kr,0:kz,0:kp),a1d(mgen)
c
      iimax=0
      do 10 k=0,n3
      do 10 j=0,n2
      do 10 i=0,n1
        ii=k*(n1+1)*(n2+1) + j*(n1+1) + i + i1
!mpk        if(ii.gt.202590)then
!mpk         print*,'FATAL ERROR IN MULTIGRID  ii='
!mpk         STOP
!mpk        endif        
        a1d(ii)=a3d(i,j,k)
   10 continue
!mpk        if(ii.gt.317000)then
!mpk          print*,'FATAL ERROR IN MULTIGRID  ii='
!mpk          STOP
!mpk         endif 
!mpk        print*,ii
      return
      end
c==================================================
      subroutine trans2(a1d,a3d,i1,n1,n2,n3)
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 a3d(0:kr,0:kz,0:kp),a1d(mgen)
c
      do 10 k=0,n3
      do 10 j=0,n2
      do 10 i=0,n1
        ii=k*(n1+1)*(n2+1) + j*(n1+1) + i + i1
        a3d(i,j,k)=a1d(ii)
   10 continue
      return
      end
c==================================================
      subroutine trans3(a3d1,a3d2,n1,n2,n3)
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=160,mgen=202590)
      real*8 a3d1(0:kr,0:kz,0:kp),a3d2(0:kr,0:kz,0:kp)
c
      do 10 k=0,n3
      do 10 j=0,n2
      do 10 i=0,n1
        a3d2(i,j,k)=a3d1(i,j,k)
   10 continue
      return
      end
c==================================================
      subroutine rufcon(a,nx,ny,nz)
      parameter(kx=33,ky=33,kz=160,mgen=202590)
      implicit real*8(a-h,o-z)
      real*8 a(0:kx,0:ky,0:kz)
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
c
c--------Z-section at Z=Zmax/2
c
      amax=biggest2(a,nx,ny,nz,ny/2)
      amin=smallest2(a,nx,ny,nz,ny/2)
      print*,' '
      print*,'At z=zmax/2'
      print*,'arraymax,min=',amax,amin
      print*,'============'
      print*,'Abscissa = r-axis, Ordinate = phi-axis'
c
      adiff=amax-amin
      if(adiff.ne.0.)adiff=1./adiff
c
      do 20 k=nz,0,-1
        do 10 i=0,nx
          x=1.5+10.*(a(i,ny/2,k)-amin)*adiff
          n=x
          c(i)=ch(n)
  10    continue
        write(6,100)(c(i),i=0,nx)
  20  continue
c
c--------Z-section at r=rmax/2
c
      amax=biggest1(a,nx,ny,nz,nx/2)
      amin=smallest1(a,nx,ny,nz,nx/2)
      print*,' '
      print*,'At r=rmax/2'
      print*,'arraymax,min=',amax,amin
      print*,'============'
      print*,'Abscissa = z-axis, Ordinate = phi-axis'
c
      adiff=amax-amin
      if(adiff.ne.0.)adiff=1./adiff
c
      do 30 k=nz,0,-1
        do 40 j=0,ny
          x=1.5+10.*(a(nx/2,j,k)-amin)*adiff
          n=x
          c(j)=ch(n)
  40    continue
        write(6,100)(c(j),j=0,ny)
  30  continue
c
c
c--------Z-section at phi=phimax/2
c
      amax=biggest3(a,nx,ny,nz,nz/2)
      amin=smallest3(a,nx,ny,nz,nz/2)
      print*,' '
      print*,'At phi=phimax/2'
      print*,'arraymax,min=',amax,amin
      print*,'============'
      print*,'Abscissa = z-axis, Ordinate = r-axis'
c
      adiff=amax-amin
      if(adiff.ne.0.)adiff=1./adiff
c
      do 50 i=nx,0,-1
        do 60 j=0,ny
          x=1.5+10.*(a(i,j,nz/2)-amin)*adiff
          n=x
          c(j)=ch(n)
  60    continue
        write(6,100)(c(j),j=0,ny)
  50  continue
c
 100  format(1x,200a1)
      return
      end
c==================================================
      double precision function smallest1(a,nx,ny,nz,ii)
      parameter(kx=33,ky=33,kz=160,mgen=202590)
      implicit real*8(a-h,o-z)
      real*8 a(0:kx,0:ky,0:kz)
c
      z=a(ii,0,0)
      do 10 j=0,ny
      do 10 k=0,nz
        if(z.gt.a(ii,j,k))z=a(ii,j,k)
  10  continue
      smallest1=z
      return
      end
c==================================================
      double precision function biggest1(a,nx,ny,nz,ii)
      parameter(kx=33,ky=33,kz=160,mgen=202590)
      implicit real*8(a-h,o-z)
      real*8 a(0:kx,0:ky,0:kz)
c
      z=a(ii,0,0)
      do 10 j=0,ny
      do 10 k=0,nz
        if(z.lt.a(ii,j,k))z=a(ii,j,k)
  10  continue
      biggest1=z
      return
      end
c==================================================
      double precision function smallest2(a,nx,ny,nz,jj)
      parameter(kx=33,ky=33,kz=160,mgen=202590)
      implicit real*8(a-h,o-z)
      real*8 a(0:kx,0:ky,0:kz)
c
      z=a(0,jj,0)
      do 10 i=0,nx
      do 10 k=0,nz
        if(z.gt.a(i,jj,k))z=a(i,jj,k)
  10  continue
      smallest2=z
      return
      end
c==================================================
      double precision function biggest2(a,nx,ny,nz,jj)
      parameter(kx=33,ky=33,kz=160,mgen=202590)
      implicit real*8(a-h,o-z)
      real*8 a(0:kx,0:ky,0:kz)
c
      z=a(0,jj,0)
      do 10 i=0,nx
      do 10 k=0,nz
        if(z.lt.a(i,jj,k))z=a(i,jj,k)
  10  continue
      biggest2=z
      return
      end
c==================================================
      double precision function smallest3(a,nx,ny,nz,kk)
      parameter(kx=33,ky=33,kz=160,mgen=202590)
      implicit real*8(a-h,o-z)
      real*8 a(0:kx,0:ky,0:kz)
c
      z=a(0,0,kk)
      do 10 i=0,nx
      do 10 j=0,ny
        if(z.gt.a(i,j,kk))z=a(i,j,kk)
  10  continue
      smallest3=z
      return
      end
c==================================================
      double precision function biggest3(a,nx,ny,nz,kk)
      parameter(kx=33,ky=33,kz=160,mgen=202590)
      implicit real*8(a-h,o-z)
      real*8 a(0:kx,0:ky,0:kz)
c
      z=a(0,0,kk)
      do 10 i=0,nx
      do 10 j=0,ny
        if(z.lt.a(i,j,kk))z=a(i,j,kk)
  10  continue
      biggest3=z
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
      subroutine ptdma(a,b,c,y,r,n)
      implicit real*8(a-h,o-z)
      real*8 a(n),b(n),c(n),y(n),r(n)
c
c     Does a periodic tdma algorithm.
c
c     if(i1.eq.5.and.i2.eq.5)then
c       do 5 i=0,n
c         print*,a(i),b(i),c(i),y(i),r(i)
c   5   continue
c     end if
c
      do 10 j=2,n-2
        fact=a(j)/b(j-1)
        b(j)=b(j)-fact*c(j-1)
        a(j)=-fact*a(j-1)
        r(j)=r(j)-fact*r(j-1)
        fact=c(n)/b(j-1)
        c(n)=-fact*c(j-1)
        b(n)=b(n)-fact*a(j-1)
        r(n)=r(n)-fact*r(j-1)
   10 continue
      fact=a(n-1)/b(n-2)
      b(n-1)=b(n-1)-fact*c(n-2)
      c(n-1)=c(n-1)-fact*a(n-2)
      r(n-1)=r(n-1)-fact*r(n-2)
      fact=c(n)/b(n-2)
      a(n)=a(n)-fact*c(n-2)
      b(n)=b(n)-fact*a(n-2)
      r(n)=r(n)-fact*r(n-2)
      fact=a(n)/b(n-1)
      b(n)=b(n)-fact*c(n-1)
      r(n)=r(n)-fact*r(n-1)
c
c     det=1.d0
c     do 15 i=1,n
c       det=det*b(i)
c  15 continue
c     print*,det
c
c
      y(n)=r(n)/b(n)
      y(n-1)=(r(n-1)-c(n-1)*y(n))/b(n-1)
c
c
      do 20 j=n-2,1,-1
        y(j)=(r(j)-c(j)*y(j+1)-a(j)*y(n))/b(j)
   20 continue
c
      return
      end
c==================================================
      subroutine pdma(a,b,c,d,e,r,x,n)
      implicit real*8(a-h,o-z)
      real*8 a(n),b(n),c(n),d(n),e(n)
      real*8 r(n),x(n)
c
c     Solves a  pentadiagonal matrix system:
c     Solves a  pentadiagonal matrix system:
c
c            /                \   /  \   /  \
c            | C1 D1 E1       |   |x1|   |r1|
c            | B2 C2 D2 E2    |   |x2|   |r2|
c            | A3 B3 C3 D3 E3 | * |x3| = |r3|
c            |    A4 B4 C4 D4 |   |x4|   |r4|
c            |       A5 B5 C5 |   |x5|   |r5|
c            \                /   \  /   \  /
c
      do 10 k=1,n-1
        d(k)=d(k)/c(k)
        e(k)=e(k)/c(k)
        r(k)=r(k)/c(k)
        c(k)=1.d0
c
        c(k+1)=c(k+1)-b(k+1)*d(k)
        if(k.ne.(n-1))d(k+1)=d(k+1)-b(k+1)*e(k)
        r(k+1)=r(k+1)-b(k+1)*r(k)
        b(k+1)=0.d0
c
        if(k.ne.(n-1))then
          b(k+2)=b(k+2)-a(k+2)*d(k)
          c(k+2)=c(k+2)-a(k+2)*e(k)
          r(k+2)=r(k+2)-a(k+2)*r(k)
        end if
c
   10 continue
c
c--------------Back substitution bit
c
      x(n)=r(n)/c(n)
c
      x(n-1)=r(n-1)-d(n-1)*x(n)
c
      do 20 j=n-2,1,-1
        x(j)=r(j)-d(j)*x(j+1)-e(j)*x(j+2)
   20 continue
c
      return
      end
c================================================
