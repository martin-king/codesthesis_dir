c23456789012345678901234567890123456789012345678901234567890123456789012
c     Program sect360v2.f                                              c
c     ====================                                             c
c!mpk located disk14/enptwl/proj14.dir/rotstag14.f                     c
c     coded originally by T.Lewis(1999)                                c
c     Martin P. King, University of Bath                               c
c     Mathematical modelling of flow in a rotating 360deg sector       c
c     with heat flux applied in the radial direction                   c
c     Three dimensional velocity/vorticity model on staggered grid     c
c         .retains original non-d scalings by Lewis                    c
c         .introduces gravitational terms(present only in Avort and    c
c          Bvort. A Grashof number(Gr) entails.                        c
c         .radial and axial grid sizes are increased.                  c
c         .hybrid scheme introduced in ABCVort and Temp(can be commented
c          out). current status: ON                                    c
c         .corser2 & corser3 added                                     c
c         .change u & w located at face-centre v at walls              c
c         .see line 115:line pert (ON)                                 c
c         .starts with lukewarm fluid(line 112)(OFF)                   c
c         .line 795 switched OFF(silly assumptions)                    c                           c
c         .output local Nu, see line 309                               c
c         .extrapolation and SOR added in subroutine POISSON?          c
c         LAST UPDATED: 12-9-00.                                       c
c----------------------------------------------------------------------c
c!mpk      real*8 ssumc(6,mgen),ssums(6,mgen)
c!mpk      common/grals/ssumc,ssums
c!mpk      common/iteration/kt
c
      INCLUDE "sect360mg.ina"
c
      real*8 httrans(1:mgen)
      character*20 fname1,fname2,fname3,fname4,fname5
      character*20 fname6,fname7,fname8,fname9,ans
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
      print*,'number of intervals in phi-direction for a 45 deg sector:'
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
c!mpk      print*,'axial temperature perturbation factor'
c!mpk      read*,pert
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
      do 10 k=0,np+1
        An(i,j,k)=0.d0
        A(i,j,k)=0.d0
        Bn(i,j,k)=0.d0
        B(i,j,k)=0.d0
        Cn(i,j,k)=0.d0
        C(i,j,k)=0.d0
        uel(i,j,k)=0.d0
        vel(i,j,k)=0.d0
        wel(i,j,k)=0.d0
        if(i.eq.0)then
          r=r0
        elseif(i.eq.nr+1)then
          r=rmax
        else
          r=r0+(dble(i)-0.5d0)*dr
        endif
!mpk        th(i,j,k)=0.5d0
!mpk        th(nr+1,j,k)=1.d0
c!mpk        z=dble(j)*dz
        phi=8.d0*datan(1.d0)*dble(k)/dble(np)
        th(i,j,k)=1.d0-dlog(r)/dlog(r0)
c!mpk    +      +pert*dsin(2.d0*pi*(r/rmax-0.5d0))*dcos(2.d0*pi*z)
     +      +amp*(rmax-r)*(r-r0)*dcos(dble(nval)*phi)
c!mpk    +      +amp*dcos(pi*z/zmax)*(rmax-r)*(r-r0)*dcos(dble(nval)*phi)
c        th(nr/4,j,np/2)=th(nr/4,j,np/2)-0.2d0
  10  continue
      do 11 j=1,nz
           th(nr/4,j,np/2)=th(nr/4,j,np/2)-0.2d0
  11  continue
c
c!mpk      do 12 i=1,nr
c!mpk      do 12 j=1,nz
c!mpk      do 12 k=0,np
c!mpk        r=r0+dble(i)*dr
c!mpk        A(i,j,k)=-(vel(i,j+1,k)-vel(i,j,k))/dz
c!mpk        A(i,0,k)=-dum10(vel,dz,i,0,k)
c!mpk        A(i,nz,k)=-dum12(vel,dz,i,nz+1,k)
c!mpk        C(i,j,k)=(vel(i+1,j,k)-vel(i,j,k))/dr
c!mpk     +           +0.5d0*(vel(i+1,j,k)+vel(i,j,k))/r
c!mpk        C(0,j,k)=dum9(vel,dr,0,j,k)
c!mpk        C(nr,j,k)=dum11(vel,dr,nr+1,j,k)+vel(nr+1,j,k)/rmax
c!mpk        A(i,0,k)=-2.d0*vel(i,1,k)/dz
c!mpk        A(i,nz,k)=2.d0*vel(i,nz,k)/dz
c!mpk        C(0,j,k)=2.d0*vel(1,j,k)/dr
c!mpk        C(nr,j,k)=-2.d0*vel(nr,j,k)/dr
c!mpk  12  continue
c
!mpk      call rufcon(th,nr+1,nz+1,np+1)
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
        read(11)(((B(i,j,k),i=0,nr),j=0,nz),k=0,np+1)
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
        read(11)(((uel(i,j,k),i=0,nr),j=0,nz+1),k=0,np+1)
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
        read(11)(((wel(i,j,k),i=0,nr+1),j=0,nz),k=0,np+1)
        close(11)
c
        print*,'filename for temperature: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)nr,nz,np
        read(11)rmax,r0,zmax,zmin
        read(11)(((th(i,j,k),i=0,nr+1),j=0,nz+1),k=0,np+1)
        close(11)
      endif
        call rufcon(th,nr+1,nz+1,np+1)
c
      do 47 i=0,nr+1
      do 47 j=0,nz
      do 47 k=0,np
        Ap(i,j,k)=A(i,j,k)
  47  continue
      do 48 i=0,nr
      do 48 j=0,nz
      do 48 k=0,np+1
        Bp(i,j,k)=B(i,j,k)
  48  continue
      do 49 i=0,nr
      do 49 j=0,nz+1
      do 49 k=0,np
        Cp(i,j,k)=C(i,j,k)
  49  continue
      do 50 i=0,nr+1
      do 50 j=0,nz+1
      do 50 k=0,np+1
        thp(i,j,k)=th(i,j,k)
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
c!mpk        call aaaaaaaa(ntot+kt)
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
c        call rufcon(th,nr+1,nz+1,np+1)
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
  30  continue
c
c----------------------------------------------------------------------c
c     Write vorticity, stream function, velocity and temperature data  c
c     to file                                                          c
c----------------------------------------------------------------------c
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
      write(11)(((B(i,j,k),i=0,nr),j=0,nz),k=0,np+1)
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
      write(11)(((uel(i,j,k),i=0,nr),j=0,nz+1),k=0,np+1)
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
      write(11)(((wel(i,j,k),i=0,nr+1),j=0,nz),k=0,np+1)
      close(11)
c
      open(11,file=fname7,form='unformatted')
      write(11)nr,nz,np
      write(11)rmax,r0,zmax,zmin
      write(11)(((th(i,j,k),i=0,nr+1),j=0,nz+1),k=0,np+1)
      close(11)
c
      if(fname8 .ne. 'no')then
        open(11,file=fname8,status='unknown')
        write(11,200)((i-1)*dt,httrans(i),i=1,ntot)
cmpk--cal of local Nu at innner & outer cyls at mid-axial--------------c
      do 31 k=0,np 
        httrans_locala=-r0*dlog(r0)*(-184.d0*th(0,nz/2,k)+225.d0*
     +   th(1,nz/2,k)-50.d0*th(2,nz/2,k)+9.d0*th(3,nz/2,k))/(60.d0*dr)
        httrans_localb=-rmax*dlog(r0)*(184.d0*th(nr+1,nz/2,k)-225.d0*
     +   th(nr,nz/2,k)+50.d0*th(nr-1,nz/2,k)-9.d0*th(nr-2,nz/2,k))/
     +   (60.d0*dr)
         write(11,*)k,httrans_locala,httrans_localb
  31  continue
c!mpk        ,ssumc(1,i),ssums(1,i),
c!mpk     +  ssumc(2,i),ssums(2,i),ssumc(3,i),ssums(3,i),ssumc(4,i),
c!mpk     +  ssums(4,i),ssumc(5,i),ssums(5,i),ssumc(6,i),ssums(6,i)
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
      INCLUDE "sect360mg.ina"
c minus real*8 urh(0:kr,0:kz,0:kp),vrh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
c      
      do 10 i=1,nr
      do 10 j=1,nz-1
      do 10 k=1,np-1
c
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        km=k-1
        ka=k+1
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
         Arra=dum5a(A,dr,i,j,k)
         Ara=dum13a(A,dr,i,j,k)
         Ar=dum13(A,dr,i,j,k)
         cnr=-2.5d0/dr**2+0.5d0/dr/r
        elseif(i.eq.nr)then
         Arra=dum7a(A,dr,i,j,k)
         Ara=dum15a(A,dr,i,j,k)
         Ar=dum15(A,dr,i,j,k)
         cnr=-2.5d0/dr**2-0.5d0/dr/r
        else
         Arra=dum1a(A,dr,i,j,k)
         Ara=dum3(A,dr,i,j,k)
         Ar=Ara
         cnr=-1.d0/dr**2
c-------------------hybrid scheme introduced by mpk--------------------c
         Pe=dabs(u*dr/Ptl)
         if(Pe.gt.2)then
c          print*,'Pe=',Pe
          if(u.gt.0)then
             Ar=(A(i,j,k)-A(i-1,j,k))/dr
          elseif(u.lt.0)then
             Ar=(A(i+1,j,k)-A(i,j,k))/dr
          endif
         endif
c-----------------hybrid scheme ends for Ar in advectice term----------c
        endif
c
        uelr=((uel(i,j,k)+uel(i,jp,k)+uel(i,j,ka)+uel(i,jp,ka))*0.25d0
     +      -(uel(im,j,k)+uel(im,jp,k)+uel(im,j,ka)+uel(im,jp,ka))
     +      *0.25d0)/dr
        Azza=dum2a(A,dz,i,j,k)
        Az=dum4(A,dz,i,j,k)
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(w*dz/Ptl)
        if(Pe.gt.2)then
c            print*,'Pe=',Pe
         if(w.gt.0)then
           Az=(A(i,j,k)-A(i,j-1,k))/dz
         elseif(w.lt.0)then
           Az=(A(i,j+1,k)-A(i,j,k))/dz
         endif
        endif
c-----------------hybrid scheme ends for Az in advectice term----------c
c
        uelz=((uel(i,jp,k)+uel(i,jp,ka)+uel(im,jp,k)+uel(im,jp,ka))*
     +        0.25d0-(uel(i,j,k)+uel(i,j,ka)+uel(im,j,k)+uel(im,j,ka))
     +       *0.25d0)/dz
        thz=(th(i,jp,k)+th(i,jp,ka)-th(i,j,k)-th(i,j,ka))*0.5d0/dz
        thh=(th(i,j,ka)+th(i,jp,ka)-th(i,j,k)-th(i,jp,k))*0.5d0/dp
c
        Ahha=(A(i,j,ka)+A(i,j,km))/dp**2
        Ah=(A(i,j,ka)-A(i,j,km))*0.5d0/dp
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(v*dp*r/Ptl)
        if(Pe.gt.2)then
c         print*,'Pe=',Pe
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
     +       0.25d0-(uel(i,j,k)+uel(i,jp,k)+uel(im,j,k)+uel(im,jp,k))
     +      *0.25d0)/dp
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
      INCLUDE "sect360mg.ina"
c minus real*8 urh(0:kr,0:kz,0:kp),vrh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
c 
      do 10 i=1,nr-1
      do 10 j=1,nz-1
      do 10 k=1,np
c
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        km=k-1
        ka=k+1
c
        r=r0+dble(i)*dr
        u=(uel(i,j,k)+uel(i,jp,k))*0.5d0
        v=(vel(i,j,km)+vel(i,jp,km)+vel(i,j,k)+vel(i,jp,k)
     +   +vel(ip,j,km)+vel(ip,jp,km)+vel(ip,j,k)+vel(ip,jp,k))*0.125d0
        w=(wel(i,j,k)+wel(ip,j,k))*0.5d0
        tem=(th(i,jp,k)+th(ip,jp,k)+th(i,j,k)+th(ip,j,k))*.25d0
        AA=(A(i,j,k)+A(ip,j,k)+A(i,j,km)+A(ip,j,km))*0.25d0
        CC=(C(i,j,k)+C(i,jp,k)+C(i,j,km)+C(i,jp,km))*0.25d0
c
        Brra=dum1a(B,dr,i,j,k)
        Bra=dum3(B,dr,i,j,k)
        Br=Bra
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(u*dr/Ptl)
        if(Pe.gt.2)then
c         print*,'Pe=',Pe
         if(u.gt.0)then
           Br=(B(i,j,k)-B(i-1,j,k))/dr
         elseif(u.lt.0)then
           Br=(B(i+1,j,k)-B(i,j,k))/dr
         endif
        endif
c-----------------hybrid scheme ends for Br in advectice term----------c
c
        velr=((vel(ip,j,km)+vel(ip,jp,km)+vel(ip,j,k)+vel(ip,jp,k))*
     +       0.25d0-(vel(i,j,km)+vel(i,jp,km)+vel(i,j,k)+vel(i,jp,k))*
     +       0.25d0)/dr
        Bzza=dum2a(B,dz,i,j,k)
        Bz=dum4(B,dz,i,j,k)
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(w*dz/Ptl)
        if(Pe.gt.2)then
c         print*,'Pe=',Pe
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
c
        if(k.eq.1)then
         Bhha=(3.2d0*B(i,j,0)+2.d0*B(i,j,2)-0.2d0*B(i,j,3))/dp**2
         Bh=(-4.d0/3.d0*B(i,j,0)+B(i,j,1)+B(i,j,2)/3.d0)/dp
         cnh=-2.5d0/dp**2/r**2
        elseif(k.eq.np)then
         Bhha=(-0.2d0*B(i,j,np-2)+2.d0*B(i,j,np-1)+3.2d0*B(i,j,np+1))/
     +        dp**2
         Bh=(-B(i,j,np-1)/3.d0-B(i,j,np)+4.d0/3.d0*B(i,j,np+1))/dp
         cnh=-2.5d0/dp**2/r**2
        else
         Bhha=(B(i,j,ka)+B(i,j,km))/dp**2
         Bh=(B(i,j,ka)-B(i,j,km))*0.5d0/dp
         cnh=-1.d0/dp**2/r**2
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(v*dp*r/Ptl)
        if(Pe.gt.2)then
c         print*,'Pe=',Pe
         if(v.gt.0)then
           Bh=(B(i,j,k)-B(i,j,km))/dp
         elseif(v.lt.0)then
           Bh=(B(i,j,ka)-B(i,j,k))/dp
         endif
        endif
c-----------------hybrid scheme ends for Bh in advectice term----------c
        endif
c
        Ah=(A(i,j,k)+A(ip,j,k)-A(i,j,km)-A(ip,j,km))*0.5d0/dp
        velh=((vel(i,j,k)+vel(ip,j,k)+vel(i,jp,k)+vel(ip,jp,k))*0.25d0
     + -(vel(i,j,km)+vel(ip,j,km)+vel(i,jp,km)+vel(ip,jp,km))*0.25d0)/dp
c
        cnr=-1.d0/dr**2
        cnz=-1.d0/dz**2
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
      INCLUDE "sect360mg.ina"
c minus real*8 urh(0:kr,0:kz,0:kp),vrh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
c
      do 10 i=1,nr-1
      do 10 j=1,nz
      do 10 k=1,np-1
c
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        km=k-1
        ka=k+1
c
        r=r0+dble(i)*dr
        u=(uel(i,j,ka)+uel(i,j,k))*0.5d0
        v=(vel(i,j,k)+vel(ip,j,k))*0.5d0
        w=(wel(i,j,k)+wel(ip,j,k)+wel(i,j,ka)+wel(ip,j,ka)
     +   +wel(i,jm,k)+wel(ip,jm,k)+wel(i,jm,ka)+wel(ip,jm,ka))*0.125d0
        AA=(A(i,j,k)+A(i,jm,k)+A(ip,j,k)+A(ip,jm,k))*0.25d0
        BB=(B(i,j,k)+B(i,jm,k)+B(i,j,ka)+B(i,jm,ka))*0.25d0
        tem=(th(i,j,k)+th(ip,j,k)+th(ip,j,ka)+th(i,j,ka))*0.25d0
c
        Crra=dum1a(C,dr,i,j,k)
        Cra=dum3(C,dr,i,j,k)
        Cr=Cra
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(u*dr/Ptl)
        if(Pe.gt.2)then
c         print*,'Pe=',Pe
         if(u.gt.0)then
           Cr=(C(i,j,k)-C(i-1,j,k))/dr
         elseif(u.lt.0)then
           Cr=(C(i+1,j,k)-C(i,j,k))/dr
         endif
        endif
c-----------------hybrid scheme ends for Cr in advectice term----------c
c
        welr=((wel(ip,j,k)+wel(ip,j,ka)+wel(ip,jm,k)+wel(ip,jm,ka))*
     +   0.25d0-(wel(i,j,k)+wel(i,j,ka)+wel(i,jm,k)+wel(i,jm,ka))*
     +   0.25d0)/dr
        uelr=(uel(ip,j,k)+uel(ip,j,ka)-uel(im,j,k)-uel(im,j,ka))*
     +   0.25d0/dr
        thr=(th(ip,j,k)+th(ip,j,ka)-th(i,j,k)-th(i,j,ka))*0.5d0/dr
        if(j.eq.1)then
          Czza=dum6a(C,dz,i,j,k)
          Cz=dum14(C,dz,i,j,k)
          cnz=-2.5d0/dz**2
        elseif(j.eq.nz)then
          Czza=dum8a(C,dz,i,j,k)
          Cz=dum16(C,dz,i,j,k)
          cnz=-2.5d0/dz**2
        else
          Czza=dum2a(C,dz,i,j,k)
          Cz=dum4(C,dz,i,j,k)
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(w*dz/Ptl)
        if(Pe.gt.2)then
         if(w.gt.0)then
           Cz=(C(i,j,k)-C(i,j-1,k))/dz
         elseif(w.lt.0)then
           Cz=(C(i,j+1,k)-C(i,j,k))/dz
         endif
        endif
c-----------------hybrid scheme ends for Cz in advectice term----------c
          cnz=-1.d0/dz**2
        endif
        welz=((wel(i,j,k)+wel(ip,j,k)+wel(i,j,ka)+wel(ip,j,ka))*
     +     0.25d0-(wel(i,jm,k)+wel(ip,jm,k)+wel(i,jm,ka)+wel(ip,jm,ka))*
     +     0.25d0)/dz
        Chha=(C(i,j,ka)+C(i,j,km))/dp**2
        Ch=(C(i,j,ka)-C(i,j,km))*0.5d0/dp
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(v*dp*r/Ptl)
        if(Pe.gt.2)then
c          print*,'Pe=',Pe
         if(v.gt.0)then
           Ch=(C(i,j,k)-C(i,j,km))/dp
         elseif(v.lt.0)then
           Ch=(C(i,j,ka)-C(i,j,k))/dp
         endif
        endif
c-----------------hybrid scheme ends for Ch in advectice term----------c
c
        welh=((wel(i,j,ka)+wel(ip,j,ka)+wel(i,jm,ka)+wel(ip,jm,ka))*
     +   0.25d0-(wel(i,j,k)+wel(ip,j,k)+wel(i,jm,k)+wel(ip,jm,k))*
     +   0.25d0)/dp
        velh=(vel(i,j,ka)+vel(ip,j,ka)-vel(i,j,km)-vel(ip,j,km))*
     +   0.25d0/dp
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
      INCLUDE "sect360mg.ina"
c minus real*8 urh(0:kr,0:kz,0:kp),vrh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
c
      do 10 i=1,nr
      do 10 j=0,nz+1
      do 10 k=0,np+1
c
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        km=k-1
        ka=k+1
c
        r=r0+(dble(i)-0.5d0)*dr
        if((j.eq.0).or.(j.eq.nz+1).or.(k.eq.0).or.(k.eq.np+1))then
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
          thrra=dum5a(th,dr,i,j,k)
          thra=dum13a(th,dr,i,j,k)
          thr=dum13(th,dr,i,j,k)
          cnr=-2.5d0/dr**2+0.5d0/dr/r
        elseif(i.eq.nr)then
          thrra=dum7a(th,dr,i,j,k)
          thra=dum15a(th,dr,i,j,k)
          thr=dum15(th,dr,i,j,k)
          cnr=-2.5d0/dr**2-0.5d0/dr/r
        else
          thrra=dum1a(th,dr,i,j,k)
          thra=dum3(th,dr,i,j,k)
          thr=thra
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(u*dr/Ptl)
        if(Pe.gt.2)then
c         print*,'Pe=',Pe
         if(u.gt.0)then
           thr=(th(i,j,k)-th(i-1,j,k))/dr
         elseif(u.lt.0)then
           thr=(th(i+1,j,k)-th(i,j,k))/dr
         endif
        endif
c-----------------hybrid scheme ends for thr in advectice term---------c
          cnr=-1.d0/dr**2
        endif
c
        if(j.eq.1)then
          thzza=dum6a(th,dz,i,j,k)
          thz=dum14(th,dz,i,j,k)
          cnz=-2.5d0/dz**2
        elseif(j.eq.nz)then
          thzza=dum8a(th,dz,i,j,k)
          thz=dum16(th,dz,i,j,k)
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
          thzza=dum2a(th,dz,i,j,k)
          thz=dum4(th,dz,i,j,k)
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(w*dz/Ptl)
        if(Pe.gt.2)then
c          print*,'Pe=',Pe
         if(w.gt.0)then
           thz=(th(i,j,k)-th(i,j-1,k))/dz
         elseif(w.lt.0)then
           thz=(th(i,j+1,k)-th(i,j,k))/dz
         endif
        endif
c-----------------hybrid scheme ends for thz in advectice term---------c
          cnz=-1.d0/dz**2
        endif
c
        if(k.eq.1)then
          thhha=(3.2d0*th(i,j,0)+2.d0*th(i,j,2)-0.2d0*th(i,j,3))/dp**2
          thh=(-4.d0*th(i,j,0)/3.d0+th(i,j,1)+th(i,j,2)/3.d0)/dp
          cnp=-2.5d0/dp**2/r**2
        elseif(k.eq.np)then
          thhha=(-0.2d0*th(i,j,np-2)+2.d0*th(i,j,np-1)+3.2d0*th(i,
     +          j,np+1))/dp**2
          thz=(-th(i,j,np-1)/3.d0-th(i,j,np)+4.d0*th(i,j,np+1)/3.d0)/dp
          cnp=-2.5d0/dp**2/r**2
        elseif(k.eq.0)then
          thhha=8.d0*th(i,j,1)/(dp**2)
          thh=0.d0
          cnp=-4.d0/(dp**2)/(r**2)
        elseif(k.eq.np+1)then
          thhha=8.d0*th(i,j,np)/(dp**2)
          thh=0.d0
          cnp=-4.d0/(dp**2)/(r**2)
        else
          thhha=(th(i,j,ka)+th(i,j,km))/dp**2
          thh=(th(i,j,ka)-th(i,j,km))/(2.d0*dp)
c-------------------hybrid scheme introduced by mpk--------------------c
        Pe=dabs(v*dp*r/Ptl)
        if(Pe.gt.2)then
c         print*,'Pe=',Pe
         if(v.gt.0)then
           thh=(th(i,j,k)-th(i,j,km))/dp
         elseif(v.lt.0)then
           thh=(th(i,j,ka)-th(i,j,k))/dp
         endif
        endif
c-----------------hybrid scheme ends for thh in advectice term---------c
         cnp=-1.d0/(dp**2)/(r**2)
        endif

c
        RHS1=thrra+thra/r+thhha/r**2+thzza+th(i,j,k)/dt  
        ULT=u*thr+v*thh/r+w*thz
        cf1=1.d0/(dt)-cnr-cnz-cnp
        cf2=cnr+cnz+cnp
        thn(i,j,k)=(cf2*thp(i,j,k)-ULT+RHS1)/cf1
c
  10  continue
c      forcing both radial walls to have the same temp distribution
c      switched on
!mpk      do 11 i=1,nr
!mpk      do 11 j=0,nz+1
!mpk         thn(i,j,0)=0.5d0*(thn(i,j,0)+thn(i,j,np+1))
!mpk         thn(i,j,np+1)=thn(i,j,0)
!mpk  11  continue         
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
      INCLUDE "sect360mg.ina"
c minus real*8 An(0:kr,0:kz,0:kp),A(0:kr,0:kz,0:kp),Ap(0:kr,0:kz,0:kp)
c       real*8 Bn(0:kr,0:kz,0:kp),B(0:kr,0:kz,0:kp),Bp(0:kr,0:kz,0:kp)
c       real*8 Cn(0:kr,0:kz,0:kp),C(0:kr,0:kz,0:kp),Cp(0:kr,0:kz,0:kp)
c       real*8 uel(0:kr,0:kz,0:kp),vel(0:kr,0:kz,0:kp),wel(0:kr,0:kz,0:kp)
c       real*8 urh(0:kr,0:kz,0:kp),vrh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
c       common/val/Ra,Re,Ptl,Gr
c       common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
c       common/sol3/uel,vel,wel
c       common/radius/r0,rmax
c
      trans1=0.d0
      do 10 k=0,np+1
      trans=dz*(-184.d0*th(0,0,k)+225.d0*th(1,0,k)
     +     -50.d0*th(2,0,k)+9.d0*th(3,0,k))/(240.d0*dr)
      trans=trans+dz*(-184.d0*th(0,1,k)+225.d0*th(1,1,k)
     +     -50.d0*th(2,1,k)+9.d0*th(3,1,k))/(80.d0*dr)
c
      do 20 j=2,nz-1
        trans=trans+dz*(-184.d0*th(0,j,k)+225.d0*th(1,j,k)-50.d0
     +       *th(2,j,k)+9.d0*th(3,j,k))/(60.d0*dr)
  20  continue
c
      trans=trans+dz*(-184.d0*th(0,nz,k)+225.d0*th(1,nz,k)
     +     -50.d0*th(2,nz,k)+9.d0*th(3,nz,k))/(80.d0*dr)
      trans=trans+dz*(-184.d0*th(0,nz+1,k)+225.d0*th(1,nz+1,k)
     +     -50.d0*th(2,nz+1,k)+9.d0*th(3,nz+1,k))/(240.d0*dr)
c
!mpk      if((k.eq.1).or.(k.eq.np+1))then
!mpk         trans1=trans1+0.25d0*(trans+transp)*dp
!mpk      else
!mpk         trans1=trans1+0.5d0*(trans+transp)*dp
!mpk      endif
!mpk      transp=trans 
      trans1=trans1+trans*dp 
  10  continue
!mpk      trans1=8.d0*trans1
c
      return
      end
c
c----------------------------------------------------------------------c
c     subroutine INT                                                   c
c     ==============                                                   c
c     Calculates the volume integrals                                  c
c----------------------------------------------------------------------c
c
c!mpk      subroutine aaaaaaaa(kt)
c
c!mpk      implicit real*8(a-h,o-z)
c!mpk      parameter(kr=65,kz=65,kp=160,mgen=100000)
c!mpk      real*8 thn(0:kr,0:kz,0:kp),th(0:kr,0:kz,0:kp),thp(0:kr,0:kz,0:kp)
c!mpk      real*8 ssumc(6,mgen),ssums(6,mgen)
c!mpk      common/grals/ssumc,ssums
c!mpk      common/sol2/thn,th,thp
c!mpk      common/num/dt,dr,dz,dp,nr,nz,np
c
c!mpk      do 10 n=1,6
c!mpk        ssumc(n,kt)=0.d0
c!mpk        ssums(n,kt)=0.d0 
c!mpk        sumc=0.d0
c!mpk        sums=0.d0
c!mpk        do 20 k=1,np
c!mpk          ang=dble(k)*dp
c!mpk          do 30 i=0,nr+1c
c!mpk            r=r0+dble(i)*dr
c!mpk          do 30 j=0,nz+1
c!mpk            fact=1.d0
c!mpk            if(i.eq.0.or.i.eq.nr+1)fact=fact/2.d0
c!mpk            if(j.eq.0.or.j.eq.nz+1)fact=fact/2.d0
c!mpk            fc=r*th(i,j,k)*dcos(dble(n)*ang)
c!mpk            fs=r*th(i,j,k)*dsin(dble(n)*ang)
c!mpk            sumc=sumc+fc*fact
c!mpk            sums=sums+fs*fact
c!mpk   30     continue
c!mpk          sumc=sumc*dr*dz
c!mpk          sums=sums*dr*dz
c
c!mpk          ssumc(n,kt)=ssumc(n,kt)+sumc*dp
c!mpk          ssums(n,kt)=ssums(n,kt)+sums*dp
c!mpk   20   continue
c!mpk   10 continue
c
c!mpk      return
c!mpk      end
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
      INCLUDE "sect360mg.ina"
c minus real*8 urh(0:kr,0:kz,0:kp),vrh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
c       common/val/Ra,Re,Ptl,Gr
c
c------left and right side walls (omitting corners & edges)
      do 10 i=1,nr
      do 10 k=1,np-1
        An(i,0,k)=-(3.d0*vel(i,1,k)-vel(i,2,k)/3.d0)/dz
        An(i,nz,k)=-(vel(i,nz-1,k)/3.d0-3.d0*vel(i,nz,k))/dz
  10  continue
      do 11 i=1,nr-1
      do 11 k=1,np
        Bn(i,0,k)=2.d0*uel(i,1,k)/dz
        Bn(i,nz,k)=-2.d0*uel(i,nz,k)/dz
  11  continue
c------inner and outer cylindrical walls (omitting corners & edges)
      do 12 j=1,nz-1
      do 12 k=1,np
        Bn(0,j,k)=-2.d0*wel(1,j,k)/dr
        Bn(nr,j,k)=2.d0*wel(nr,j,k)/dr
  12  continue
      do 13 j=1,nz
      do 13 k=1,np-1
        Cn(0,j,k)=(3.d0*vel(1,j,k)-vel(2,j,k)/3.d0)/dr
        Cn(nr,j,k)=(vel(nr-1,j,k)/3.d0-3.d0*vel(nr,j,k))/dr
  13  continue    
c------radial barriers (omitting corners & edges)
      do 14 i=1,nr
      do 14 j=1,nz
        An(i,j,0)=1.d0/(r0+(dble(i)-0.5d0)*dr)*(2.d0*wel(i,j,1))/dp
        An(i,j,np)=1.d0/(r0+(dble(i)-0.5d0)*dr)*(-2.d0*wel(i,j,np))/dp
  14  continue  
      do 15 i=1,nr-1
      do 15 j=1,nz    
        Cn(i,j,0)=1.d0/(r0+(dble(i)-0.5d0)*dr)*(-2.d0*uel(i,j,1))/dp
        Cn(i,j,np)=1.d0/(r0+(dble(i)-0.5d0)*dr)*(2.d0*uel(i,j,np))/dp
  15  continue
c------
      do 16 j=0,nz+1
      do 16 k=0,np+1
       thn(0,j,k)=0.d0
       thn(nr+1,j,k)=1.d0
  16  continue
c
      do 40 i=0,nr+1
      do 40 j=0,nz+1
      do 40 k=0,np+1
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
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine rhsw(wrh)
c
      INCLUDE "sect360mg.ina"
c minus real*8 urh(0:kr,0:kz,0:kp),vrh(0:kr,0:kz,0:kp)
c       common/val/Ra,Re,Ptl,Gr
c
      do 10 i=1,nr
      do 10 j=1,nz-1
      do 10 k=1,np
        km=k-1
        r=r0+(dble(i)-0.5d0)*dr
        Ah=(An(i,j,k)-An(i,j,km))/dp
        Br=(Bn(i,j,k)-Bn(i-1,j,k))/dr
        if(i.eq.1)then
          wrh(i,j,k)=-Bn(1,j,k)/dr-0.5d0*Bn(1,j,k)/r
        elseif(i.eq.nr)then
          wrh(i,j,k)=Bn(nr-1,j,k)/dr-0.5d0*Bn(nr-1,j,k)/r
        else
          wrh(i,j,k)=-Br-0.5d0*(Bn(i,j,k)+Bn(i-1,j,k))/r
        endif
c
        if(k.eq.1)then
          wrh(i,j,k)=wrh(i,j,k)+An(i,j,k)/r/dp
        elseif(k.eq.np)then
          wrh(i,j,k)=wrh(i,j,k)-An(i,j,km)/r/dp
        else
          wrh(i,j,k)=wrh(i,j,k)+Ah/r
        endif
c    
   10 continue
c
      return 
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine rhsv(vrh)
c
      INCLUDE "sect360mg.ina"
c minus real*8 urh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
c       common/val/Ra,Re,Ptl,Gr
c
      do 10 i=1,nr
      do 10 j=1,nz
      do 10 k=1,np-1
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
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine rhsu(urh)
c
      INCLUDE "sect360mg.ina"
c minus real*8 vrh(0:kr,0:kz,0:kp),wrh(0:kr,0:kz,0:kp)
c       common/val/Ra,Re,Ptl,Gr
c
      do 10 i=1,nr-1
      do 10 j=1,nz
      do 10 k=1,np
        km=k-1
        r=r0+dble(i)*dr
        Bz=(Bn(i,j,k)-Bn(i,j-1,k))/dz
        welz=0.5d0*(wel(i,j,k)+wel(i+1,j,k)-wel(i,j-1,k)
     +      -wel(i+1,j-1,k))/dz
        Ch=(Cn(i,j,k)-Cn(i,j,km))/dp
c
        if(j.eq.1)then
          urh(i,j,k)=Bn(i,1,k)/dz-2.d0*welz/r
        elseif(j.eq.nz)then
          urh(i,j,k)=-Bn(i,nz-1,k)/dz-2.d0*welz/r
        else
          urh(i,j,k)=Bz-2.d0*welz/r
        endif
c
        if(k.eq.1)then
          urh(i,j,k)=urh(i,j,k)-Cn(i,j,k)/r/dp
        elseif(k.eq.np)then
          urh(i,j,k)=urh(i,j,k)+Cn(i,j,k-1)/r/dp
        else
          urh(i,j,k)=urh(i,j,k)-Ch/r
        endif
   10 continue
c
      return
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      real*8 function dum1a(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum1a=(v(i-1,j,k)+v(i+1,j,k))/dx**2
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum2a(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum2a=(v(i,j-1,k)+v(i,j+1,k))/dx**2
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum3(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum3=(v(i+1,j,k)-v(i-1,j,k))*0.5d0/dx
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum4(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum4=(v(i,j+1,k)-v(i,j-1,k))*0.5d0/dx
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum5a(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum5a=(3.2d0*v(i-1,j,k)+2.d0*v(i+1,j,k)-0.2d0*v(i+2,j,k))/dx**2
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum6a(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum6a=(3.2d0*v(i,j-1,k)+2.d0*v(i,j+1,k)-0.2d0*v(i,j+2,k))/dx**2
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum7a(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum7a=(3.2d0*v(i+1,j,k)+2.d0*v(i-1,j,k)-0.2d0*v(i-2,j,k))/dx**2
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum8a(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum8a=(3.2d0*v(i,j+1,k)+2.d0*v(i,j-1,k)-0.2d0*v(i,j-2,k))/dx**2
c
      return
      end
c----------------------------------------------------------------------c
c!mpk      real*8 function dum9(v,dx,i,j,k)
c!mpkc
c!mpk      implicit real*8(a-h,o-z)
c!mpk      parameter(kr=33,kz=33,kp=161)
c!mpk      real*8 v(0:kr,0:kz,0:kp)
c!mpkc
c!mpk      dum9=(-8.d0*v(i,j,k)/3.d0+3.d0*v(i+1,j,k)-v(i+2,j,k)/3.d0)/dx
c!mpkc
c!mpk      return
c!mpk      end
c----------------------------------------------------------------------c
c!mpk      real*8 function dum10(v,dx,i,j,k)
c
c!mpk      implicit real*8(a-h,o-z)
c!mpk      parameter(kr=65,kz=65,kp=161)
c!mpk      real*8 v(0:kr,0:kz,0:kp)
c
c!mpk      dum10=(-8.d0*v(i,j,k)/3.d0+3.d0*v(i,j+1,k)-v(i,j+2,k)/3.d0)/dx
c
c!mpk      return
c!mpk      end
c!mpkc----------------------------------------------------------------------c
c!mpk      real*8 function dum11(v,dx,i,j,k)
c!mpkc
c!mpk      implicit real*8(a-h,o-z)
c!mpk      parameter(kr=65,kz=65,kp=161)
c!mpk      real*8 v(0:kr,0:kz,0:kp)
c!mpkc
c!mpk      dum11=(v(i-2,j,k)/3.d0-3.d0*v(i-1,j,k)+8.d0*v(i,j,k)/3.d0)/dx
c!mpkc
c!mpk      return
c!mpk      end
c!mpkc----------------------------------------------------------------------c
c!mpk      real*8 function dum12(v,dx,i,j,k)
c!mpkc
c!mpk      implicit real*8(a-h,o-z)
c!mpk      parameter(kr=65,kz=65,kp=161)
c!mpk      real*8 v(0:kr,0:kz,0:kp)
c!mpkc
c!mpk      dum12=(v(i,j-2,k)/3.d0-3.d0*v(i,j-1,k)+8.d0*v(i,j,k)/3.d0)/dx
c!mpkc
c!mpk      return
c!mpk      end
c----------------------------------------------------------------------c
      real*8 function dum13(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum13=(-4.d0*v(i-1,j,k)/3.d0+v(i,j,k)+v(i+1,j,k)/3.d0)/dx
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum13a(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum13a=(-4.d0*v(i-1,j,k)/3.d0+v(i+1,j,k)/3.d0)/dx
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum14(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum14=(-4.d0*v(i,j-1,k)/3.d0+v(i,j,k)+v(i,j+1,k)/3.d0)/dx
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum15(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum15=(-v(i-1,j,k)/3.d0-v(i,j,k)+4.d0*v(i+1,j,k)/3.d0)/dx
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum15a(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum15a=(-v(i-1,j,k)/3.d0+4.d0*v(i+1,j,k)/3.d0)/dx
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum16(v,dx,i,j,k)
c
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161)
      real*8 v(0:kr,0:kz,0:kp)
c
      dum16=(-v(i,j-1,k)/3.d0-v(i,j,k)+4.d0*v(i,j+1,k)/3.d0)/dx
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine poissoninit(nr,nz,np,dr,dz,dp,nlev)
c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c-----------------------------------------------------------------------
c
      INCLUDE "sect360mg.inb"
c minus real*8 psi(0:kr,0:kz,0:kp),rhs(0:kr,0:kz,0:kp)
c       real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
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
      if(nlev.gt.4)nlev=4
c
      ii=1
      do 20 i=1,nlev
        n1(i)=nr*2/(2**i)
        n2(i)=nz*2/(2**i)
        n3(i)=np*2/(2**i)
        hr(i)=dr*dble(2**i)/2.d0
        hz(i)=dz*dble(2**i)/2.d0
        hp(i)=dp*dble(2**i)/2.d0
        hhr(i)=hr(i)**2
        hhz(i)=hz(i)**2
        hhp(i)=hp(i)**2
!mpk        print*,ii
        ibeg(i)=ii
        ii=ii+(n1(i)+1)*(n2(i)+1)*(n3(i)+1)+1
!mpk        print*,ii
c       print*,i
c       print*,n1(i),n2(i),n3(i),ibeg(i)
c       print*,hr(i),hz(i),hp(i)
   20 continue
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine poisson1(psi,psip,rhs,nr,nz,np,dr,dz,dp,nlev)
c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c-----------------------------------------------------------------------
c
      INCLUDE "sect360mg.inb"
c minus real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
c
      nvs=10000
      nrel1=2
      nrel2=2
c
      do 10 i=1,mgen2
        v(i)=0.d0
        f(i)=0.d0
   10 continue
c
      do 30 i=0,nr+1
      do 30 j=0,nz
      do 30 k=0,np+1
        psic(i,j,k)=0.d0
!mpk        psif(i,j,k)=psi(i,j,k)
        psif(i,j,k)=2.d0*psi(i,j,k)-psip(i,j,k)
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
        do 31 ks=0,np+1
        psif(is,js,ks)=psip(is,js,ks)+1.15d0*(psif(is,js,ks)-
     +    psip(is,js,ks))
        psip(is,js,ks)=psif(is,js,ks)
   31   continue
        call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
        endif
        if(i.eq.1.or.nlev.eq.1)then
          call relax1(nrel1,1,nconv,nlev)
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser1(ilev)
            call relax1(nrel1,ilev+1,nconv,nlev)
   51     continue
          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax1(nrel2,ilev-1,nconv,nlev)
   52     continue
        end if
        if(nconv.eq.1)then
c         print*,'Converged'
          do 59 ii=0,nr+1
          do 59 jj=0,nz
          do 59 kk=0,np+1
            psip(ii,jj,kk)=psi(ii,jj,kk)
   59     continue
c       internal points updated
          do 60 ii=1,nr
          do 60 jj=1,nz-1
          do 60 kk=1,np
            psi(ii,jj,kk)=psif(ii,jj,kk)
   60     continue
c       endwalls
          do 61 ii=0,nr+1
          do 61 kk=0,np+1
            psi(ii,0,kk)=0.d0
            psi(ii,nz,kk)=0.d0
   61     continue
c       cylinders
          do 62 jj=1,nz-1
          do 62 kk=0,np+1
            if(kk.eq.0)then
              psi(0,jj,0)=-psi(1,jj,1)
              psi(nr+1,jj,0)=-psi(nr,jj,1)
            elseif(kk.eq.np+1)then
              psi(0,jj,np+1)=-psi(1,jj,np)
              psi(nr+1,jj,np+1)=-psi(nr,jj,np)
            else
              psi(0,jj,kk)=-psi(1,jj,kk)
              psi(nr+1,jj,kk)=-psi(nr,jj,kk)
            endif
   62     continue
c       radial barriers
          do 63 ii=1,nr
          do 63 jj=1,nz-1
            psi(ii,jj,0)=-psi(ii,jj,1)
            psi(ii,jj,np+1)=-psi(ii,jj,np)
   63     continue
          print*,'w vel; sweeps=',i
c         call rufcon(psi,nr,nz,np)
          return
        end if
   50 continue
      return
      end
c
c-----------------------------------------------------------------------
      subroutine poisson2(psi,psip,rhs,nr,nz,np,dr,dz,dp,nlev)
c
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c-----------------------------------------------------------------------
c
      INCLUDE "sect360mg.inb"
c minus real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
c
      nvs=10000
      nrel1=2
      nrel2=2
c
      do 10 i=1,mgen2
        v(i)=0.d0
        f(i)=0.d0
   10 continue
c
      do 30 i=0,nr+1
      do 30 j=0,nz+1
      do 30 k=0,np
        psic(i,j,k)=0.d0
!mpk        psif(i,j,k)=psi(i,j,k)
        psif(i,j,k)=2.d0*psi(i,j,k)-psip(i,j,k)
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
        psif(is,js,ks)=psip(is,js,ks)+1.15d0*(psif(is,js,ks)-
     +   psip(is,js,ks))
        psip(is,js,ks)=psif(is,js,ks)
   31   continue
        call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
      endif
        if(i.eq.1.or.nlev.eq.1)then
          call relax2(nrel1,1,nconv,nlev)
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser2(ilev)
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
c       internal points
          do 60 ii=1,nr
          do 60 jj=1,nz
          do 60 kk=1,np-1
            psi(ii,jj,kk)=psif(ii,jj,kk)
   60     continue
c       radial barriers
          do 61 ii=0,nr+1
          do 61 jj=0,nz+1
            psi(ii,jj,0)=0.d0
            psi(ii,jj,np)=0.d0
   61     continue
c       cylinders
          do 62 jj=0,nz+1
          do 62 kk=1,np-1
            psi(0,jj,kk)=0.d0        
            psi(nr+1,jj,kk)=0.d0
   62     continue
c       endwalls
          do 63 ii=1,nr
          do 63 kk=1,np-1
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
c-----------------------------------------------------------------------
c
      INCLUDE "sect360mg.inb"
c minus real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
c
      nvs=10000
      nrel1=2
      nrel2=2
c
      do 10 i=1,mgen2
        v(i)=0.d0
        f(i)=0.d0
   10 continue
c
      do 30 i=0,nr
      do 30 j=0,nz+1
      do 30 k=0,np+1
        psic(i,j,k)=0.d0
!mpk        psif(i,j,k)=psi(i,j,k)
        psif(i,j,k)=2.d0*psi(i,j,k)-psip(i,j,k)
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
        do 31 ks=0,np+1
        psif(is,js,ks)=psip(is,js,ks)+1.15d0*(psif(is,js,ks)-
     +    psip(is,js,ks))
        psip(is,js,ks)=psif(is,js,ks)
   31   continue
        call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
       endif
        if(i.eq.1.or.nlev.eq.1)then
          call relax3(nrel1,1,nconv,nlev)
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser3(ilev)
            call relax3(nrel1,ilev+1,nconv,nlev)
   51     continue
          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax3(nrel2,ilev-1,nconv,nlev)
   52     continue
        end if
        if(nconv.eq.1)then
c         print*,'Converged'
          do 59 ii=0,nr
          do 59 jj=0,nz+1
          do 59 kk=0,np+1
            psip(ii,jj,kk)=psi(ii,jj,kk)
   59     continue
          do 60 ii=1,nr-1
          do 60 jj=1,nz
          do 60 kk=1,np
            psi(ii,jj,kk)=psif(ii,jj,kk)
   60     continue
c       cylinders
          do 61 jj=0,nz+1
          do 61 kk=0,np+1
            psi(0,jj,kk)=0.d0
            psi(nr,jj,kk)=0.d0
   61     continue
c       endwalls
          do 62 ii=1,nr-1
          do 62 kk=0,np+1
            if(kk.eq.0)then
             psi(ii,0,0)=-psi(ii,1,1)
             psi(ii,nz+1,0)=-psi(ii,nz,1)
            elseif(kk.eq.np+1)then
             psi(ii,0,np+1)=-psi(ii,1,np)
             psi(ii,nz+1,np+1)=-psi(ii,nz,np)
            else
             psi(ii,0,kk)=-psi(ii,1,kk)
             psi(ii,nz+1,kk)=-psi(ii,nz,kk)
            endif
   62     continue
c       radial barriers
          do 63 ii=1,nr-1
          do 63 jj=1,nz
            psi(ii,jj,0)=-psi(ii,jj,1)
            psi(ii,jj,np+1)=-psi(ii,jj,np)
   63     continue
          print*,'u vel; sweeps=',i
c         call rufcon(psi,nr,nz,np)
          return
        end if
   50 continue
      return
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine relax1(nrel,klev,nconv,nlev)
c
      INCLUDE "sect360mg.inb"
c minus real*8 psi(0:kr,0:kz,0:kp),rhs(0:kr,0:kz,0:kp)
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
        do 10 j=1,n2(klev)-1
        do 10 k=1,n3(klev)
          do 11 i=1,n1(klev)
          rval=r0+(dble(i)-0.5d0)*hr(klev)
c
            if(i.eq.1)then
              tdmaa(i)=0.d0
              tdmab(i)=-1.d0/hhr(klev)-0.5d0/hr(klev)/rval-
     +                 2.d0/hhz(klev)
              tdmac(i)=1.d0/hhr(klev)+0.5d0/hr(klev)/rval
              tdmar(i)=-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)+
     +                 rf(i,j,k)
            elseif(i.eq.n1(klev))then
              tdmaa(i)=1.d0/hhr(klev)-0.5d0/hr(klev)/rval
              tdmab(i)=-1.d0/hhr(klev)+0.5d0/hr(klev)/rval-
     +                 2.d0/hhz(klev)
              tdmac(i)=0.d0
              tdmar(i)=-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)+
     +                rf(i,j,k)
            else
              tdmaa(i)=1.d0/hhr(klev)-0.5d0/hr(klev)/rval
              tdmab(i)=-(2.d0/hhr(klev)+2.d0/hhz(klev))
              tdmac(i)=1.d0/hhr(klev)+0.5d0/hr(klev)/rval
             
              tdmar(i)=-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)+
     +                 rf(i,j,k)
            endif
c
            if(k.eq.1)then
              tdmab(i)=tdmab(i)-1.d0/hhp(klev)/rval**2
              tdmar(i)=tdmar(i)-psif(i,j,2)/hhp(klev)/rval**2
            elseif(k.eq.n3(klev))then
              tdmab(i)=tdmab(i)-1.d0/hhp(klev)/rval**2
              tdmar(i)=tdmar(i)-psif(i,j,k-1)/hhp(klev)/rval**2
            else
              tdmab(i)=tdmab(i)-2.d0/hhp(klev)/rval**2
              tdmar(i)=tdmar(i)-(psif(i,j,k+1)+psif(i,j,k-1))/
     +                 hhp(klev)/rval**2
            endif
   11     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n1(klev))
          do 12 i=1,n1(klev)
            psif(i,j,k)=tdmas(i)
   12     continue
   10   continue
c                               line relax in z-dir
        do 13 i=1,n1(klev)
           rval=r0+(dble(i)-0.5d0)*hr(klev)
        do 13 k=1,n3(klev)
          do 14 j=1,n2(klev)-1
            tdmaa(j)=1.d0/hhz(klev)
            tdmac(j)=1.d0/hhz(klev)
            if(j.eq.1)then
              tdmaa(j)=0.d0
            elseif(j.eq.n2(klev)-1)then
              tdmac(j)=0.d0
            endif
c            
            if(i.eq.1)then
              tdmab(j)=-2.d0/hhz(klev)-1.d0/hhr(klev)-0.5d0/hr(klev)/
     +                  rval
              tdmar(j)=-(psif(2,j,k)/hhr(klev)+0.5d0*psif(2,j,k)/
     +                  hr(klev)/rval)+rf(i,j,k)
            elseif(i.eq.n1(klev))then
              tdmab(j)=-2.d0/hhz(klev)-1.d0/hhr(klev)+0.5d0/hr(klev)/
     +                 rval
              tdmar(j)=-(psif(i-1,j,k)/hhr(klev)-0.5d0*psif(i-1,j,k)/
     +                 hr(klev)/rval)+rf(i,j,k)
            else
              tdmab(j)=-2.d0/hhz(klev)-2.d0/hhr(klev)
              tdmar(j)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +                  -0.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)/
     +                  rval+rf(i,j,k)
            endif
c
            if(k.eq.1)then
              tdmab(j)=tdmab(j)-1.d0/hhp(klev)/rval**2
              tdmar(j)=tdmar(j)-psif(i,j,2)/hhp(klev)/rval**2
            elseif(k.eq.n3(klev))then
              tdmab(j)=tdmab(j)-1.d0/hhp(klev)/rval**2
              tdmar(j)=tdmar(j)-psif(i,j,k-1)/hhp(klev)/rval**2
            else
              tdmab(j)=tdmab(j)-2.d0/hhp(klev)/rval**2
              tdmar(j)=tdmar(j)-(psif(i,j,k+1)+psif(i,j,k-1)
     +                 )/hhp(klev)/rval**2 
            endif
   14     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n2(klev)-1)
          do 15 j=1,n2(klev)-1
            psif(i,j,k)=tdmas(j)
   15     continue
   13   continue
c                             line relax in phi-dir
        do 16 i=1,n1(klev)
          rval=r0+(dble(i)-0.5d0)*hr(klev)
        do 16 j=1,n2(klev)-1
          do 17 k=1,n3(klev)
c
            if(k.eq.1)then
              tdmaa(k)=0.d0
              tdmab(k)=-(2.d0/hhz(klev)+1.d0/hhp(klev)/rval**2)
              tdmac(k)=1.d0/hhp(klev)/rval**2
              tdmar(k)=-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)+
     +                 rf(i,j,k)
            elseif(k.eq.n3(klev))then
              tdmaa(k)=1.d0/hhp(klev)/rval**2
              tdmab(k)=-(2.d0/hhz(klev)+1.d0/hhp(klev)/rval**2)
              tdmac(k)=0.d0
              tdmar(k)=-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)+
     +                 rf(i,j,k)
            else
              tdmaa(k)=1.d0/hhp(klev)/rval**2
              tdmab(k)=-(2.d0/hhz(klev)+2.d0/hhp(klev)/rval**2)
              tdmac(k)=1.d0/hhp(klev)/rval**2
              tdmar(k)=-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)+
     +                 rf(i,j,k)
            endif
c
            if(i.eq.1)then
              tdmab(k)=tdmab(k)-0.5d0/hr(klev)/rval-1.d0/hhr(klev)
              tdmar(k)=tdmar(k)-psif(i+1,j,k)/hhr(klev)-0.5d0*
     +                 psif(i+1,j,k)/hr(klev)/rval
            elseif(i.eq.n1(klev))then
              tdmab(k)=tdmab(k)+0.5d0/hr(klev)/rval-1.d0/hhr(klev)
              tdmar(k)=tdmar(k)-psif(i-1,j,k)/hhr(klev)+0.5d0*
     +                 psif(i-1,j,k)/hr(klev)/rval
            else
              tdmab(k)=tdmab(k)-2.d0/hhr(klev)
              tdmar(k)=tdmar(k)-(psif(i+1,j,k)+psif(i-1,j,k))/
     +                 hhr(klev)-0.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/
     +                 hr(klev)/rval           
            endif
   17     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n3(klev))
          do 18 k=1,n3(klev)
              psif(i,j,k)=tdmas(k)
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
c 
         if(i.eq.1)then
            temp=(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)-
     +         psif(i,j,k)*2.d0/hhz(klev)+(-psif(i,j,k)+psif(i+1,j,k))/
     +         hhr(klev)+0.5d0*(-psif(i,j,k)+psif(i+1,j,k))/hr(klev)/
     +         rval
         elseif(i.eq.n1(klev))then
            temp=(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)-
     +        psif(i,j,k)*2.d0/hhz(klev)+(-psif(i,j,k)+psif(i-1,j,k))/
     +        hhr(klev)+0.5d0*(psif(i,j,k)-psif(i-1,j,k))/hr(klev)/
     +        rval
         else
            temp=(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)-
     +        psif(i,j,k)*2.d0/hhz(klev)+(psif(i-1,j,k)-
     +        2.d0*psif(i,j,k)+psif(i+1,j,k))/hhr(klev)+
     +        0.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)/rval
         endif  
c
         if(k.eq.1)then
           temp=temp+(-psif(i,j,k)+psif(i,j,k+1))/hhp(klev)/rval**2
         elseif(k.eq.n3(klev))then
           temp=temp+(psif(i,j,k-1)-psif(i,j,k))/hhp(klev)/rval**2
         else
           temp=temp+(psif(i,j,k-1)+psif(i,j,k+1))/hhp(klev)/
     +          rval**2-psif(i,j,k)*2.d0/hhp(klev)/rval**2
         endif
c
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
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine relax2(nrel,klev,nconv,nlev)
c
      INCLUDE "sect360mg.inb"
      real*8 pdmaa(kp),pdmab(kp),pdmac(kp),pdmad(kp),pdmae(kp),pdmas(kp)
      real*8 pdmar(kp)
c minus real*8 psi(0:kr,0:kz,0:kp),rhs(0:kr,0:kz,0:kp)
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
        do 10 j=1,n2(klev)
        do 10 k=1,n3(klev)-1
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
     +               /hhz(klev)-(psif(i,j,k+1)+psif(i,j,k-1))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            elseif(j.eq.n2(klev))then
             pdmac(i)=pdmac(i)-2.d0/hhz(klev)
             pdmar(i)=-(5.d0/3.d0*psif(i,j-1,k)-0.2*psif(i,j-2,k))
     +               /hhz(klev)-(psif(i,j,k+1)+psif(i,j,k-1))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            else
             pdmac(i)=pdmac(i)-2.d0/hhz(klev)
             pdmar(i)=-(psif(i,j-1,k)+psif(i,j+1,k))
     +               /hhz(klev)-(psif(i,j,k+1)+psif(i,j,k-1))/hhp(klev)
     +               /rval**2+rf(i,j,k)
            endif
c
   11     continue
          call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n1(klev))
          do 12 i=1,n1(klev)
            psif(i,j,k)=pdmas(i)
   12     continue
   10   continue
c                                 line relax in z-dir
        do 13 i=1,n1(klev)
            rval=r0+(dble(i)-0.5d0)*hr(klev)
        do 13 k=1,n3(klev)-1
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
     +                -(psif(i,j,k+1)+psif(i,j,k-1))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            elseif(i.eq.n1(klev))then
            pdmac(j)=pdmac(j)-2.d0/hhr(klev)      !mpk -1.d0/rval/hr(klev)
            pdmar(j)=-(5.d0/3.d0*psif(i-1,j,k)-0.2d0*psif(i-2,j,k))/
     +                hhr(klev)+4.d0/3.d0*psif(i-1,j,k)/hr(klev)/rval
     +                -(psif(i,j,k+1)+psif(i,j,k-1))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            else
              pdmac(j)=pdmac(j)-2.d0/hhr(klev)
              pdmar(j)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)-1.5d0*
     +                (psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)/rval
     +                -(psif(i,j,k+1)+psif(i,j,k-1))/hhp(klev)/rval**2
     +                +rf(i,j,k)
            endif
   14     continue
          call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n2(klev))
          do 15 j=1,n2(klev)
            psif(i,j,k)=pdmas(j)
   15     continue
   13   continue
c                             line relax in phi-dir
        do 16 i=1,n1(klev)
          rval=r0+(dble(i)-0.5d0)*hr(klev)
        do 16 j=1,n2(klev)
          do 17 k=1,n3(klev)-1
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
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n3(klev)-1)
          do 18 k=1,n3(klev)-1
              psif(i,j,k)=tdmas(k)
   18     continue
   16   continue
c
c                                  Calculating residual.
c
          resmax=0.d0
          do 20 i=1,n1(klev)
             rval=r0+(dble(i)-0.5d0)*hr(klev)
          do 20 j=1,n2(klev)
          do 20 k=1,n3(klev)-1
          temp=(psif(i,j,k+1)-2.d0*psif(i,j,k)+psif(i,j,k-1))
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
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine relax3(nrel,klev,nconv,nlev)
c!mpk      common/iteration/kt
c
      INCLUDE "sect360mg.inb"
c minus real*8 psi(0:kr,0:kz,0:kp),rhs(0:kr,0:kz,0:kp)
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
c                                   line relax in r-dir
        do 10 j=1,n2(klev)
        do 10 k=1,n3(klev)
          do 11 i=1,n1(klev)-1
            rval=dble(i)*hr(klev) + r0
            tdmaa(i)=1.d0/hhr(klev)-1.5d0/hr(klev)/rval
            tdmab(i)=-(2.d0/hhr(klev)-1.d0/rval**2)
            tdmac(i)=1.d0/hhr(klev)+1.5d0/hr(klev)/rval
            if(i.eq.1)then
              tdmaa(i)=0.d0
            elseif(i.eq.n1(klev)-1)then
              tdmac(i)=0.d0
            endif
c            
            if(j.eq.1)then
              tdmab(i)=tdmab(i)-1.d0/hhz(klev)
              tdmar(i)=-psif(i,j+1,k)/hhz(klev)+rf(i,j,k)
            elseif(j.eq.n2(klev))then
              tdmab(i)=tdmab(i)-1.d0/hhz(klev)
              tdmar(i)=-psif(i,j-1,k)/hhz(klev)+rf(i,j,k)
            else
              tdmab(i)=tdmab(i)-2.d0/hhz(klev)
              tdmar(i)=-(psif(i,j-1,k)+psif(i,j+1,k))/hhz(klev)+
     +                 rf(i,j,k)
            endif
c
            if(k.eq.1)then
              tdmab(i)=tdmab(i)-1.d0/hhp(klev)/rval**2
              tdmar(i)=tdmar(i)-psif(i,j,k+1)/hhp(klev)/rval**2
            elseif(k.eq.n3(klev))then
              tdmab(i)=tdmab(i)-1.d0/hhp(klev)/rval**2
              tdmar(i)=tdmar(i)-psif(i,j,k-1)/hhp(klev)/rval**2
            else
              tdmab(i)=tdmab(i)-2.d0/hhp(klev)/rval**2
              tdmar(i)=tdmar(i)-(psif(i,j,k+1)+psif(i,j,k-1))/
     +                 hhp(klev)/rval**2
            endif      
   11     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n1(klev)-1)
          do 12 i=1,n1(klev)-1
            psif(i,j,k)=tdmas(i)
   12     continue
   10   continue
c                                 line relax in z-dir
        do 13 i=1,n1(klev)-1
           rval=dble(i)*hr(klev) + r0
        do 13 k=1,n3(klev)
          do 14 j=1,n2(klev)
            if(j.eq.1)then
              tdmaa(j)=0.d0
              tdmab(j)=-1.d0/hhz(klev)-2.d0/hhr(klev)+1.d0/rval**2
              tdmac(j)=1.d0/hhz(klev)
            elseif(j.eq.n2(klev))then
              tdmaa(j)=1.d0/hhz(klev)
              tdmab(j)=-1.d0/hhz(klev)-2.d0/hhr(klev)+1.d0/rval**2
              tdmac(j)=0.d0
            else
              tdmaa(j)=1.d0/hhz(klev)
              tdmab(j)=-(2.d0/hhz(klev)+2.d0/hhr(klev))+1.d0/
     +                 rval**2
              tdmac(j)=1.d0/hhz(klev)
            endif
            tdmar(j)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +               /rval+rf(i,j,k)
c
            if(k.eq.1)then
              tdmab(j)=tdmab(j)-1.d0/hhp(klev)/rval**2
              tdmar(j)=tdmar(j)-psif(i,j,k+1)/hhp(klev)/rval**2
            elseif(k.eq.n3(klev))then
              tdmab(j)=tdmab(j)-1.d0/hhp(klev)/rval**2
              tdmar(j)=tdmar(j)-psif(i,j,k-1)/hhp(klev)/rval**2
            else
              tdmab(j)=tdmab(j)-2.d0/hhp(klev)/rval**2
              tdmar(j)=tdmar(j)-(psif(i,j,k-1)+psif(i,j,k+1))/
     +                 hhp(klev)/rval**2
            endif
   14     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n2(klev))
          do 15 j=1,n2(klev)
            psif(i,j,k)=tdmas(j)
   15     continue
   13   continue
c                                line relax in phi-dir
        do 16 i=1,n1(klev)-1
           rval=dble(i)*hr(klev) + r0
        do 16 j=1,n2(klev)
          do 17 k=1,n3(klev)
            if(k.eq.1)then
              tdmaa(k)=0.d0
              tdmab(k)=-1.d0/rval**2/hhp(klev)-2.d0/hhr(klev)+
     +                 1/rval**2
              tdmac(k)=1.d0/rval**2/hhp(klev)
            elseif(k.eq.n3(klev))then
              tdmaa(k)=1.d0/rval**2/hhp(klev)
              tdmab(k)=-1.d0/rval**2/hhp(klev)-2.d0/hhr(klev)+
     +                 1/rval**2
              tdmac(k)=0.d0
            else
              tdmaa(k)=1.d0/rval**2/hhp(klev)
              tdmab(k)=-2.d0/rval**2/hhp(klev)-2.d0/hhr(klev)+
     +                 1/rval**2
              tdmac(k)=1.d0/rval**2/hhp(klev)
            endif
c
            if(j.eq.1)then
              tdmab(k)=tdmab(k)-1.d0/hhz(klev)
              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +               /rval-psif(i,j+1,k)/hhz(klev)+rf(i,j,k)
            elseif(j.eq.n2(klev))then
              tdmab(k)=tdmab(k)-1.d0/hhz(klev)
              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +               /rval-psif(i,j-1,k)/hhz(klev)+rf(i,j,k)
            else
              tdmab(k)=tdmab(k)-2.d0/hhz(klev)
              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
     +               /rval-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)
     +               +rf(i,j,k)
            endif
   17     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n3(klev))
          do 18 k=1,n3(klev)
              psif(i,j,k)=tdmas(k)
   18     continue
   16   continue
c
c                                 Calculating residual.
c
          resmax=0.d0
          do 20 i=1,n1(klev)-1
!mpk            ii=i*levfact
          rval=dble(i)*hr(klev) + r0
          do 20 j=1,n2(klev)
          do 20 k=1,n3(klev)
          temp=(psif(i+1,j,k)-psif(i,j,k)*2.d0+psif(i-1,j,k))/
     +         hhr(klev)+1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/
     +         hr(klev)/rval+psif(i,j,k)/rval**2
          if(j.eq.1)then
            temp=temp+(-psif(i,j,k)+psif(i,j+1,k))/hhz(klev)
          elseif(j.eq.n2(klev))then
            temp=temp+(-psif(i,j,k)+psif(i,j-1,k))/hhz(klev)
          else
            temp=temp+(psif(i,j+1,k)-2.d0*psif(i,j,k)+psif(i,j-1,k))
     +          /hhz(klev)
          endif
c
          if(k.eq.1)then
            temp=temp+(-psif(i,j,k)+psif(i,j,k+1))/rval**2/hhp(klev)
          elseif(k.eq.n3(klev))then
            temp=temp+(-psif(i,j,k)+psif(i,j,k-1))/rval**2/hhp(klev)
          else
            temp=temp+(-2.d0*psif(i,j,k)+psif(i,j,k-1)+
     +           psif(i,j,k+1))/rval**2/hhp(klev) 
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
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine corser1(klev)
c
      INCLUDE "sect360mg.inb"
c minus real*8 psi(0:kr,0:kz,0:kp),rhs(0:kr,0:kz,0:kp)
c       real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
c       common /radius/r0,rmax
c
c------------------Start off with restriction
c
      call trans3(rc,rf,n1(klev),n2(klev),n3(klev))
c          at endwalls
      do 3 i=0,n1(klev)+1
      do 3 k=0,n3(klev)+1
         rf(i,0,k)=0.d0
         rf(i,n2(klev),k)=0.d0
  3   continue
c          at cylinders
      do 4 j=1,n2(klev)-1
      do 4 k=0,n3(klev)+1
         if(k.eq.0)then
          rf(n1(klev)+1,j,0)=-rf(n1(klev),j,1)
          rf(0,j,0)=-rf(1,j,1)
         elseif(k.eq.n3(klev))then
          rf(n1(klev)+1,j,k)=-rf(n1(klev),j,k-1)
          rf(0,j,k)=-rf(1,j,k-1)
         else
          rf(n1(klev)+1,j,k)=-rf(n1(klev),j,k)
          rf(0,j,k)=-rf(1,j,k)
         endif
  4   continue
c          at radial barriers
      do 5 i=1,n1(klev)
      do 5 j=1,n2(klev)-1
         rf(i,j,n3(klev)+1)=-rf(i,j,n3(klev))
         rf(i,j,0)=-rf(i,j,1)
  5   continue
c
!      do 10 i=1,n1(klev+1)-1
      do 10 i=1,n1(klev+1)
      do 10 j=1,n2(klev+1)-1
!      do 10 k=1,n3(klev+1)-1
      do 10 k=1,n3(klev+1)
        ii=i+i
        jj=j+j
        kk=k+k
        iip=ii+1
        iim=ii-1
        jjp=jj+1
        jjm=jj-1
        kkp=kk+1
        kkm=kk-1          
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
        temp=temp/64.d0
        rc(i,j,k)=temp
   10 continue
c
c-----------setting faces to zero
c     at endwalls
      do 20 i=0,n1(klev+1)+1
      do 20 k=0,n3(klev+1)+1   
!      do 20 i=0,n1(klev+1)+1
!      do 20 k=0,n3(klev+1)+1
        rc(i,0,k)=0.d0
        rc(i,n2(klev+1),k)=0.d0
        psic(i,0,k)=0.d0
        psic(i,n2(klev+1),k)=0.d0
   20 continue
c     at cylinders
      do 30 j=1,n2(klev+1)-1
      do 30 k=0,n3(klev+1)+1
        if(k.eq.0)then
          rc(0,j,0)=-rc(1,j,1)
          rc(n1(klev+1)+1,j,0)=-rc(n1(klev+1),j,1)
          psic(0,j,0)=-psic(1,j,1)
          psic(n1(klev+1)+1,j,0)=-psic(n1(klev+1),j,1)
        elseif(k.eq.n3(klev+1)+1)then
          rc(0,j,k)=-rc(1,j,k-1)
          rc(n1(klev+1)+1,j,k)=-rc(n1(klev+1),j,k-1)
          psic(0,j,k)=-psic(1,j,k-1)
          psic(n1(klev+1)+1,j,k)=-psic(n1(klev+1),j,k-1)
        else
          rc(0,j,k)=-rc(1,j,k)
          rc(n1(klev+1)+1,j,k)=-rc(n1(klev+1),j,k)
          psic(0,j,k)=-psic(1,j,k)
          psic(n1(klev+1)+1,j,k)=-psic(n1(klev+1),j,k)
        endif
   30 continue
!      do 30 j=0,n2(klev+1)
!      do 30 k=0,n3(klev+1)
!           rc(0,j,k)=0.d0
!           rc(n1(klev+1),j,k)=0.d0
!           psic(0,j,k)=0.d0
!           psic(n1(klev+1),j,k)=0.d0
!   30 continue
c     at radial barriers
      do 40 i=1,n1(klev+1)
      do 40 j=1,n2(klev+1)-1
        rc(i,j,0)=-rc(i,j,1)
        rc(i,j,n3(klev+1)+1)=-rc(i,j,n3(klev+1))
        psic(i,j,0)=-psic(i,j,1)
        psic(i,j,n3(klev+1)+1)=-psic(i,j,n3(klev+1))
   40 continue
!       do 40 i=0,n1(klev+1)
!       do 40 j=0,n2(klev+1)
!         rc(i,j,0)=0.d0
!         rc(i,j,n3(klev+1))=0.d0
!         psic(i,j,0)=0.d0
!         psic(i,j,n3(klev+1))=0.d0
!   40 continue
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
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine corser2(klev)
c
      include "sect360mg.inb"
c minus psi(0:kr,0:kz,0:kp),rhs(0:kr,0:kz,0:kp)
c       tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
c       common /radius/r0,rmax
c
c------------------Start off with restriction
c
      call trans3(rc,rf,n1(klev),n2(klev),n3(klev))
c          radial barriers
      do 2 i=0,n1(klev)+1
      do 2 j=0,n2(klev)+1
          rf(i,j,0)=0.d0
          rf(i,j,n3(klev))=0.d0
  2   continue
c          cylinders
      do 3 j=0,n2(klev)+1
      do 3 k=1,n3(klev)-1
          rf(n1(klev)+1,j,k)=0.d0
          rf(0,j,k)=0.d0
  3   continue
c          endwalls
      do 4 i=1,n1(klev)
      do 4 k=1,n3(klev)-1
         rf(i,n2(klev)+1,k)=0.d0
         rf(i,0,k)=0.d0
  4   continue
c
      do 10 i=1,n1(klev+1)
!      do 10 i=0,n1(klev+1)-1
      do 10 j=1,n2(klev+1)
!      do 10 j=0,n2(klev+1)-1
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
        temp=temp/64.d0
        rc(i,j,k)=temp
   10 continue
c
c-----------setting faces to zero
c     at radial barriers
!      do 20 i=0,n1(klev+1)
      do 20 i=0,n1(klev+1)+1
      do 20 j=0,n2(klev+1)+1
        rc(i,j,0)=0.d0
        rc(i,j,n3(klev+1))=0.d0
        psic(i,j,0)=0.d0
        psic(i,j,n3(klev+1))=0.d0
   20 continue
c     at cylinders
      do 30 j=0,n2(klev+1)+1
      do 30 k=1,n3(klev+1)-1
         rc(0,j,k)=0.d0
         rc(n1(klev+1)+1,j,k)=0.d0
         psic(0,j,k)=0.d0
         psic(n1(klev+1)+1,j,k)=0.d0
   30 continue
!       do 30 j=0,n2(klev+1)
!       do 30 k=0,n3(klev+1)
!          rc(0,j,k)=0.d0
!          rc(n1(klev+1),j,k)=0.d0
!          psic(0,j,k)=0.d0
!          psic(n1(klev+1),j,k)=0.d0
!   30 continue
c     at endwalls   
      do 40 i=1,n1(klev+1)
      do 40 k=1,n3(klev+1)-1
        rc(i,0,k)=0.d0
        rc(i,n2(klev+1)+1,k)=0.d0
        psic(i,0,k)=0.d0
        psic(i,n2(klev+1)+1,k)=0.d0
   40 continue
!       do 40 i=1,n1(klev+1)
!       do 40 k=1,n3(klev+1)-1
!         rc(i,0,k)=0.d0
!         rc(i,n2(klev+1),k)=0.d0
!         psic(i,0,k)=0.d0
!         psic(i,n2(klev+1),k)=0.d0
!   40 continue
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
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine corser3(klev)
c
      INCLUDE "sect360mg.inb"
c minus real*8 psi(0:kr,0:kz,0:kp),rhs(0:kr,0:kz,0:kp)
c       real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
c       common /radius/r0,rmax
c
c------------------Start off with restriction
c
      call trans3(rc,rf,n1(klev),n2(klev),n3(klev))
c          cylinders
      do 2 j=0,n2(klev)+1
      do 2 k=0,n3(klev)+1
         rf(0,j,k)=0.d0
         rf(n1(klev),j,k)=0.d0
  2   continue
c          endwalls
      do 3 i=1,n1(klev)-1
      do 3 k=0,n3(klev)+1
         if(k.eq.0)then
          rf(i,n2(klev)+1,0)=-rf(i,n2(klev),1)
          rf(i,0,0)=-rf(i,1,1)
         elseif(k.eq.n3(klev)+1)then
          rf(i,n2(klev)+1,k)=-rf(i,n2(klev),k-1)
          rf(i,0,k)=-rf(i,1,k-1)
         else
          rf(i,n2(klev)+1,k)=-rf(i,n2(klev),k)
          rf(i,0,k)=-rf(i,1,k)
         endif
  3   continue
c          radial barriers
      do 4 i=1,n1(klev)-1
      do 4 j=1,n2(klev)
         rf(i,j,n3(klev)+1)=-rf(i,j,n3(klev))
         rf(i,j,0)=-rf(i,j,1)
  4   continue
c
      do 10 i=1,n1(klev+1)-1
      do 10 j=1,n2(klev+1)
!      do 10 j=1,n2(klev+1)-1   
      do 10 k=1,n3(klev+1)
        ii=i+i
        jj=j+j
        kk=k+k
        iip=ii+1
        iim=ii-1
        jjp=jj+1
        jjm=jj-1
        kkp=kk+1
        kkm=kk-1          
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
        temp=temp/64.d0
        rc(i,j,k)=temp
   10 continue
c
c-----------setting faces to zero
c     at cylinders
!      do 20 j=0,n2(klev+1)
!      do 20 k=0,n3(klev+1)
      do 20 j=0,n2(klev+1)+1
      do 20 k=0,n3(klev+1)+1
        rc(0,j,k)=0.d0
        rc(n1(klev+1),j,k)=0.d0
        psic(0,j,k)=0.d0
        psic(n1(klev+1),j,k)=0.d0
   20 continue
c     at endwalls   
      do 30 i=1,n1(klev+1)-1
      do 30 k=0,n3(klev+1)+1
        if(k.eq.0)then
         rc(i,0,0)=-rc(i,1,1)
         rc(i,n2(klev+1)+1,0)=-rc(i,n2(klev+1),1)
         psic(i,0,0)=-psic(i,1,1)
         psic(i,n2(klev+1)+1,0)=-psic(i,n2(klev+1),1)
        elseif(k.eq.n3(klev+1)+1)then
         rc(i,0,k)=-rc(i,1,k-1)
         rc(i,n2(klev+1)+1,k)=-rc(i,n2(klev+1),k-1)
         psic(i,0,k)=-psic(i,1,k-1)
         psic(i,n2(klev+1)+1,k)=-psic(i,n2(klev+1),k-1)
        else
         rc(i,0,k)=-rc(i,1,k)
         rc(i,n2(klev+1)+1,k)=-rc(i,n2(klev+1),k)
         psic(i,0,k)=-psic(i,1,k)
         psic(i,n2(klev+1)+1,k)=-psic(i,n2(klev+1),k)
        endif
   30 continue
!      do 30 i=0,n1(klev+1)
!      do 30 k=0,n3(klev+1)
!         rc(i,0,k)=0.d0
!         rc(i,n2(klev+1),k)=0.d0
!         psic(i,0,k)=0.d0
!         psic(i,n2(klev+1),k)=0.d0
!   30 continue
c     at radial barriers
      do 40 i=1,n1(klev+1)-1
      do 40 j=1,n2(klev+1)
        rc(i,j,0)=-rc(i,j,1)
        rc(i,j,n3(klev+1)+1)=-rc(i,j,n3(klev+1))
        psic(i,j,0)=-psic(i,j,1)
        psic(i,j,n3(klev+1)+1)=-psic(i,j,n3(klev+1))
   40 continue
!      do 40 i=0,n1(klev+1)
!      do 40 j=0,n2(klev+1)
!        rc(i,j,0)=0.d0
!        rc(i,j,n3(klev+1))=0.d0
!        psic(i,j,0)=0.d0
!        psic(i,j,n3(klev+1))=0.d0
!   40 continue
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
c==================================================
      subroutine finer(klev)
c
      INCLUDE "sect360mg.inb"
c minus real*8 psi(0:kr,0:kz,0:kp),rhs(0:kr,0:kz,0:kp)
c       real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
c       common /radius/r0,rmax   
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
      parameter(kr=33,kz=33,kp=161,mgen2=202590)
      real*8 a3d(0:kr,0:kz,0:kp),a1d(mgen2)
c
      do 10 k=0,n3
      do 10 j=0,n2
      do 10 i=0,n1
        ii=k*(n1+1)*(n2+1) + j*(n1+1) + i + i1
        if(ii.gt.202590)then
          print*,'FATAL ERROR IN MULTIGRID  ii=',ii
          stop
        endif
        a1d(ii)=a3d(i,j,k)
   10 continue
c!mpk      print*,ii
      return
      end
c==================================================
      subroutine trans2(a1d,a3d,i1,n1,n2,n3)
      implicit real*8(a-h,o-z)
      parameter(kr=33,kz=33,kp=161,mgen2=202590)
      real*8 a3d(0:kr,0:kz,0:kp),a1d(mgen2)
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
      parameter(kr=33,kz=33,kp=161,mgen2=202590)
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
      parameter(kx=33,ky=33,kz=161,mgen2=202590)
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
!mpk          print*,'n=',n
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
      parameter(kx=33,ky=33,kz=161,mgen2=202590)
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
      parameter(kx=33,ky=33,kz=161,mgen2=202590)
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
      parameter(kx=33,ky=33,kz=161,mgen2=202590)
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
      parameter(kx=33,ky=33,kz=161,mgen2=202590)
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
      parameter(kx=33,ky=33,kz=161,mgen2=202590)
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
      parameter(kx=33,ky=33,kz=161,mgen2=202590)
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
