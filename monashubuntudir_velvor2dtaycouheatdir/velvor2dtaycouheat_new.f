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
c        .edited 15-3-05 exclude derivatives in z                      c  
c        .9-5-05 work begins to implement new time and space           c
c         discretisation schemes in CVort only because 2D.             c 
c         Dufort-Frankel for diffusion term, Leapfrog for advection    c 
c         term, 4th order accurate compact space FD for advection term.c
c         Leapfrog has Robert filter. Starts first timestep with       c  
c         forward Euler                                                c    
c        .11-5-05 work begins to use compact FD for del_r T and        c
c         del_theta T terms in Cvort.                                  c
c        .20-5-05 changed these to improve convergence                 c
c         nrel1=1 nrel2=1; introducing Full Multigrid (Briggs et al);  c
c         removing line relaxation in phi direction.                   c 
c        .30-5-05 put in standard central 4th order for other terms.   c
c         (25 April 06 above is commented out)                         c
c         put in order 8 Shapiro filter for high wavenumbers.          c                  
c         if p=2pidx/L, the response is (1-sin^16 p/2)U_j, so          c
c         eliminates L=2dx (highest wave number) and damps 4dx         c
c         see Kalnay's book p. 102 (27 April 06 commented out for a    c
c         check)                                                       c
c        . A B An Bn Ap Bp welp are loaned to do something else now    c     
c        .weln vrms defined to stored new kinds of calc. concerning    c
c         b.l's                                                        c
c        .look for aug06 for removing unnecessary repeated calculationsc
c         in relax2 and relax3                                         c
c        .look for aug06 for for the removed trans1 calls in relax2,   c 
c         relax3. corser and finer                                     c
c        .14aug06 put in OpenMP directives to solve Poisson2 and       c
c         Poisson3 separately in two threads                           c
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=290,kz=65,kp=1610,mgen=100000)
      real*8 An(0:kr,1,0:kp),A(0:kr,1,0:kp),Ap(0:kr,1,0:kp)
      real*8 Bn(0:kr,1,0:kp),B(0:kr,1,0:kp),Bp(0:kr,1,0:kp)
      real*8 Cn(0:kr,1,0:kp),C(0:kr,1,0:kp),Cp(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 uelp(0:kr,1,0:kp),velp(0:kr,1,0:kp)
      real*8 welp(0:kr,1,0:kp),weln(0:kr,1,0:kp),vrms(0:kr,1,0:kp)
      real*8 urh(0:kr,1,0:kp),vrh(0:kr,1,0:kp),wrh(0:kr,1,0:kp)
      real*8 httrans(0:mgen), httrans2(0:mgen)
!mpk2d for output 
      real*4 r4out(1:82)
      real*8 ssumc(6,mgen),ssums(6,mgen)
      character*20 fname1,fname2,fname3,fname4,fname5
      character*20 fname6,fname7,fname8,fname9,ans
      common/grals/ssumc,ssums
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
!mpk2d      print*,'number of intervals in z-direction: '
!mpk2d      read*,nz
      nz=1
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
      print*,'Taylor number: '
      read*,Ta
!mpktc Gr is loaned to Ta
       Gr=Ta
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
!mpk2d      dz=zmax/dble(nz)
      dz=0.d0
      dp=2.d0*pi/dble(np)
c
c----------------------------------------------------------------------c
c     Set initial conditions                                           c
c----------------------------------------------------------------------c
c
      do 10 i=0,nr+1
!mpk2d      do 10 j=0,nz+1
      j=1
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
!mpktc    th(i,j,k)=1.d0-dlog(r)/dlog(r0)
        th(i,j,k)=1.d0-dlog(r)/dlog(r0)
c    +      +pert*dsin(2.d0*pi*(r/rmax-0.5d0))*dcos(2.d0*pi*z)
!mpk2d     +      +amp*(rmax-r)*(r-r0)*dcos(dble(nval)*phi)
c    +      +amp*dcos(pi*z/zmax)*(rmax-r)*(r-r0)*dcos(dble(nval)*phi)
!mpk2d         th(i,j,k)=th(i,j,k)+amp*(rmax-r)*(r-r0)*dcos(dble(nval)*phi)
!mpk2d     +             *th(i,j,k)

         if(i.eq.0)then 
!mpktc           th(i,j,k)=0.d0
           th(i,j,k)=0.d0
         elseif(i.eq.nr+1)then
!mpktc           th(i,j,k)=1.d0
           th(i,j,k)=1.d0
!         if(i.gt.12 .and. i.lt.101)then
         else
!mpk g77 doesn't have this           call random_number(randno)
!mpk but it doesn't matter i can do something like this:
!mpk for nz=1200, .../150 will give 8 azimuthal-waves
       randno=0.5d0+0.5d0*sin(dble(k)/150.d0*2.d0*pi)
           th(i,j,k)=th(i,j,k)+(r-rmax)/rmax*((randno*2.d0)-1.d0)
     &               *0.7d0*th(i,j,k)
         endif
  10  continue

!mpktc
!       do  k=0,np
!        Cn(nr,j,k)=(8.d0/3.d0*vel(nr+1,j,k)+vel(nr-1,j,k)/3.d0-
!     &      3.d0*vel(nr,j,k))/dr+vel(nr+1,j,k)/rmax
!       enddo
!mpktc

!mpk      do 11 j=1,nz-1
!mpk      th(nr/4,j,np/2)=th(nr/4,j,np/2)-0.2d0
!mpk  11  continue
!mpk2d         th(nr/4,1,np/2)=th(nr/4,1,np/2)-0.5d0*th(nr/4,1,np/2)

         ymoi = float(nr)/2.
         ymok = 3.*float(np)/4.
         ymwi = 5.d0
         ymri = 5.d0
         ymwk = 20.d0
         ymrk = 50.d0
!         call random_number(randno)
         ymamp = -0.d0
         do i = 0, nr+1 
            ymdi = abs(float(i) - ymoi)         
            ymii = max(0.d0, ymdi-ymri)/ymwi
            do k = 0, np
              ymdk = abs(float(k) - ymok)       
              ymkk = max(0.d0, ymdk-ymrk)/ymwk
              th(i,1,k) = th(i,1,k) + 
     &                    th(i,1,k)*ymamp*exp(-ymii**2-ymkk**2)
            enddo
         enddo

         ymoi = float(nr)/2.
         ymok = float(np)/4.
         ymwi = 5.d0
         ymri = 5.d0
         ymwk = 20.d0
         ymrk = 50.d0
!         call random_number(randno)
         ymamp = -0.d0
         do i = 0, nr+1 
            ymdi = abs(float(i) - ymoi)         
            ymii = max(0.d0, ymdi-ymri)/ymwi
            do k = 0, np
              ymdk = abs(float(k) - ymok)       
              ymkk = max(0.d0, ymdk-ymrk)/ymwk
              th(i,1,k) = th(i,1,k) +
     &                    th(i,1,k)*ymamp*exp(-ymii**2-ymkk**2)
            enddo
         enddo

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
!mpk2d        read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
        read(11)(((A(i,j,k),i=0,nr+1),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for vorticity B: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
!mpk2d         read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
        read(11)(((B(i,j,k),i=0,nr),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for vorticity C: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
!mpk2d         read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
        read(11)(((C(i,j,k),i=0,nr),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for velocity u: '
        read*,fname9
         open(11,file=fname9,form='unformatted')
!mpk2d         read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
        read(11)(((uel(i,j,k),i=0,nr),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for velocity v: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
!mpk2d         read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
        read(11)(((vel(i,j,k),i=0,nr+1),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for velocity w: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
!mpk2d         read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
        read(11)(((wel(i,j,k),i=0,nr+1),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for temperature: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
!mpk2d         read(11)nr,nz,np
!mpk2d         read(11)rmax,r0,zmax,zmin
        read(11)(((th(i,j,k),i=0,nr+1),j=1,1),k=0,np)
        close(11)
c
      endif
!mpk      call rufcon(th,nr+1,nz+1,np)
c
      do 50 i=0,nr+1
!mpk2d      do 50 j=0,nz+1
       j=1
      do 50 k=0,np
        Ap(i,j,k)=A(i,j,k)
        Bp(i,j,k)=B(i,j,k)
        Cp(i,j,k)=C(i,j,k)
        thp(i,j,k)=th(i,j,k)
        uelp(i,j,k)=uel(i,j,k)
        velp(i,j,k)=vel(i,j,k)
!mpk welp and weln are used for something else 
!mpk search for mpkzmtm
        welp(i,j,k)=0.d0
        weln(i,j,k)=0.d0
  50  continue
c
c----------------------------------------------------------------------c
c     Solve for vorticity, stream function, axial velocity and         c
c     temperature on internal points                                   c
c----------------------------------------------------------------------c
c
      ntot=0
c
cmpk don't want this now      open(12,file='thcvorttseriesfiner.dat',status='unknown')
      
       kt=1
       call Cvorticity1
        call temperature1
        call heat(trans1,trans2)
!mpk don't need this right now        call aaaaaaaa(ntot+kt)
!mpk2d        call rhsw(wrh)
!mpk2d        call poisson1(wel,welp,wrh,nr,nz,np,dr,dz,dp,nlev)
        call rhsv(vrh)
        call poisson2(vel,velp,vrh,nr,nz,np,dr,dz,dp,nlev)
        call rhsu(urh)
        call poisson3(uel,uelp,urh,nr,nz,np,dr,dz,dp,nlev)
        call update
        httrans(ntot+kt)=-r0*dlog(r0/rmax)*trans1/(2.d0*pi)
        httrans2(ntot+kt)=-rmax*dlog(r0/rmax)*trans2/(2.d0*pi)
      
      do 15 ku=1,9999

      do 20 kt=2,nts
!mpk2d        call rufcon(th,nr+1,nz+1,np)
!mpk2d        call Avorticity
!mpk2d        call Bvorticity
        call Cvorticity
        call temperature
        call heat(trans1,trans2)
!mpk don't need this right now        call aaaaaaaa(ntot+kt)
!mpk2d        call rhsw(wrh)
!mpk2d        call poisson1(wel,welp,wrh,nr,nz,np,dr,dz,dp,nlev)

!$OMP PARALLEL    

!$OMP sections
        call rhsv(vrh)
        call poisson2(vel,velp,vrh,nr,nz,np,dr,dz,dp,nlev)
!$OMP section
        call rhsu(urh)
        call poisson3(uel,uelp,urh,nr,nz,np,dr,dz,dp,nlev)
!$OMP end sections      
        
!$OMP END PARALLEL 

        call update
c
!mpk2d  httrans(ntot+kt)=-r0*dlog(r0/rmax)*trans1/(zmax*2.d0*pi)
        httrans(ntot+kt)=-r0*dlog(r0/rmax)*trans1/(2.d0*pi)
        httrans2(ntot+kt)=-rmax*dlog(r0/rmax)*trans2/(2.d0*pi)
        
!mpk2d save time        print*,httrans(ntot+kt),httrans2(ntot+kt),ntot+kt
c
       if(mod(kt,50).eq.0)then
        print*,httrans(ntot+kt),httrans2(ntot+kt),ntot+kt
        print*,'CFLMAXU=',cflmaxu
        print*,'CFLMAXV=',cflmaxv
!mpk2d        call rufcon(th,nr+1,nz+1,np)
       endif
cmpk don't want this now see above too       if(mod(kt,10).eq.0)then
c        write(12,201)th(15,1,2),th(64,1,2),th(112,1,2),
c     +               th(15,1,240),th(64,1,240),th(112,1,240),  
c     +               C(15,1,2),C(64,1,2),C(112,1,2),
c     +               C(15,1,240),C(64,1,240),C(112,1,240)
c       endif

!mpkzmtm calculating zonal-mean time-mean of the temperature variance
!mpk [\bar{ (T')^2 }] = [ \bar{T^2} - \bar{T}^2 ] see work note on
!mpk 7-6-05. borrow arrays weln (mean-sq T) and welp (mean T)
!mpk 16-6-05 Bp for mean-sq uel; vrms for mean-sq vel
       tkti=1.d0/dble(kt-1)
       do i=0,nr+1
        do k=0,np
         weln(i,1,k)=tkti*weln(i,1,k)*dble(kt-2)+tkti*th(i,1,k)**2.d0
         welp(i,1,k)=tkti*welp(i,1,k)*dble(kt-2)+tkti*th(i,1,k)
         Bp(i,1,k)=tkti*Bp(i,1,k)*dble(kt-2)+tkti*uel(i,1,k)**2.d0
         vrms(i,1,k)=tkti*vrms(i,1,k)*dble(kt-2)+tkti*vel(i,1,k)**2.d0
        enddo
       enddo
!mpk now see after the time integration loop
       
       
!mpk checking CFL conditions-------------------------------------------
       cflmaxu=0.d0
       do i=1,nr
        do k=0,np
          cfl=dabs(uel(i,1,k)*dt/dr)
          if(cfl.gt.cflmaxu)then 
             cflmaxu=cfl
             ucfl=dabs(uel(i,1,k))
          endif
        enddo
       enddo
       if(cflmaxu.gt.1.)then
        print*,'WARNING! CFL EXCEEDED!'
       endif
c       print*,'CFLMAXU=',cflmaxu

       cflmaxv=0.d0
       do i=1,nr-1
        if(i.eq.0)then
          r=r0
        elseif(i.eq.nr+1)then
          r=rmax
        else
          r=r0+(dble(i)-0.5d0)*dr
        endif
        do k=0,np
          cfl=dabs(vel(i,1,k)*dt/(r*dp))
          if(cfl.gt.cflmaxv)then 
            cflmaxv=cfl
            vcfl=dabs(vel(i,1,k))
            rcfl=r
          endif
        enddo
       enddo
       if(cflmaxv.gt.1.)then
        print*,'WARNING! CFL EXCEEDED!'
       endif
c       print*,'CFLMAXV=',cflmaxv
       
!mpk       if (kt.ge.10)then
!mpk        if(cflmaxu.gt.cflmaxv)then
!mpk         dt=0.1d0*dr/ucfl
!mpk         print*,'dt=',dt
!mpk        else
!mpk         dt=0.1d0*dp*rcfl/vcfl
!mpk         print*,'dt=',dt
!mpk        endif
!mpk       endif
   

  20  continue
c
      ntot=ntot+nts
c-------------------------------------------writing to files begins
       open(14,file='zmtmtempvar.dat',status='unknown')
      do i=0,nr+1
       zmtmsqtemp=0.d0
       zmsqtmtemp=0.d0
       zmtmsquel=0.d0
       zmtmsqwel=0.d0
       do k=0,np
         zmtmsqtemp=zmtmsqtemp+weln(i,1,k)
         zmsqtmtemp=zmsqtmtemp+welp(i,1,k)**2.d0
         zmtmsquel=zmtmsquel+Bp(i,1,k)
         zmtmsqwel=zmtmsqwel+vrms(i,1,k)
       enddo
         whatever=(zmtmsqtemp-zmsqtmtemp)/dble(np+1)
         whatever2=(zmtmsquel/dble(np+1))**0.5d0
         whatever3=(zmtmsqwel/dble(np+1))**0.5d0
         write(14,202)whatever,whatever2,whatever3
      enddo
      close(14)    
c
c---------------------------------------------------------------------c
c     Write vorticity, stream function, velocity and temperature data c
c     to file                                                         c
c---------------------------------------------------------------------c
c
      open(11,file=fname1,form='unformatted')
!mpk2d      write(11)nr,nz,np
!mpk2d      write(11)rmax,r0,zmax,zmin
      write(11)(((A(i,j,k),i=0,nr+1),j=1,1),k=0,np)
      close(11)
c
      open(11,file=fname2,form='unformatted')
!mpk2d      write(11)nr,nz,np
!mpk2d      write(11)rmax,r0,zmax,zmin
      write(11)(((B(i,j,k),i=0,nr),j=1,1),k=0,np)
      close(11)
c
      open(11,file=fname3,form='unformatted')
!mpk2d      write(11)nr,nz,np
!mpk2d      write(11)rmax,r0,zmax,zmin
      write(11)(((C(i,j,k),i=0,nr),j=1,1),k=0,np)
      close(11)
c
      open(11,file=fname4,form='unformatted')
!mpk2d      write(11)nr,nz,np
!mpk2d      write(11)rmax,r0,zmax,zmin
      write(11)(((uel(i,j,k),i=0,nr),j=1,1),k=0,np)
      close(11)
c
      open(11,file=fname5,form='unformatted')
!mpk2d      write(11)nr,nz,np
!mpk2d      write(11)rmax,r0,zmax,zmin
      write(11)(((vel(i,j,k),i=0,nr+1),j=1,1),k=0,np)
      close(11)
c
      open(11,file=fname6,form='unformatted')
!mpk2d      write(11)nr,nz,np
!mpk2d      write(11)rmax,r0,zmax,zmin
      write(11)(((wel(i,j,k),i=0,nr+1),j=1,1),k=0,np)
      close(11)
c
!mpk2d       open(11,file=fname7,form='unformatted',access='sequential'
!mpk2d      +     ,status='new')
            open(11,file=fname7,form='unformatted')
!mpk2d      write(11)nr,nz,np
!mpk2d      write(11)rmax,r0,zmax,zmin
!mpk2d      write(11,201)(((th(i,j,k),i=0,nr+1),j=1,1),k=0,np)
            write(11)(((th(i,j,k),i=0,nr+1),j=1,1),k=0,np)
!mpk2d        do ken=0,np
!mpk2d          do ien=0,nr+1
!mpk2d            r4out(ien+1)=th(ien,1,ken)
!mpk2d          enddo
!mpk2d           write(11)r4out 
!mpk2d        enddo
      close(11)
c
      if(fname8 .ne. 'no')then
        open(11,file=fname8,status='unknown')
        write(11,200)((i-1)*dt,httrans(i),ssumc(1,i),ssums(1,i),
     +  ssumc(2,i),ssums(2,i),ssumc(3,i),ssums(3,i),ssumc(4,i),
     +  ssums(4,i),ssumc(5,i),ssums(5,i),ssumc(6,i),ssums(6,i),i=1,ntot)
        close(11)
      endif

c
  200  format(14(1x,e16.8))
  201  format(12(1x,e16.8))
  202  format(3(1x,e16.8))
c---------------------------------------------writing to files ends
      
      print*,'more iterations? (y/n): '
      read*,ans
      if(ans.eq.'n'.or.ans.eq.'no')goto 30
      print*,'no. of iterations: '
      read*,nts
c
  15  continue
c
  30  continue

      close(12)
      stop
      end
cc----------------------------------------------------------------------c
c     Subroutine CVORTICITY                                            c
c     =====================                                            c
c     Calculates values of vorticity using the duFort-Frankel method   c
c----------------------------------------------------------------------c
c
      subroutine Cvorticity1
c   
      implicit real*8(a-h,o-z)
      parameter(kr=290,kz=33,kp=1610,mgen=100000)
      real*8 An(0:kr,1,0:kp),A(0:kr,1,0:kp),Ap(0:kr,1,0:kp)
      real*8 Bn(0:kr,1,0:kp),B(0:kr,1,0:kp),Bp(0:kr,1,0:kp)
      real*8 Cn(0:kr,1,0:kp),C(0:kr,1,0:kp),Cp(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 uelp(0:kr,1,0:kp),velp(0:kr,1,0:kp)
      real*8 welp(0:kr,1,0:kp)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      dri = 1.d0/dr
      dpi = 1.d0/dp
      j = 1

!mpk2d compact FD, finding Cr    
      do k=0,np
       do i=1,nr-1
         tdmaa(i)=1.d0
         tdmab(i)=4.d0
         tdmac(i)=1.d0
         if(i.eq.1)then
          tdmar(i)=3.d0*(C(2,j,k)-C(0,j,k))*dri-(-1.5d0*C(0,j,k)+2.d0*
     +             C(1,j,k)-0.5d0*C(2,j,k))*dri
         elseif(i.eq.nr-1)then
          tdmar(i)=3.d0*(C(nr,j,k)-C(nr-2,j,k))*dri-(0.5d0*C(nr-2,j,k)-
     +             2.d0*C(nr-1,j,k)+1.5d0*C(nr,j,k))*dri
         else
          tdmar(i)=3.d0*(C(i+1,j,k)-C(i-1,j,k))*dri
         endif   
       enddo    
         call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,nr-1)
         do ii=1,nr-1
!mpk2d borrowing A to store Cr
          A(ii,j,k)=tdmas(ii)
         enddo
      enddo

!mpk2d compact FD, finding Ch
      do i=1,nr-1
       do k=0,np
         tdmaa(k+1)=1.d0
         tdmab(k+1)=4.d0
         tdmac(k+1)=1.d0
         if(k.eq.np)then
           tdmar(np+1)=3.d0*(C(i,j,1)-C(i,j,np-1))*dpi
         elseif(k.eq.0)then
           tdmar(1)=3.d0*(C(i,j,1)-C(i,j,np-1))*dpi
         else
           tdmar(k+1)=3.d0*(C(i,j,k+1)-C(i,j,k-1))*dpi
         endif 
       enddo
       call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,np+1)
       do kk=0,np
!mpk2d borrowing B to store Ch
          B(i,j,kk)=tdmas(kk+1)
       enddo
      enddo
!mpk2d------------------compact FD. finding thr
      do k=0,np
       do i=1,nr-1
         tdmaa(i)=1.d0
         tdmab(i)=4.d0
         tdmac(i)=1.d0
         ka=k+1
         if(k.eq.np)ka=1
         if(i.eq.1)then
          tdmar(i)=1.5d0*(th(2,j,k)+th(2,j,ka)-th(1,j,k)-th(1,j,ka))
     +      *dri-(-2./3.*(th(0,j,k)+th(0,j,ka)) +0.5d0*(th(1,j,k)+
     +      th(1,j,ka)) + 1./6.*(th(2,j,k)+th(2,j,ka)))*dri
         elseif(i.eq.nr-1)then
          tdmar(i)=1.5d0*(th(nr,j,k)+th(nr,j,ka)-th(nr-1,j,k)-
     +      th(nr-1,j,ka))*dri-(-1./6.*(th(nr-1,j,k)+th(nr-1,j,ka))-
     +      0.5d0*(th(nr,j,k)+th(nr,j,ka))+2./3.*(th(nr+1,j,k)+
     +      th(nr+1,j,ka) ) )*dri
         else
          tdmar(i)=1.5d0*(th(i+1,j,k)+th(i+1,j,ka)-th(i,j,k)
     +             -th(i,j,ka))*dri
         endif   
       enddo    
         call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,nr-1)
         do ii=1,nr-1
!mpk2d borrowing An to store thr
          An(ii,j,k)=tdmas(ii)
         enddo
      enddo

!mpk2d------------------ compact FD, finding thh
      do i=1,nr-1
       do k=0,np
         tdmaa(k+1)=1.d0
         tdmab(k+1)=4.d0
         tdmac(k+1)=1.d0
         ka=k+1
         if(k.eq.np)ka=1
           tdmar(k+1)=1.5d0*(th(i+1,j,ka)+th(i,j,ka)-th(i+1,j,k)-
     +                th(i,j,k))*dpi
       enddo
       call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,np+1)
       do kk=0,np
!mpk2d borrowing B to store Ch
          Bn(i,j,kk)=tdmas(kk+1)
       enddo
      enddo

!mpk2d-----------------------------------------------------------------
      do 10 i=1,nr-1
!mpk2d      do 10 j=1,nz
      do 10 k=0,np
c
        ip=i+1
        im=i-1
!mpk2d        jp=j+1
!mpk2d        jm=j-1
        km=k-1
        kmm=k-2
        ka=k+1
        kaa=k+2
        if(k.eq.1)kmm=np-1
        if(k.eq.np-1)kaa=1
        if(k.eq.0)then 
         km=np-1
         kmm=km-1
        endif
        if(k.eq.np)then
          ka=1
          kaa=ka+1
        endif
c
        r=r0+dble(i)*dr
        ri=1.d0/r
        u=(uel(i,j,ka)+uel(i,j,k))*0.5d0
        v=(vel(i,j,k)+vel(ip,j,k))*0.5d0       
        tem=(th(i,j,k)+th(ip,j,k)+th(ip,j,ka)+th(i,j,ka))*0.25d0

!mpk2d Compact FD        thr=(th(ip,j,k)+th(ip,j,ka)-th(i,j,k)-th(i,j,ka))*0.5d0*dri
        thr=An(i,j,k)
c
!mpk2d now the 4th order central difference  
         Chha=(C(i,j,ka)+C(i,j,km))*dpi**2
!22         Chha=(-0.0833333333d0*C(i,j,kmm)+1.3333333333d0*C(i,j,km)+
!22     +        1.3333333333d0*C(i,j,ka)-0.0833333333d0*C(i,j,kaa))*
!22     +        dpi**2
!mpk2d Compact FD        Ch=(C(i,j,ka)-C(i,j,km))*0.5d0*dpi
        Ch=B(i,j,k)
c
!22         if(i.eq.1 .or. i.eq.nr-1)then
          Crra=(C(i-1,j,k)+C(i+1,j,k))*dri**2
          cnr=-1.d0*dri**2
          Cra=(C(i+1,j,k)-C(i-1,j,k))*0.5d0*dri
          Cr=A(i,j,k)
!22         else
!mpk2d 4th order CF away from boundaries
!22          Crra=(-0.0833333333d0*C(i-2,j,k)+1.3333333333d0*C(i-1,j,k)+
!22     +        1.3333333333d0*C(i+1,j,k)-0.0833333333d0*C(i+2,j,k))*
!22     +        dri**2
!22          cnr=-1.25d0*dri**2
!22          Cra=(0.0833333333d0*C(i-2,j,k)-0.6666666666d0*C(i-1,j,k)+
!22     +       0.6666666666d0*C(i+1,j,k)-0.0833333333d0*C(i+2,j,k))*
!22     +        dri
!22          Cr=A(i,j,k)
!22         endif
!mpk2d compact FD        Cra=(C(i+1,j,k)-C(i-1,j,k))*0.5d0*dri
!mpk2d Compact FD       Cr=Cra
c
        welh=0.d0
!mpk2d Compact FD    thh=(th(i,j,ka)+th(ip,j,ka)-th(i,j,k)-th(ip,j,k))*0.5d0*dpi
        thh=Bn(i,j,k)
c
!mpk2d        cnr=-1.25d0*dri**2
!22         cnh=(-1.25d0*dpi**2)*(ri**2)
        cnh=(-1.d0*dpi**2)*(ri**2)
c forward timestepping just for first timestep
        RHS1=Crra+Cra*ri+Chha*(ri**2)+C(i,j,k)/(Ptl*dt)
        ULT=(u*Cr+v*Ch*ri)/Ptl
        RHS2=-Ra*thh+2.d0*Ra*(-thr*u-thh*v*ri)/Re/Ptl
        cf1=1.d0/(Ptl*dt)-cnr-cnh
        cf2=cnr+cnh
        Cn(i,j,k)=(cf2*Cp(i,j,k)-ULT-RHS2+RHS1)/cf1
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
      parameter(kr=290,kz=33,kp=1610,mgen=100000)
      real*8 An(0:kr,1,0:kp),A(0:kr,1,0:kp),Ap(0:kr,1,0:kp)
      real*8 Bn(0:kr,1,0:kp),B(0:kr,1,0:kp),Bp(0:kr,1,0:kp)
      real*8 Cn(0:kr,1,0:kp),C(0:kr,1,0:kp),Cp(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 uelp(0:kr,1,0:kp),velp(0:kr,1,0:kp)
      real*8 welp(0:kr,1,0:kp)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      dri = 1.d0/dr
      dpi = 1.d0/dp
      j = 1

!mpk2d------------------- compact FD, finding Cr    
      do k=0,np
       do i=1,nr-1
         tdmaa(i)=1.d0
         tdmab(i)=4.d0
         tdmac(i)=1.d0
         if(i.eq.1)then
          tdmar(i)=3.d0*(C(2,j,k)-C(0,j,k))*dri-(-1.5d0*C(0,j,k)+2.d0*
     +             C(1,j,k)-0.5d0*C(2,j,k))*dri
         elseif(i.eq.nr-1)then
          tdmar(i)=3.d0*(C(nr,j,k)-C(nr-2,j,k))*dri-(0.5d0*C(nr-2,j,k)-
     +             2.d0*C(nr-1,j,k)+1.5d0*C(nr,j,k))*dri
         else
          tdmar(i)=3.d0*(C(i+1,j,k)-C(i-1,j,k))*dri
         endif   
       enddo    
         call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,nr-1)
         do ii=1,nr-1
!mpk2d borrowing A to store Cr
          A(ii,j,k)=tdmas(ii)
         enddo
      enddo

!mpk2d------------------ compact FD, finding Ch
      do i=1,nr-1
       do k=0,np
         tdmaa(k+1)=1.d0
         tdmab(k+1)=4.d0
         tdmac(k+1)=1.d0
         if(k.eq.np)then
           tdmar(np+1)=3.d0*(C(i,j,1)-C(i,j,np-1))*dpi
         elseif(k.eq.0)then
           tdmar(1)=3.d0*(C(i,j,1)-C(i,j,np-1))*dpi
         else
           tdmar(k+1)=3.d0*(C(i,j,k+1)-C(i,j,k-1))*dpi
         endif 
       enddo
       call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,np+1)
       do kk=0,np
!mpk2d borrowing B to store Ch
          B(i,j,kk)=tdmas(kk+1)
       enddo
      enddo

!mpk2d------------------compact FD. finding thr
      do k=0,np
       do i=1,nr-1
         tdmaa(i)=1.d0
         tdmab(i)=4.d0
         tdmac(i)=1.d0
         ka=k+1
         if(k.eq.np)ka=1
         if(i.eq.1)then
          tdmar(i)=1.5d0*(th(2,j,k)+th(2,j,ka)-th(1,j,k)-th(1,j,ka))
     +      *dri-(-2./3.*(th(0,j,k)+th(0,j,ka)) +0.5d0*(th(1,j,k)+
     +      th(1,j,ka)) + 1./6.*(th(2,j,k)+th(2,j,ka)))*dri
         elseif(i.eq.nr-1)then
          tdmar(i)=1.5d0*(th(nr,j,k)+th(nr,j,ka)-th(nr-1,j,k)-
     +      th(nr-1,j,ka))*dri-(-1./6.*(th(nr-1,j,k)+th(nr-1,j,ka))-
     +      0.5d0*(th(nr,j,k)+th(nr,j,ka))+2./3.*(th(nr+1,j,k)+
     +      th(nr+1,j,ka) ) )*dri
         else
          tdmar(i)=1.5d0*(th(i+1,j,k)+th(i+1,j,ka)-th(i,j,k)
     +             -th(i,j,ka))*dri
         endif   
       enddo    
         call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,nr-1)
         do ii=1,nr-1
!mpk2d borrowing An to store thr
          An(ii,j,k)=tdmas(ii)
         enddo
      enddo
!mpk2d------------------ compact FD, finding thh
      do i=1,nr-1
       do k=0,np
         tdmaa(k+1)=1.d0
         tdmab(k+1)=4.d0
         tdmac(k+1)=1.d0
         ka=k+1
         if(k.eq.np)ka=1
           tdmar(k+1)=1.5d0*(th(i+1,j,ka)+th(i,j,ka)-th(i+1,j,k)-
     +                th(i,j,k))*dpi
       enddo
       call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,np+1)
       do kk=0,np
!mpk2d borrowing B to store Ch
          Bn(i,j,kk)=tdmas(kk+1)
       enddo
      enddo

!mpk2d----------------------------------------timestep begins 
      do 10 i=1,nr-1
!mpk2d      do 10 j=1,nz
      do 10 k=0,np
c
        ip=i+1
        im=i-1
!mpk2d        jp=j+1
!mpk2d        jm=j-1
        km=k-1
        kmm=km-1
        ka=k+1
        kaa=k+2
        if(k.eq.1)kmm=np-1
        if(k.eq.np-1)kaa=1
        if(k.eq.0)then 
         km=np-1
         kmm=km-1
        endif
        if(k.eq.np)then
          ka=1
          kaa=ka+1
        endif
c
        r=r0+dble(i)*dr
        ri=1.d0/r
        u=(uel(i,j,ka)+uel(i,j,k))*0.5d0
        v=(vel(i,j,k)+vel(ip,j,k))*0.5d0       
        tem=(th(i,j,k)+th(ip,j,k)+th(ip,j,ka)+th(i,j,ka))*0.25d0

!mpk2d Compact FD        thr=(th(ip,j,k)+th(ip,j,ka)-th(i,j,k)-th(i,j,ka))*0.5d0*dri
        thr=An(i,j,k)
c
!mpk2d now the 4th order central difference  
         Chha=(C(i,j,ka)+C(i,j,km))*dpi**2
!22         Chha=(-0.0833333333d0*C(i,j,kmm)+1.3333333333d0*C(i,j,km)+
!22     +        1.3333333333d0*C(i,j,ka)-0.0833333333d0*C(i,j,kaa))*
!22     +        dpi**2
!mpk2d Compact FD        Ch=(C(i,j,ka)-C(i,j,km))*0.5d0*dpi
        Ch=B(i,j,k)
c
!22        if(i.eq.1 .or. i.eq.nr-1)then
         Crra=(C(i-1,j,k)+C(i+1,j,k))*dri**2
         cnr=-1.d0*dri**2
         Cra=(C(i+1,j,k)-C(i-1,j,k))*0.5d0*dri
         Cr=A(i,j,k)
!22        else
!mpk2d 4th order CF away from boundaries
!22         Crra=(-0.0833333333d0*C(i-2,j,k)+1.3333333333d0*C(i-1,j,k)+
!22    +        1.3333333333d0*C(i+1,j,k)-0.0833333333d0*C(i+2,j,k))*
!22     +        dri**2
!22          cnr=-1.25d0*dri**2
!22          Cra=(0.0833333333d0*C(i-2,j,k)-0.6666666666d0*C(i-1,j,k)+
!22     +       0.6666666666d0*C(i+1,j,k)-0.0833333333d0*C(i+2,j,k))*
!22     +        dri
!22          Cr=A(i,j,k)
!22        endif

!mpk2d Compact FD          Cra=(C(i+1,j,k)-C(i-1,j,k))*0.5d0*dri
!mpk2d Compact FD       Cr=Cra
c
        welh=0.d0
!mpk2d Compact FD        thh=(th(i,j,ka)+th(ip,j,ka)-th(i,j,k)-th(ip,j,k))*0.5d0*dpi
        thh=Bn(i,j,k)
c
!mpk2d doing 4th order for diffusion term so adding -2.5/2 to central term
!mpk2d         cnr=-1.25d0*dri**2
!22         cnh=(-1.25d0*dpi**2)*(ri**2)
        cnh=(-1.d0*dpi**2)*(ri**2)
c back to truly leapfrog/Dufort-Frankel
        RHS1=Crra+Cra*ri+Chha*(ri**2)+0.5d0*Cp(i,j,k)/(Ptl*dt)
        ULT=(u*Cr+v*Ch*ri)/Ptl
        RHS2=-Ra*thh+2.d0*Ra*(-thr*u-thh*v*ri)/Re/Ptl
        cf1=0.5d0/(Ptl*dt)-cnr-cnh
        cf2=cnr+cnh
        Cn(i,j,k)=(cf2*Cp(i,j,k)-ULT-RHS2+RHS1)/cf1
        
!mpk also stored in loaned Ap for the Shapiro filter below
        Ap(i,j,k)=Cn(i,j,k)
c
  10  continue

!mpk Shapiro filter of order 2n=16, see Kalnay p. 102
!mpk \bar U^{2n}_{j} = [1-(-D^n)]U_j, where diffusion operator 
!mpk DU_j=0.25(U_{j+1}-2U_{j}+U_{j-1}) is applied to the original 
!mpk field n times. 2n = 16 is good. reduce norder for high resolution 


c        do norder=1,6

!mpk first applied in r direction and stored in An
c         do k=0,np
c          do i=1,nr-1
c            An(i,j,k)=-0.25d0*(Ap(i+1,j,k)-2.d0*Ap(i,j,k)+Ap(i-1,j,k))
c          enddo
c         enddo

c         do k=0,np
c          do i=1,nr-1
c           Ap(i,j,k)=An(i,j,k)
c          enddo
c         enddo

c        enddo

c        do i=1,nr-1
c         do k=0,np
c          Ap(i,j,k)=Cn(i,j,k)-Ap(i,j,k)
c         enddo
c        enddo
        
!mpk then applied in tangential direction and stored back to Ap  
c        do norder=1,8
      
c         do i=1,nr-1
c          do k=0,np
c           km=k-1
c           ka=k+1
c           if(k.eq.0)km=np-1
c           if(k.eq.np)ka=1          
c            An(i,j,k)=-0.25d0*(Ap(i,j,ka)-2.d0*Ap(i,j,k)+Ap(i,j,km))
c          enddo
c         enddo

c         do k=0,np
c          do i=1,nr-1
c           Ap(i,j,k)=An(i,j,k)
c           enddo
c         enddo

c        enddo
      
c        do i=1,nr-1
c         do k=0,np
c          Cn(i,j,k)=Cn(i,j,k)-Ap(i,j,k)
c         enddo
c        enddo

!mpk Robert filter
        do i=1,nr-1
         do k=0,np
          C(i,j,k)=C(i,j,k)+0.02d0*(Cn(i,j,k)-2.d0*C(i,j,k)+Cp(i,j,k))
         enddo
        enddo
   
      return
      end
c
c----------------------------------------------------------------------c
c     subroutine TEMPERATURE                                           c
c     ======================                                           c
c     Calculates temperature using the duFort-Frankel mehod            c
c----------------------------------------------------------------------c
c
      subroutine temperature1
c
      implicit real*8(a-h,o-z)
      parameter(kr=290,kz=33,kp=1610,mgen=100000)
      real*8 An(0:kr,1,0:kp),A(0:kr,1,0:kp),Ap(0:kr,1,0:kp)
      real*8 Bn(0:kr,1,0:kp),B(0:kr,1,0:kp),Bp(0:kr,1,0:kp)
      real*8 Cn(0:kr,1,0:kp),C(0:kr,1,0:kp),Cp(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 uelp(0:kr,1,0:kp),velp(0:kr,1,0:kp)
      real*8 welp(0:kr,1,0:kp)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      dri = 1.d0/dr
      dpi = 1.d0/dp
      j = 1

!mpk2d compact FD, finding thr    
      do k=0,np
       do i=1,nr
         tdmaa(i)=1.d0
         tdmab(i)=4.d0
         tdmac(i)=1.d0
         if(i.eq.1)then
          tdmab(i)=1.d0
          tdmac(i)=2.d0
          tdmar(i)=dri*(0.533333333d0*th(0,j,k)-3.6d0*th(1,j,k)+
     +             2.966666667d0*th(2,j,k)+0.1d0*th(4,j,k))
         elseif(i.eq.nr)then
          tdmaa(i)=2.d0
          tdmab(i)=1.d0
          tdmar(i)=dri*(-0.1d0*th(nr-3,j,k)-2.966666667d0*th(nr-1,j,k)+
     +             3.6d0*th(nr,j,k)-0.533333333d0*th(nr+1,j,k))
         else
          tdmar(i)=3.d0*(th(i+1,j,k)-th(i-1,j,k))*dri
         endif   
       enddo    
         call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,nr)
         do ii=1,nr
!mpk2d borrowing A to store thr
          A(ii,j,k)=tdmas(ii)
         enddo
      enddo

!mpk2d compact FD, finding thh
      do i=1,nr
       do k=0,np
         tdmaa(k+1)=1.d0
         tdmab(k+1)=4.d0
         tdmac(k+1)=1.d0
         if(k.eq.np)then
           tdmar(np+1)=3.d0*(th(i,j,1)-th(i,j,np-1))*dpi
         elseif(k.eq.0)then
           tdmar(1)=3.d0*(th(i,j,1)-th(i,j,np-1))*dpi
         else
           tdmar(k+1)=3.d0*(th(i,j,k+1)-th(i,j,k-1))*dpi
         endif 
       enddo
       call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,np+1)
       do kk=0,np
!mpk2d borrowing B to store thh
          B(i,j,kk)=tdmas(kk+1)
       enddo
      enddo

!mpk2d-----------------------------timestepping begins
      do 10 i=1,nr
!mpk2d      do 10 j=0,nz+1
      do 10 k=0,np
c
        ip=i+1
        im=i-1
!mpk2d        jp=j+1
!mpk2d        jm=j-1
        km=k-1
        kmm=km-1
        ka=k+1
        kaa=ka+1
        if(k.eq.1)kmm=np-1
        if(k.eq.np-1)kaa=1
        if(k.eq.0)then 
         km=np-1
         kmm=km-1
        endif
        if(k.eq.np)then
          ka=1
          kaa=ka+1
        endif
c
        r=r0+(dble(i)-0.5d0)*dr
        ri=1.d0/r
          u=(uel(i,j,k)+uel(im,j,k))*0.5d0
          v=(vel(i,j,k)+vel(i,j,km))*0.5d0
c
        if(i.eq.1)then
         thrra=(3.2d0*th(i-1,j,k)+2.d0*th(i+1,j,k)-0.2d0*th(i+2,j,k))*
     +         dri**2
         thra=(-4.d0*th(i-1,j,k)/3.d0+th(i+1,j,k)/3.d0)*dri
!mpk2d         thr=(-4.d0*th(i-1,j,k)/3.d0+th(i,j,k)+th(i+1,j,k)/3.d0)*dri
          thr=A(i,j,k)
!mpk 11aug07 is this bug serious?         cnr=-2.5d0*dri**2+0.5d0*dr*ri
!mpk 11aug07 corrected here
           cnr=-2.5d0*dri**2+0.5d0*dri*ri
        elseif(i.eq.nr)then
          thrra=(3.2d0*th(i+1,j,k)+2.d0*th(i-1,j,k)-0.2d0*th(i-2,j,k))*
     +          dri**2
          thra=(-th(i-1,j,k)/3.d0+4.d0*th(i+1,j,k)/3.d0)*dri
!mpk2d          thr=(-th(i-1,j,k)/3.d0-th(i,j,k)+4.d0*th(i+1,j,k)/3.d0)*dri
          thr=A(i,j,k)
          cnr=-2.5d0*dri**2-0.5d0*dri*ri
!22        elseif(i.eq.2 .or. i.eq.nr-1)then
        else
          thrra=(th(i-1,j,k)+th(i+1,j,k))*dri**2
          thra=(th(i+1,j,k)-th(i-1,j,k))*0.5d0*dri
!mpk2d          thr=thra
          thr=A(i,j,k)
c
          cnr=-1.d0*dri**2
!22         else
!mpk2d 4th order central diffrence
!22            thrra=(-0.083333333d0*th(i-2,j,k)+1.333333333d0*th(i-1,j,k)+
!22     +        1.333333333d0*th(i+1,j,k)-0.083333333d0*th(i+2,j,k))*
!22     +        dri**2
!22           thra=(0.0833333333d0*th(i-2,j,k)-0.6666666666d0*th(i-1,j,k)+
!22     +       0.6666666666d0*th(i+1,j,k)-0.0833333333d0*th(i+2,j,k))*
!22     +         dri
!22            thr=A(i,j,k)
!22            cnr=-1.25*dri**2
        endif
c
!mpk2d 4th order central diffrence   
         thhha=(th(i,j,ka)+th(i,j,km))*dpi**2
!22         thhha=(-0.083333333d0*th(i,j,kmm)+1.333333333d0*th(i,j,km)+
!22     +        1.333333333d0*th(i,j,ka)-0.083333333d0*th(i,j,kaa))*
!22     +        dpi**2
!mpk2d Compact FD       thh=(th(i,j,ka)-th(i,j,km))*0.5d0*dpi
        thh=B(i,j,k)
c
!22        cnp=-(1.25d0*dpi**2)*ri**2
        cnp=-(1.d0*dpi**2)*ri**2
c just for the first timestep, forward timestepping
        RHS1=thrra+thra*ri+thhha*ri**2+th(i,j,k)/dt  
        ULT=u*thr+v*thh*ri
        cf1=1.d0/dt-cnr-cnp
        cf2=cnr+cnp
        thn(i,j,k)=(cf2*thp(i,j,k)-ULT+RHS1)/cf1
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
      parameter(kr=290,kz=33,kp=1610,mgen=100000)
      real*8 An(0:kr,1,0:kp),A(0:kr,1,0:kp),Ap(0:kr,1,0:kp)
      real*8 Bn(0:kr,1,0:kp),B(0:kr,1,0:kp),Bp(0:kr,1,0:kp)
      real*8 Cn(0:kr,1,0:kp),C(0:kr,1,0:kp),Cp(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 uelp(0:kr,1,0:kp),velp(0:kr,1,0:kp)
      real*8 welp(0:kr,1,0:kp)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      common/val/Ra,Re,Ptl,Gr
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      dri = 1.d0/dr
      dpi = 1.d0/dp
      j = 1

!mpk2d compact FD, finding thr    
      do k=0,np
       do i=1,nr
         tdmaa(i)=1.d0
         tdmab(i)=4.d0
         tdmac(i)=1.d0
         if(i.eq.1)then
          tdmab(i)=1.d0
          tdmac(i)=2.d0
          tdmar(i)=dri*(0.533333333d0*th(0,j,k)-3.6d0*th(1,j,k)+
     +             2.966666667d0*th(2,j,k)+0.1d0*th(4,j,k))
         elseif(i.eq.nr)then
          tdmaa(i)=2.d0
          tdmab(i)=1.d0
          tdmar(i)=dri*(-0.1d0*th(nr-3,j,k)-2.966666667d0*th(nr-1,j,k)+
     +             3.6d0*th(nr,j,k)-0.533333333d0*th(nr+1,j,k))
         else
          tdmar(i)=3.d0*(th(i+1,j,k)-th(i-1,j,k))*dri
         endif   
       enddo    
         call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,nr)
         do ii=1,nr
!mpk2d borrowing A to store thr
          A(ii,j,k)=tdmas(ii)
         enddo
      enddo

!mpk2d compact FD, finding thh
      do i=1,nr
       do k=0,np
         tdmaa(k+1)=1.d0
         tdmab(k+1)=4.d0
         tdmac(k+1)=1.d0
         if(k.eq.np)then
           tdmar(np+1)=3.d0*(th(i,j,1)-th(i,j,np-1))*dpi
         elseif(k.eq.0)then
           tdmar(1)=3.d0*(th(i,j,1)-th(i,j,np-1))*dpi
         else
           tdmar(k+1)=3.d0*(th(i,j,k+1)-th(i,j,k-1))*dpi
         endif 
       enddo
       call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,np+1)
       do kk=0,np
!mpk2d borrowing B to store thh
          B(i,j,kk)=tdmas(kk+1)
       enddo
      enddo

!mpk2d------------------------------------------timestep begins
      do 10 i=1,nr
!mpk2d      do 10 j=0,nz+1
      do 10 k=0,np
c
        ip=i+1
        im=i-1
!mpk2d        jp=j+1
!mpk2d        jm=j-1
        km=k-1
        kmm=km-1
        ka=k+1
        kaa=ka+1
        if(k.eq.1)kmm=np-1
        if(k.eq.np-1)kaa=1
        if(k.eq.0)then 
         km=np-1
         kmm=km-1
        endif
        if(k.eq.np)then
          ka=1
          kaa=ka+1
        endif
c
        r=r0+(dble(i)-0.5d0)*dr
        ri=1.d0/r
          u=(uel(i,j,k)+uel(im,j,k))*0.5d0
          v=(vel(i,j,k)+vel(i,j,km))*0.5d0
c
        if(i.eq.1)then
         thrra=(3.2d0*th(i-1,j,k)+2.d0*th(i+1,j,k)-0.2d0*th(i+2,j,k))*
     +         dri**2
         thra=(-4.d0*th(i-1,j,k)/3.d0+th(i+1,j,k)/3.d0)*dri
!mpk2d         thr=(-4.d0*th(i-1,j,k)/3.d0+th(i,j,k)+th(i+1,j,k)/3.d0)*dri
         thr=A(i,j,k)
!mpk 11aug07 is this bug serious? cnr=-2.5d0*dri**2+0.5d0*dr*ri
!mpk 11aug07 corrected here
         cnr=-2.5d0*dri**2+0.5d0*dri*ri
        elseif(i.eq.nr)then
          thrra=(3.2d0*th(i+1,j,k)+2.d0*th(i-1,j,k)-0.2d0*th(i-2,j,k))*
     +          dri**2
          thra=(-th(i-1,j,k)/3.d0+4.d0*th(i+1,j,k)/3.d0)*dri
!mpk2d          thr=(-th(i-1,j,k)/3.d0-th(i,j,k)+4.d0*th(i+1,j,k)/3.d0)*dri
          thr=A(i,j,k)
          cnr=-2.5d0*dri**2-0.5d0*dri*ri
!22        elseif(i.eq.2 .or. i.eq.nr-1)then
         else
          thrra=(th(i-1,j,k)+th(i+1,j,k))*dri**2
          thra=(th(i+1,j,k)-th(i-1,j,k))*0.5d0*dri
!mpk2d          thr=thra
          thr=A(i,j,k)
c
          cnr=-1.d0*dri**2
!22        else
!mpk2d 4th order central diffrence
!22           thrra=(-0.083333333d0*th(i-2,j,k)+1.333333333d0*th(i-1,j,k)+
!22     +        1.333333333d0*th(i+1,j,k)-0.083333333d0*th(i+2,j,k))*
!22     +        dri**2
!22           thra=(0.0833333333d0*th(i-2,j,k)-0.6666666666d0*th(i-1,j,k)+
!22     +       0.6666666666d0*th(i+1,j,k)-0.0833333333d0*th(i+2,j,k))*
!22     +         dri
!22           thr=A(i,j,k)
!22           cnr=-1.25*dri**2
        endif
c
!mpk2d 4th order central diffrence  
        thhha=(th(i,j,ka)+th(i,j,km))*dpi**2
!22        thhha=(-0.083333333d0*th(i,j,kmm)+1.333333333d0*th(i,j,km)+
!22     +        1.333333333d0*th(i,j,ka)-0.083333333d0*th(i,j,kaa))*
!22     +        dpi**2
!mpk2d        thh=(th(i,j,ka)-th(i,j,km))*0.5d0*dpi
        thh=B(i,j,k)
c
!22         cnp=-(1.25d0*dpi**2)*ri**2
        cnp=-(1.d0*dpi**2)*ri**2
c
c  back to truly leapfrog/Dufort-Frankel
        RHS1=thrra+thra*ri+thhha*ri**2+0.5d0*thp(i,j,k)/dt  
        ULT=u*thr+v*thh*ri
        cf1=0.5d0/dt-cnr-cnp
        cf2=cnr+cnp
        thn(i,j,k)=(cf2*thp(i,j,k)-ULT+RHS1)/cf1

!mpk for Shapiro filter below    
        Ap(i,j,k)=thn(i,j,k)
c
  10  continue

!mpk Shapiro filter of order 2n=16, see Kalnay p. 102
!mpk \bar U^{2n}_{j} = [1-(-D^n)]U_j, where diffusion operator 
!mpk DU_j=0.25(U_{j+1}-2U_{j}+U_{j-1}) is applied to the original 
!mpk field n times, you may reduce n for high grid resolution run
 
!mpk commented to check 24 April 2006 
!       do norder=1,6
!mpk first applied it in r direction 
!         do k=0,np
!          do i=1,nr
!            An(i,j,k)=-0.25d0*(Ap(i+1,j,k)-2.d0*Ap(i,j,k)+Ap(i-1,j,k))
!          enddo
!         enddo

!         do k=0,np
!          do i=1,nr
!            Ap(i,j,k)=An(i,j,k)
!          enddo
!         enddo
 
!        enddo

!mpk stored filtered field in Ap
!        do i=1,nr
!         do k=0,np
!          Ap(i,j,k)=thn(i,j,k)-Ap(i,j,k)
!         enddo
!        enddo

!mpk then applied in tangential direction    
!        do norder=1,6        
!         do i=1,nr
!          do k=0,np
!           km=k-1
!           ka=k+1
!           if(k.eq.0)km=np-1
!           if(k.eq.np)ka=1          
!            An(i,j,k)=-0.25d0*(Ap(i,j,ka)-2.d0*Ap(i,j,k)+Ap(i,j,km))
!          enddo
!         enddo

!         do k=0,np
!          do i=1,nr
!            Ap(i,j,k)=An(i,j,k)
!          enddo
!         enddo

!        enddo
      
!        do i=1,nr
!         do k=0,np
!          thn(i,j,k)=thn(i,j,k)-Ap(i,j,k)
!         enddo
!        enddo

!mpk Robert filter
      do i=1,nr
       do k=0,np
        th(i,j,k)=th(i,j,k)+0.02d0*(thn(i,j,k)-2.d0*th(i,j,k)
     +            +thp(i,j,k))
       enddo
      enddo
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
      subroutine heat(trans1,trans2)
c
      implicit real*8(a-h,o-z)
      parameter(kr=290,kz=33,kp=1610,mgen=100000)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      common/sol2/thn,th,thp
      common/num/dt,dr,dz,dp,nr,nz,np
c
!mpk2d inner cylinder heat transfer
      trans1=0.d0
!mpk2d outer cylinder heat transfer
      trans2=0.d0
      dri=1.d0/dr
      j=1
      do 10 k=1,np
!mpk2d      trans=dz*(-184.d0*th(0,0,k)+225.d0*th(1,0,k)
!mpk2d     +     -50.d0*th(2,0,k)+9.d0*th(3,0,k))/(240.d0*dr)
!mpk2d      trans=trans+dz*(-184.d0*th(0,1,k)+225.d0*th(1,1,k)
!mpk2d     +     -50.d0*th(2,1,k)+9.d0*th(3,1,k))/(80.d0*dr)
!mpk2d      trans=trans+dz*(-184.d0*th(0,nz,k)+225.d0*th(1,nz,k)
!mpk2d     +     -50.d0*th(2,nz,k)+9.d0*th(3,nz,k))/(80.d0*dr)
!mpk2d      trans=trans+dz*(-184.d0*th(0,nz+1,k)+225.d0*th(1,nz+1,k)
!mpk2d     +     -50.d0*th(2,nz+1,k)+9.d0*th(3,nz+1,k))/(240.d0*dr)
c
!mpk2d      do 20 j=2,nz-1
!mpk2d        trans=trans+dz*(-184.d0*th(0,j,k)+225.d0*th(1,j,k)-50.d0
!mpk2d     +       *th(2,j,k)+9.d0*th(3,j,k))/(60.d0*dr)
!mpk2d  20  continue
!mpk2d        trans1=trans1+trans*dp
        ka=k+1
        if(k.eq.np)ka=1
        trans1=trans1+0.5d0*dp*(-3.35238095d0*th(0,j,k)+
     +         4.375d0*th(1,j,k)-1.45833333d0*th(2,j,k)+
     +         0.525d0*th(3,j,k)-0.08928571d0*th(4,j,k))*dri+
     +         0.5d0*dp*(-3.35238095d0*th(0,j,ka)+
     +         4.375d0*th(1,j,ka)-1.45833333d0*th(2,j,ka)+
     +         0.525d0*th(3,j,ka)-0.08928571d0*th(4,j,ka))*dri

         trans2=trans2+0.5d0*dp*(-3.35238095d0*th(nr+1,j,k)+
     +         4.375d0*th(nr,j,k)-1.45833333d0*th(nr-1,j,k)+
     +         0.525d0*th(nr-2,j,k)-0.08928571d0*th(nr-3,j,k))*dri+
     +         0.5d0*dp*(-3.35238095d0*th(nr+1,j,ka)+
     +         4.375d0*th(nr,j,ka)-1.45833333d0*th(nr-1,j,ka)+
     +         0.525d0*th(nr-2,j,ka)-0.08928571d0*th(nr-3,j,ka))*dri
         
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
c      subroutine aaaaaaaa(kt)
!mpkc
c      implicit real*8(a-h,o-z)
c      parameter(kr=161,kz=33,kp=800,mgen=100000)
c      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
c      real*8 ssumc(6,mgen),ssums(6,mgen)
c      common/grals/ssumc,ssums
c      common/sol2/thn,th,thp
c      common/num/dt,dr,dz,dp,nr,nz,np
!mpkc
c      do 10 n=1,6
c        ssumc(n,kt)=0.d0
c        ssums(n,kt)=0.d0 
c        sumc=0.d0
c        sums=0.d0
c        do 20 k=1,np
c          ang=dble(k)*dp
c          do 30 i=0,nr+1
c            r=r0+dble(i)*dr
c!mpk          do 30 j=0,nz+1
c            fact=1.d0
c            if(i.eq.0.or.i.eq.nr+1)fact=fact/2.d0
c!mpk            if(j.eq.0.or.j.eq.nz+1)fact=fact/2.d0
c            fc=r*th(i,1,k)*dcos(dble(n)*ang)
c            fs=r*th(i,1,k)*dsin(dble(n)*ang)
c            sumc=sumc+fc*fact
c            sums=sums+fs*fact
c   30     continue
c!mpk2d          sumc=sumc*dr*dz
c!mpk2d         sums=sums*dr*dz
c          sumc=sumc*dr
c          sums=sums*dr
c!mpkc
c          ssumc(n,kt)=ssumc(n,kt)+sumc*dp
c          ssums(n,kt)=ssums(n,kt)+sums*dp
c   20   continue
c   10 continue
c!mpkc
c      return
c      end
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
      parameter(kr=290,kz=33,kp=1610,mgen=100000)
      real*8 An(0:kr,1,0:kp),A(0:kr,1,0:kp),Ap(0:kr,1,0:kp)
      real*8 Bn(0:kr,1,0:kp),B(0:kr,1,0:kp),Bp(0:kr,1,0:kp)
      real*8 Cn(0:kr,1,0:kp),C(0:kr,1,0:kp),Cp(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 uelp(0:kr,1,0:kp),velp(0:kr,1,0:kp)
      real*8 welp(0:kr,1,0:kp)
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c------left and right side walls (omitting edges)
!mpk2d      do 10 i=1,nr
!mpk2d      do 10 k=0,np
!mpk2d         An(i,0,k)=-(3.d0*vel(i,1,k)-vel(i,2,k)/3.d0)/dz
!mpk2d         An(i,nz,k)=-(vel(i,nz-1,k)/3.d0-3.d0*vel(i,nz,k))/dz
!mpk2d  10  continue
!mpk2d       do 11 i=1,nr-1
!mpk2d       do 11 k=0,np
!mpk2d         Bn(i,0,k)=(3.d0*uel(i,1,k)-uel(i,2,k)/3.d0)/dz
!mpk2d         Bn(i,nz,k)=(-3.d0*uel(i,nz,k)+uel(i,nz-1,k)/3.d0)/dz
!mpk2d   11  continue
c------inner and outer cylindrical walls (omitting edges)
      j = 1
!mpk2d      do 12 j=1,nz-1
      do 12 k=0,np
        Bn(0,j,k)=-(3.d0*wel(1,j,k)-wel(2,j,k)/3.d0)/dr
        Bn(nr,j,k)=-(wel(nr-1,j,k)/3.d0-3.d0*wel(nr,j,k))/dr
  12  continue 
!mpk2d      do 13 j=1,nz
      do 13 k=0,np
        Cn(0,j,k)=(3.d0*vel(1,j,k)-vel(2,j,k)/3.d0)/dr
!mpktc        Cn(nr,j,k)=(vel(nr-1,j,k)/3.d0-3.d0*vel(nr,j,k))/dr
        Cn(nr,j,k)=(8.d0/3.d0*vel(nr+1,j,k)+vel(nr-1,j,k)/3.d0-3.d0*
     &            vel(nr,j,k))/dr+vel(nr+1,j,k)/rmax
  13  continue
c
!mpk2d      do 30 j=0,nz+1
      do 30 k=0,np
!mpktc        thn(0,j,k)=0.d0
         thn(0,j,k)=0.d0
!mpktc        thn(nr+1,j,k)=1.d0
         thn(nr+1,j,k)=1.d0
  30  continue
c
      j=1
      do 40 i=0,nr+1
!mpk2d      do 40 j=0,nz+1
      do 40 k=0,np
!mpk2d         Ap(i,j,k)=A(i,1,k)
!mpk2d         A(i,j,k)=An(i,1,k)
!mpk2d         Bp(i,j,k)=B(i,1,k)
!mpk2d         B(i,j,k)=Bn(i,1,k)
        Cp(i,j,k)=C(i,1,k)
        C(i,j,k)=Cn(i,1,k)
        thp(i,j,k)=th(i,1,k)
        if (thn(i,1,k) .gt. 1.d0) then 
          thn(i,1,k)=0.99d0
          print*,'WARNING! Temp damped'
        endif
        if (thn(i,1,k) .lt. 0.d0) then
          thn(i,1,k)=0.01d0
          print*,'WARNING! Temp damped'
        endif
        th(i,j,k)=thn(i,1,k)
  40  continue
c
      return
      end
c
c----------------------------------------------------------------------c
      subroutine rhsv(vrh)
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(kr=290,kz=33,kp=1610,mgen=100000)
      real*8 An(0:kr,1,0:kp),A(0:kr,1,0:kp),Ap(0:kr,1,0:kp)
      real*8 Bn(0:kr,1,0:kp),B(0:kr,1,0:kp),Bp(0:kr,1,0:kp)
      real*8 Cn(0:kr,1,0:kp),C(0:kr,1,0:kp),Cp(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 uelp(0:kr,1,0:kp),velp(0:kr,1,0:kp)
      real*8 welp(0:kr,1,0:kp)
      real*8 vrh(0:kr,1,0:kp)
!mpktc
      common/val/Ra,Re,Ptl,Gr
!mpktc
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      j=1
      do 10 i=1,nr
!mpk2d       do 10 j=1,nz
!mpktc
       hri=1.d0/dr
       hhri=hri**2.d0
!mpktc
      do 10 k=0,np
        r=r0+(dble(i)-0.5d0)*dr
!mpk2d        Az=(An(i,j,k)-An(i,j-1,k))/dz
        Cr=(Cn(i,j,k)-Cn(i-1,j,k))/dr
        if(i.eq.1)then
          vrh(i,j,k)=Cn(1,j,k)/dr+Cn(1,j,k)/r
        elseif(i.eq.nr)then
          vrh(i,j,k)=-Cn(nr-1,j,k)/dr+Cn(nr-1,j,k)/r
!mpktc
          velo=-(Gr**0.5d0)*Ptl*0.5d0*rmax
          vrh(i,j,k)=vrh(i,j,k)-8.d0/15.d0*velo*hhri+rvali*hri*velo-
     +     4.d0/3.d0*velo*hri*rvali+velo*rvali**2.d0
!mpktc
        else
          vrh(i,j,k)=Cr+(Cn(i,j,k)+Cn(i-1,j,k))/r
        endif
!mpk2d        if(j.eq.1)then
!mpk2d          vrh(i,j,k)=vrh(i,j,k)-An(i,1,k)/dz
!mpk2d        elseif(j.eq.nz)then
!mpk2d          vrh(i,j,k)=vrh(i,j,k)+An(i,nz-1,k)/dz
!mpk2d        else
!mpk2d          vrh(i,j,k)=vrh(i,j,k)-Az
!mpk2d        endif
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
      parameter(kr=290,kz=33,kp=1610,mgen=100000)
      real*8 An(0:kr,1,0:kp),A(0:kr,1,0:kp),Ap(0:kr,1,0:kp)
      real*8 Bn(0:kr,1,0:kp),B(0:kr,1,0:kp),Bp(0:kr,1,0:kp)
      real*8 Cn(0:kr,1,0:kp),C(0:kr,1,0:kp),Cp(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp),th(0:kr,1,0:kp),thp(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 uelp(0:kr,1,0:kp),velp(0:kr,1,0:kp)
      real*8 welp(0:kr,1,0:kp)
      real*8 urh(0:kr,1,0:kp)
      common/sol1/An,A,Ap,Bn,B,Bp,Cn,C,Cp
      common/sol2/thn,th,thp
      common/sol3/uel,vel,wel,uelp,velp,welp
      common/num/dt,dr,dz,dp,nr,nz,np
      common/radius/r0,rmax
c
      j=1
      do 10 i=1,nr-1
!mpk2d      do 10 j=1,nz
      do 10 k=0,np
        km=k-1
        if(k.eq.0)km=np-1
        r=r0+dble(i)*dr
!mpk2d        Bz=(Bn(i,j,k)-Bn(i,j-1,k))/dz
!mpk2d        welz=0.5d0*(wel(i,j,k)+wel(i+1,j,k)-wel(i,j-1,k)
!mpk2d     +      -wel(i+1,j-1,k))/dz
        Ch=(Cn(i,j,k)-Cn(i,j,km))/dp
!mpk2d        if(j.eq.1)then
!mpk2d          urh(i,j,k)=Bn(i,1,k)/dz-Ch/r-2.d0*welz/r
!mpk2d        elseif(j.eq.nz)then
!mpk2d          urh(i,j,k)=-Bn(i,nz-1,k)/dz-Ch/r-2.d0*welz/r
!mpk2d        else
!mpk2d          urh(i,j,k)=Bz-Ch/r-2.d0*welz/r
!mpk2d        endif
         urh(i,j,k)=-Ch/r
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
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,1,0:kp),psif(0:kr,1,0:kp)
      real*8 rc(0:kr,1,0:kp),rf(0:kr,1,0:kp),tempa(0:kr,1,0:kp)
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
!mpk2d      jtest=nz
!mpk2d      jtest2=nz
c     print*,itest,itest2,jtest,jtest2
      do 1 i=2,8
        itest=itest/2
!mpk2d        jtest=jtest/2
!mpk2d        if(itest*2-itest2.eq.0.and.jtest*2-jtest2.eq.0)then
         if(itest*2-itest2.eq.0)then
          if(itest2.ge.4)nlev=nlev+1
          itest2=itest
!mpk2d          jtest2=jtest
        else
          goto 2
        end if
    1 continue
    2 continue
      print*,'possible number of MG levels=',nlev
!mpk2d high resolution so more levels should be fine but i try to limit
!it because there may be a numerical instability issue.
      if(nlev.ge.4)nlev=4
      print*,'but i limit the number of MG levels=',nlev
c
      ii=1
      do 20 i=1,nlev
        n1(i)=nr*2/(2**i)
!mpk2d         n2(i)=nz*2/(2**i)
        n2(i)=nz
        n3(i)=np*2/(2**i)
        hr(i)=dr*dble(2**i)*0.5d0
!mpk2d        hz(i)=dz*dble(2**i)*0.5d0
        hz(i)=dz
        hp(i)=dp*dble(2**i)*0.5d0
        hhr(i)=hr(i)**2
!mpk2d        hhz(i)=hz(i)**2
        hhz(i)=hz(1)**2
        hhp(i)=hp(i)**2
        ibeg(i)=ii
!mpk2d        ii=ii+(n1(i)+1)*(n2(i)+1)*(n3(i)+1)+1
        ii=ii+(n1(i)+1)*(n3(i)+1)+1  
        
c       print*,i
c       print*,n1(i),n2(i),n3(i),ibeg(i)
c       print*,hr(i),hz(i),hp(i)
   20 continue
c
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
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 psi(0:kr,1,0:kp),psip(0:kr,1,0:kp)
      real*8 rhs(0:kr,1,0:kp)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,1,0:kp),psif(0:kr,1,0:kp)
      real*8 rc(0:kr,1,0:kp),rf(0:kr,1,0:kp),tempa(0:kr,1,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
!$OMP THREADPRIVATE(/backup/,/arrays/)

c
      nvs=10000
      nrel1=1
      nrel2=1
c
      do 10 i=1,mgen
        v(i)=0.d0
        f(i)=0.d0
   10 continue
c
      j=1
      do 30 i=0,nr+1
!mpk2d      do 30 j=0,nz+1
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
!mpk2d        do 31 js=0,nz+1
        do 31 ks=0,np
        js=1
        psif(is,js,ks)=psip(is,js,ks)+1.2d0*(psif(is,js,ks)-
     +   psip(is,js,ks))
        psip(is,js,ks)=psif(is,js,ks)
   31   continue
        call trans1(psif,v,ibeg(1),n1(1),n2(1),n3(1))
      endif
        if(i.eq.1.or.nlev.eq.1)then
          call relax2(nrel1,1,nconv,nlev)
          do 55 ii=0,nr+1
!mpk2d            do 55 jj=0,nz+1
             jj=1
            do 55 kk=0,np
              psip(ii,jj,kk)=psif(ii,jj,kk)
   55     continue
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser(ilev)
c            print*,'Sweep beg',ilev
            call relax2(nrel1,ilev+1,nconv,nlev)
c            print*,'Sweep end',ilev
   51     continue

!mpk2d introducing FMG, works well, see Briggs et al.'s Tutorial
          do ifmg=0,nlev-2

           do ilev=nlev,nlev-ifmg,-1
             call finer(ilev)
             call relax2(nrel2,ilev-1,nconv,nlev)
           enddo 
           do ilev=nlev-ifmg-1,nlev-1
              call corser(ilev)
              call relax2(nrel1,ilev+1,nconv,nlev)
           enddo

          enddo

          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax2(nrel2,ilev-1,nconv,nlev)
   52     continue
        end if
        if(nconv.eq.1)then
c         print*,'Converged'
          do 64 ii=0,nr+1
!mpk2d          do 64 jj=0,nz+1
          jj=1
          do 64 kk=0,np
            psip(ii,jj,kk)=psi(ii,jj,kk)
   64     continue
c       UPDATING
c       internal points
          do 60 ii=1,nr
!mpk2d          do 60 jj=1,nz
          jj=1
          do 60 kk=0,np
            psi(ii,jj,kk)=psif(ii,jj,kk)
   60     continue
c       cylinders
!mpk2d          do 62 jj=0,nz+1
          jj=1
          do 62 kk=0,np
            psi(0,jj,kk)=0.d0        
            psi(nr+1,jj,kk)=0.d0     
   62     continue
c       endwalls
!mpk2d          do 63 ii=1,nr
!mpk2d          do 63 kk=0,np
!mpk2d            psi(ii,0,kk)=0.d0        
!mpk2d            psi(ii,nz+1,kk)=0.d0
!mpk2d   63     continue    
!mpk2d         print*,'v vel; sweeps=',i
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
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 psi(0:kr,1,0:kp),psip(0:kr,1,0:kp)
      real*8 rhs(0:kr,1,0:kp)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,1,0:kp),psif(0:kr,1,0:kp)
      real*8 rc(0:kr,1,0:kp),rf(0:kr,1,0:kp),tempa(0:kr,1,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
!$OMP THREADPRIVATE(/backup/,/arrays/)

c
      nvs=10000
      nrel1=1
      nrel2=1
c
      do 10 i=1,mgen
        v(i)=0.d0
        f(i)=0.d0
   10 continue
c
      do 30 i=0,nr
!mpk2d      do 30 j=0,nz+1 
      j=1
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
!mpk2d        do 31 js=0,nz+1
         js=1
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
!mpk2d            do 54 jj=0,nz+1
            jj=1
            do 54 kk=0,np
             psip(ii,jj,kk)=psif(ii,jj,kk)
   54       continue
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser(ilev)
            call relax3(nrel1,ilev+1,nconv,nlev)
   51     continue

!mpk2d FMG see Briggs et al. Tutorial
         do ifmg=0,nlev-2

           do ilev=nlev,nlev-ifmg,-1
             call finer(ilev)
             call relax3(nrel2,ilev-1,nconv,nlev)
           enddo 
           do ilev=nlev-ifmg-1,nlev-1
              call corser(ilev)
              call relax3(nrel1,ilev+1,nconv,nlev)
           enddo

         enddo

          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax3(nrel2,ilev-1,nconv,nlev)
   52     continue
        end if
        if(nconv.eq.1)then
c         print*,'Converged'
          do 63 ii=0,nr
!mpk2d          do 63 jj=0,nz+1
          jj=1
          do 63 kk=0,np
            psip(ii,jj,kk)=psi(ii,jj,kk)
   63     continue
c       UPDATING
c       internal points
          do 60 ii=1,nr-1
!mpk2d          do 60 jj=1,nz
           jj=1
          do 60 kk=0,np
            psi(ii,jj,kk)=psif(ii,jj,kk)
   60     continue
c       cylinders
!mpk2d          do 61 jj=0,nz+1
           jj=1
          do 61 kk=0,np
            psi(0,jj,kk)=0.d0
            psi(nr,jj,kk)=0.d0
   61     continue
c       endwalls
!mpk2d          do 62 ii=0,nr
!mpk2d          do 62 kk=0,np
!mpk2d             psi(ii,0,kk)=0.d0
!mpk2d             psi(ii,nz+1,kk)=0.d0
!mpk2d   62     continue
!mpk2d          print*,'u vel; sweeps=',i
c         call rufcon(psi,nr,nz,np)
          return
        end if
   50 continue
      return
      end
c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      subroutine relax2(nrel,klev,nconv,nlev)
c
      implicit real*8(a-h,o-z)
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,1,0:kp),psif(0:kr,1,0:kp)
      real*8 rc(0:kr,1,0:kp),rf(0:kr,1,0:kp),tempa(0:kr,1,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      real*8 pdmaa(kp),pdmab(kp),pdmac(kp),pdmad(kp),pdmae(kp),pdmas(kp)
      real*8 pdmar(kp)
c8aug06
      real*8 rvali(kr),rvalisq(kr)
c8aug06
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
!$OMP THREADPRIVATE(/backup/,/arrays/)

c
      tolcon=1.d-5
      nconv=0
c
      call trans2(f,rf,ibeg(klev),n1(klev),n2(klev),n3(klev))
c      print*,'before ',psif(1,1,1)
c14aug06      call trans2(v,psif,ibeg(klev),n1(klev),n2(klev),n3(klev))
c       if(klev.eq.2)then
c       do k=1,n3(klev)
c         do i=1,n1(klev)
c           print*,i,k,klev,psic(i,1,k),psif(i,1,k)
c          if(psic(i,1,k).ne.psif(i,1,k))then
c            print*,i,k,klev,psic(i,1,k),psif(i,1,k)
c            STOP
c          endif
c         enddo

c        enddo
c       endif
c      print*,'afterr ',psif(1,1,1)
c
c------------------------------relaxing with some G & S.
c
      hhri=1./hhr(klev)
      hri=1./hr(klev)
      hhpi=1./hhp(klev)
      do 200 irel=1,nrel
c                                 line relax in r-dir
        do 10 k=1,n3(klev)
!mpk2d        do 10 j=1,n2(klev)
           j=1
          ka=k+1
          km=k-1
          if(k.eq.0)km=n3(klev)-1
          if(k.eq.n3(klev))ka=1
            do 11 i=1,n1(klev)
c8aug06            rval=r0+(dble(i)-0.5d0)*hr(klev)  
c8aug06            rvali=1.d0/rval
            if(k.eq.1)then
              rval=r0+(dble(i)-0.5d0)*hr(klev)
              rvali(i)=1.d0/rval
              rvalisq(i)=rvali(i)**2.
            if(i.eq.1)then
              pdmaa(i)=0.d0
              pdmab(i)=0.d0
              pdmac(i)=-2.d0*hhri-2.d0*hhpi*rvalisq(i)
     +                 +rvalisq(i)       !mpk +1.d0/rval/hr(klev)
              pdmad(i)=5.d0/3.d0*hhri+4.d0/3.d0*hri*rvali(i)
              pdmae(i)=-0.2d0*hhri
            elseif(i.eq.n1(klev))then
              pdmaa(i)=-0.2d0*hhri
              pdmab(i)=5.d0/3.d0*hhri-4.d0/3.d0*hri*rvali(i)
              pdmac(i)=-2.d0*hhri-2.d0*hhpi*rvalisq(i)
     +                +rvalisq(i)        !mpk -1.d0/rval/hr(klev)
              pdmad(i)=0.d0
              pdmae(i)=0.d0
            else
              pdmaa(i)=0.d0
              pdmab(i)=hhri-1.5d0*hri*rvali(i)
              pdmac(i)=-(2.d0*hhri+2.d0*hhpi*rvalisq(i))
     +                +rvalisq(i)
              pdmad(i)=hhri+1.5d0*hri*rvali(i)
              pdmae(i)=0.d0
            endif
            endif
              tdmaa(i)=pdmaa(i)
              tdmab(i)=pdmab(i)
              tdmac(i)=pdmac(i)
              tdmar(i)=pdmad(i)
              tdmas(i)=pdmae(i) 
c
!mpk2d            if(j.eq.1)then
!mpk2d             pdmac(i)=pdmac(i)-2.d0/hhz(klev)
!mpk2d             pdmar(i)=-(5.d0/3.d0*psif(i,j+1,k)-0.2*psif(i,j+2,k))
!mpk2d     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
!mpk2d     +               /rval**2+rf(i,j,k)
!mpk2d            elseif(j.eq.n2(klev))then
!mpk2d             pdmac(i)=pdmac(i)-2.d0/hhz(klev)
!mpk2d             pdmar(i)=-(5.d0/3.d0*psif(i,j-1,k)-0.2*psif(i,j-2,k))
!mpk2d     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
!mpk2d     +               /rval**2+rf(i,j,k)
!mpk2d            else
!mpk2d             pdmac(i)=pdmac(i)-2.d0/hhz(klev)
!mpk2d             pdmar(i)=-(psif(i,j-1,k)+psif(i,j+1,k))
!mpk2d     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
!mpk2d     +               /rval**2+rf(i,j,k)
!mpk2d            endif
c
c8aug06                 pdmac(i)=pdmac(i)
               pdmar(i)=-(psif(i,j,ka)+psif(i,j,km))*hhpi
     +               *rvalisq(i)+rf(i,j,k)

   11     continue
          call pdma(tdmaa,tdmab,tdmac,tdmar,tdmas,pdmar,pdmas,n1(klev))
          do 12 i=1,n1(klev)
            psif(i,j,k)=pdmas(i)
             if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
   12     continue
   10   continue
c                                 line relax in z-dir
!mpk2d        do 13 k=1,n3(klev)
!mpk2d        do 13 i=1,n1(klev)
!mpk2d            rval=r0+(dble(i)-0.5d0)*hr(klev)
!mpk2d            rvali=1.d0/rval
!mpk2d          ka=k+1
!mpk2d          km=k-1
!mpk2d          if(k.eq.0)km=n3(klev)-1
!mpk2d          if(k.eq.n3(klev))ka=1
!mpk2d          do 14 j=1,n2(klev)
!mpk2d            if(j.eq.1)then
!mpk2d              pdmaa(j)=0.d0
!mpk2d              pdmab(j)=0.d0
!mpk2d              pdmac(j)=-2.d0/hhz(klev)-2.d0/hhp(klev)/rval**2
!mpk2d     +                +1.d0/rval**2
!mpk2d              pdmad(j)=5.d0/3.d0/hhz(klev)
!mpk2d              pdmae(j)=-0.2d0/hhz(klev)
!mpk2d            elseif(j.eq.n2(klev))then
!mpk2d              pdmaa(j)=-0.2d0/hhz(klev)
!mpk2d              pdmab(j)=5.d0/3.d0/hhz(klev)
!mpk2d              pdmac(j)=-2.d0/hhz(klev)-2.d0/hhp(klev)/rval**2
!mpk2d     +                +1.d0/rval**2
!mpk2d              pdmad(j)=0.d0
!mpk2d              pdmae(j)=0.d0
!mpk2d            else
!mpk2d              pdmaa(j)=0.d0
!mpk2d              pdmab(j)=1.d0/hhz(klev)
!mpk2d              pdmac(j)=-(2.d0/hhz(klev)+2.d0/hhp(klev)/rval**2)
!mpk2d     +                +1.d0/rval**2
!mpk2d              pdmad(j)=1.d0/hhz(klev)
!mpk2d              pdmae(j)=0.d0
!mpk2d            endif
c
!mpk2d            if(i.eq.1)then
!mpk2d            pdmac(j)=pdmac(j)-2.d0/hhr(klev)     !mpk +1.d0/rval/hr(klev)
!mpk2d            pdmar(j)=-(5.d0/3.d0*psif(i+1,j,k)-0.2d0*psif(i+2,j,k))/hhr
!mpk2d     +                (klev)-4.d0/3.d0*psif(i+1,j,k)/hr(klev)/rval
!mpk2d     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
!mpk2d     +                +rf(i,j,k)
!mpk2d            elseif(i.eq.n1(klev))then
!mpk2d            pdmac(j)=pdmac(j)-2.d0/hhr(klev)      !mpk -1.d0/rval/hr(klev)
!mpk2d            pdmar(j)=-(5.d0/3.d0*psif(i-1,j,k)-0.2d0*psif(i-2,j,k))/
!mpk2d     +                hhr(klev)+4.d0/3.d0*psif(i-1,j,k)/hr(klev)/rval
!mpk2d     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
!mpk2d     +                +rf(i,j,k)
!mpk2d            else
!mpk2d              pdmac(j)=pdmac(j)-2.d0/hhr(klev)
!mpk2d              pdmar(j)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)-1.5d0*
!mpk2d     +                (psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)/rval
!mpk2d     +                -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)/rval**2
!mpk2d     +                +rf(i,j,k)
!mpk2d            endif
!mpk2d   14     continue
!mpk2d          call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n2(klev))
!mpk2d          do 15 j=1,n2(klev)
!mpk2d            psif(i,j,k)=pdmas(j)
!mpk2d            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
!mpk2d   15     continue
!mpk2d   13   continue

!mpk2d apparently is not influecing MG sweeps. so comment out to save cost
c                             line relax in phi-dir
!mpk2d         do 16 j=1,n2(klev)
!mpk2d        j=1   
!mpk2d        do 16 i=1,n1(klev)
!mpk2d            rval=r0+(dble(i)-0.5d0)*hr(klev)
!mpk2d            rvali=1.d0/rval
!mpk2d          do 17 k=1,n3(klev)
!mpk2d            tdmaa(k)=hhpi*rvali**2
!mpk2d            tdmac(k)=hhpi*rvali**2
!mpk2d            tdmab(k)=-(2.d0*hhpi*rvali**2)+rvali**2
!mpk2d            if(i.eq.1)then
!mpk2d            tdmab(k)=tdmab(k)-2.d0*hhri      !mpk +1.d0/rval/hr(klev)
!mpk2d            tdmar(k)=-(5.d0/3.d0*psif(i+1,j,k)-0.2d0*psif(i+2,j,k))*
!mpk2d     +                hhri-4.d0/3.d0*psif(i+1,j,k)*hri*rvali
!mpk2d            elseif(i.eq.n1(klev))then
!mpk2d            tdmab(k)=tdmab(k)-2.d0*hhri       !mpk -1.d0/rval/hr(klev)
!mpk2d            tdmar(k)=-(5.d0/3.d0*psif(i-1,j,k)-0.2d0*psif(i-2,j,k))
!mpk2d     +              *hhri+4.d0/3.d0*psif(i-1,j,k)*hri*rvali
!mpk2d           else
!mpk2d             tdmab(k)=tdmab(k)-2.d0*hhri
!mpk2d             tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))*hhri-1.5d0*
!mpk2d     +                (psif(i+1,j,k)-psif(i-1,j,k))*hri*rvali
!mpk2d            endif
c
!mpk2d!mpk2d            if(j.eq.1)then
!mpk2d!mpk2d             tdmab(k)=tdmab(k)-2.d0/hhz(klev)
!mpk2d!mpk2d             tdmar(k)=tdmar(k)-(5.d0/3.d0*psif(i,j+1,k)
!mpk2d!mpk2d     +               -0.2d0*psif(i,j+2,k))/hhz(klev)+rf(i,j,k)
!mpk2d!mpk2d            elseif(j.eq.n2(klev))then
!mpk2d!mpk2d             tdmab(k)=tdmab(k)-2.d0/hhz(klev)
!mpk2d!mpk2d             tdmar(k)=tdmar(k)-(5.d0/3.d0*psif(i,j-1,k)
!mpk2d!mpk2d     +               -0.2d0*psif(i,j-2,k))/hhz(klev)+rf(i,j,k)
!mpk2d!mpk2d            else
!mpk2d!mpk2d             tdmab(k)=tdmab(k)-2.d0/hhz(klev)
!mpk2d!mpk2d             tdmar(k)=tdmar(k)-(psif(i,j-1,k)+psif(i,j+1,k))
!mpk2d!mpk2d     +               /hhz(klev)+rf(i,j,k)
!mpk2d!mpk2d            endif
!mpk2d                  tdmab(k)=tdmab(k)
!mpk2d                  tdmar(k)=tdmar(k)+rf(i,j,k)
c
!mpk2d   17     continue
!mpk2d          call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n3(klev))
!mpk2d          do 18 k=0,n3(klev)
!mpk2d            if(k.eq.0)then
!mpk2d              psif(i,j,k)=tdmas(n3(klev))
!mpk2d            else
!mpk2d              psif(i,j,k)=tdmas(k)
!mpk2d            endif
!mpk2d   18     continue
!mpk2d   16   continue
c
c                                  Calculating residual.
c
          resmax=0.d0
          do 20 k=1,n3(klev)
          do 20 i=1,n1(klev)
c8aug06             rval=r0+(dble(i)-0.5d0)*hr(klev)
c8aug06             rvali=1.d0/rval
             j=1
!mpk2d          do 20 j=1,n2(klev)
c8aug06         do 20 k=1,n3(klev)
          ka=k+1
          km=k-1
          if(k.eq.0)km=n3(klev)-1
          if(k.eq.n3(klev))ka=1
          temp=(psif(i,j,ka)-2.d0*psif(i,j,k)+psif(i,j,km))
     +        *hhpi*rvalisq(i)+psif(i,j,k)*rvalisq(i)
          if(i.eq.1)then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0/3.d0*psif(i+1,j,k)-
     +           0.2d0*psif(i+2,j,k))*hhri+
     +           4.d0/3.d0*psif(i+1,j,k)*hri*rvali(i)
          elseif(i.eq.n1(klev))then
            temp=temp+(-2.d0*psif(i,j,k)+5.d0/3.d0*psif(i-1,j,k)
     +          -0.2d0*psif(i-2,j,k))*hhri-
     +          4.d0/3.d0*psif(i-1,j,k)*hri*rvali(i)
          else
            temp=temp+(psif(i-1,j,k)-2.d0*psif(i,j,k)+psif(i+1,j,k))
     +          *hhri+1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))
     +          *hri*rvali(i)
          endif
c
!mpk2d         if(j.eq.1)then
!mpk2d            temp=temp+(-2.d0*psif(i,j,k)+5.d0/3.d0*psif(i,j+1,k)
!mpk2d     +          -0.2d0*psif(i,j+2,k))/hhz(klev)
!mpk2d          elseif(j.eq.n2(klev))then
!mpk2d            temp=temp+(-2.d0*psif(i,j,k)+5.d0/3.d0*psif(i,j-1,k)
!mpk2d     +          -0.2d0*psif(i,j-2,k))/hhz(klev)
!mpk2d          else
!mpk2d            temp=temp+(psif(i,j+1,k)-2.d0*psif(i,j,k)+psif(i,j-1,k))
!mpk2d     +          /hhz(klev)
!mpk2d          endif
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
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,1,0:kp),psif(0:kr,1,0:kp)
      real*8 rc(0:kr,1,0:kp),rf(0:kr,1,0:kp),tempa(0:kr,1,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      real*8 tdmaa(kp),tdmab(kp),tdmac(kp),tdmar(kp),tdmas(kp)
      real*8 pdmaa(kp),pdmab(kp),pdmac(kp),pdmad(kp),pdmae(kp),pdmas(kp)
      real*8 pdmar(kp)
c8aug06
      real*8 rvali(kr),rvalisq(kr)
c8aug06
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
      common /radius/r0,rmax
      common/iteration/kt
!$OMP THREADPRIVATE(/backup/,/arrays/)

c
      tolcon=1.d-5
      nconv=0
c
      call trans2(f,rf,ibeg(klev),n1(klev),n2(klev),n3(klev))
c14aug06      call trans2(v,psif,ibeg(klev),n1(klev),n2(klev),n3(klev))
c
c------------------------------relaxing with some G & S.
c
      hhri=1./hhr(klev)
      hri=1./hr(klev)
      hhpi=1./hhp(klev)
      do 200 irel=1,nrel
c                                 line relax in r-dir        
        do 10 k=1,n3(klev)
!mpk2d        do 10 j=1,n2(klev)
          j=1
          ka=k+1
          km=k-1
          if(k.eq.n3(klev))ka=1
          do 11 i=1,n1(klev)-1
c8aug06            rval=dble(i)*hr(klev)+r0
c8aug06            rvali=1./rval
            if(k.eq.1)then
             rval=dble(i)*hr(klev)+r0
             rvali(i)=1./rval
             rvalisq(i)=rvali(i)**2.
             pdmaa(i)=hhri-1.5d0*hri*rvali(i)
             pdmab(i)=-(2.d0*hhri-rvalisq(i)+2.d0*hhpi*rvalisq(i))
             pdmac(i)=hhri+1.5d0*hri*rvali(i)
            endif
             tdmaa(i)=pdmaa(i)
             tdmab(i)=pdmab(i)
             tdmac(i)=pdmac(i)
!mpk2d            if(j.eq.1)then
!mpk2d              tdmab(i)=tdmab(i)-2.d0/hhz(klev)
!mpk2d              tdmar(i)=-(5.d0*psif(i,j+1,k)/3.d0-0.2*psif(i,j+2,k))
!mpk2d     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
!mpk2d     +               /rval**2+rf(i,j,k)
!mpk2d            elseif(j.eq.n2(klev))then
!mpk2d              tdmab(i)=tdmab(i)-2.d0/hhz(klev)
!mpk2d              tdmar(i)=-(5.d0*psif(i,j-1,k)/3.d0-0.2*psif(i,j-2,k))
!mpk2d     +               /hhz(klev)-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
!mpk2d     +               /rval**2+rf(i,j,k)
!mpk2d            else
!mpk2d              tdmab(i)=tdmab(i)-2.d0/hhz(klev)
!mpk2d              tdmar(i)=-(psif(i,j-1,k)+psif(i,j+1,k))/hhz(klev)
!mpk2d     +               -(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
!mpk2d     +               /rval**2+rf(i,j,k)
!mpk2d            endif
c8aug06             tdmab(i)=tdmab(i)
              tdmar(i)=-(psif(i,j,ka)+psif(i,j,km))*hhpi*rvalisq(i)
     +                 +rf(i,j,k)
   11     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n1(klev)-1)
          do 12 i=1,n1(klev)-1
            psif(i,j,k)=tdmas(i)
            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
   12     continue
   10   continue
c                                 line relax in z-dir
!mpk2d        do 13 k=1,n3(klev)
!mpk2d        do 13 i=1,n1(klev)-1
!mpk2d            rval=dble(i)*hr(klev)+r0
!mpk2d          ka=k+1
!mpk2d          km=k-1
!mpk2d          if(k.eq.0)km=n3(klev)-1
!mpk2d          if(k.eq.n3(klev))ka=1
!mpk2d          do 14 j=1,n2(klev)
!mpk2d            if(j.eq.1)then
!mpk2d              pdmaa(j)=0.d0
!mpk2d              pdmab(j)=0.d0
!mpk2d              pdmac(j)=-2.d0/hhz(klev)-2.d0/hhr(klev)+1.d0/rval**2
!mpk2d     +                -2.d0/hhp(klev)/rval**2
!mpk2d              pdmad(j)=5.d0/hhz(klev)/3.d0
!mpk2d              pdmae(j)=-0.2d0/hhz(klev)
!mpk2d            elseif(j.eq.n2(klev))then
!mpk2d              pdmaa(j)=-0.2d0/hhz(klev)
!mpk2d              pdmab(j)=5.d0/hhz(klev)/3.d0
!mpk2d              pdmac(j)=-2.d0/hhz(klev)-2.d0/hhr(klev)+1.d0/rval**2
!mpk2d     +                -2.d0/hhp(klev)/rval**2
!mpk2d              pdmad(j)=0.d0
!mpk2d              pdmae(j)=0.d0
!mpk2d            else
!mpk2d              pdmaa(j)=0.d0
!mpk2d              pdmab(j)=1.d0/hhz(klev)
!mpk2d              pdmac(j)=-(2.d0/hhz(klev)+2.d0/hhr(klev)+2.d0/hhp(klev)
!mpk2d     +              /rval**2)+1.d0/rval**2
!mpk2d              pdmad(j)=1.d0/hhz(klev)
!mpk2d              pdmae(j)=0.d0
!mpk2d            endif
!mpk2d            pdmar(j)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
!mpk2d     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
!mpk2d     +               /rval-(psif(i,j,ka)+psif(i,j,km))/hhp(klev)
!mpk2d     +               /rval**2+rf(i,j,k)
!mpk2d   14     continue
!mpk2d          call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n2(klev))
!mpk2d          do 15 j=1,n2(klev)
!mpk2d            psif(i,j,k)=pdmas(j)
!mpk2d            if(k.eq.n3(klev))psif(i,j,0)=psif(i,j,k)
!mpk2d   15     continue
!mpk2d   13   continue

!mpk2d apparently is not influecing MG sweeps. so comment out to save cost
c                                line relax in phi-dir
!mpk2d        do 16 j=1,n2(klev)
!mpk2d        j=1
!mpk2d        do 16 i=1,n1(klev)-1
!mpk2d            rval=dble(i)*hr(klev)+r0
!mpk2d            rvali=1.d0/rval
!mpk2d          do 17 k=1,n3(klev)
!mpk2d            tdmaa(k)=hhpi*rvali**2
!mpk2d            tdmab(k)=-(2.d0*hhri-rvali**2+2.d0*hhpi*rvali**2)
!mpk2d            tdmac(k)=hhpi*rvali**2
!mpk2d!mpk2d           if(j.eq.1)then
!mpk2d!mpk2d              tdmab(k)=tdmab(k)-2.d0/hhz(klev)
!mpk2d!mpk2d              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
!mpk2d!mpk2d     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
!mpk2d!mpk2d     +               /rval-(5.d0*psif(i,j+1,k)/3.d0
!mpk2d!mpk2d     +               -0.2*psif(i,j+2,k))/hhz(klev)+rf(i,j,k)
!mpk2d!mpk2d            elseif(j.eq.n2(klev))then
!mpk2d!mpk2d              tdmab(k)=tdmab(k)-2.d0/hhz(klev)
!mpk2d!mpk2d              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
!mpk2d!mpk2d     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
!mpk2d!mpk2d     +               /rval-(5.d0*psif(i,j-1,k)/3.d0
!mpk2d!mpk2d     +               -0.2*psif(i,j-2,k))/hhz(klev)+rf(i,j,k)
!mpk2d!mpk2d            else
!mpk2d!mpk2d             tdmab(k)=tdmab(k)-2.d0/hhz(klev)
!mpk2d!mpk2d              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))/hhr(klev)
!mpk2d!mpk2d     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))/hr(klev)
!mpk2d!mpk2d     +               /rval-(psif(i,j+1,k)+psif(i,j-1,k))/hhz(klev)
!mpk2d!mpk2d     +               +rf(i,j,k)
!mpk2d!mpk2d            endif
!mpk2d              tdmab(k)=tdmab(k)
!mpk2d              tdmar(k)=-(psif(i+1,j,k)+psif(i-1,j,k))*hhri
!mpk2d     +               -1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))*hri
!mpk2d     +               *rvali+rf(i,j,k)
!mpk2d   17     continue
!mpk2d          call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n3(klev))
!mpk2d          do 18 k=0,n3(klev)
!mpk2d            if(k.eq.0)then
!mpk2d              psif(i,j,k)=tdmas(n3(klev))
!mpk2d            else
!mpk2d              psif(i,j,k)=tdmas(k)
!mpk2d            endif
!mpk2d   18     continue
!mpk2d   16   continue
c
c                                 Calculating residual.
c
          resmax=0.d0
          do 20 k=1,n3(klev)
          do 20 i=1,n1(klev)-1
c8aug06            rval=dble(i)*hr(klev)+r0
c8aug06            rvali=1./rval
!mpk2d          do 20 j=1,n2(klev)
           j=1
c8aug06          do 20 k=1,n3(klev)
          ka=k+1
          km=k-1
          if(k.eq.n3(klev))ka=1
          temp=(psif(i+1,j,k)+psif(i-1,j,k))*hhri
     +        +1.5d0*(psif(i+1,j,k)-psif(i-1,j,k))*hri
     +        *rvali(i)+psif(i,j,k)*rvalisq(i)
     +        +(psif(i,j,ka)+psif(i,j,km))*hhpi*rvalisq(i)
     +        -psif(i,j,k)*(2.d0*hhri+2.d0*hhpi*rvalisq(i))
!mpk2d          if(j.eq.1)then
!mpk2d            temp=temp+(-2.d0*psif(i,j,k)+5.d0*psif(i,j+1,k)/3.d0
!mpk2d     +          -0.2d0*psif(i,j+2,k))/hhz(klev)
!mpk2d          elseif(j.eq.n2(klev))then
!mpk2d            temp=temp+(-2.d0*psif(i,j,k)+5.d0*psif(i,j-1,k)/3.d0
!mpk2d     +          -0.2d0*psif(i,j-2,k))/hhz(klev)
!mpk2d          else
!mpk2d            temp=temp+(psif(i,j+1,k)-2.d0*psif(i,j,k)+psif(i,j-1,k))
!mpk2d     +          /hhz(klev)
!mpk2d          endif
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
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,1,0:kp),psif(0:kr,1,0:kp)
      real*8 rc(0:kr,1,0:kp),rf(0:kr,1,0:kp),tempa(0:kr,1,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
!$OMP THREADPRIVATE(/backup/,/arrays/)

c
c------------------Start off with restriction
c
      call trans3(rc,rf,n1(klev),n2(klev),n3(klev))
c
      do 10 k=1,n3(klev+1)-1
      do 10 i=1,n1(klev+1)-1
!mpk2d      do 10 j=1,n2(klev+1)-1
       j=1
c8aug06      do 10 k=1,n3(klev+1)-1
        ii=i+i
!mpk2d     jj=j+j
        jj=1
        kk=k+k
        iip=ii+1
        iim=ii-1
!mpk2d        jjp=jj+1
        jjp=1
!mpk2d        jjm=jj-1
        jjm=1
        kkp=kk+1
        kkm=kk-1
!mpk2d a bug in original code? not used anyway        if(k.eq.0)jjm=n3(klev)-1
!mpk2d        if(k.eq.n3(klev+1))kkp=1
c
        psic(i,j,k)=0.d0
c13aug06 beg
        psif(i,j,k)=0.d0
c13aug06 end
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
!mpk2d      do 20 j=0,n2(klev+1)
      j=1
      do 20 k=0,n3(klev+1)
        rc(0,j,k)=0.d0
        rc(n1(klev+1),j,k)=0.d0
        psic(0,j,k)=0.d0
        psic(n1(klev+1),j,k)=0.d0
c13aug06 beg 
        psif(0,j,k)=0.d0
        psif(n1(klev+1),j,k)=0.d0
c13aug06 end

   20 continue
c
!mpk2d      do 30 i=0,n1(klev+1)
!mpk2d      do 30 k=0,n3(klev+1)
!mpk2d        rc(i,0,k)=0.d0
!mpk2d        rc(i,n2(klev+1),k)=0.d0
!mpk2d        psic(i,0,k)=0.d0
!mpk2d        psic(i,n2(klev+1),k)=0.d0
!mpk2d   30 continue
c
c      do 40 i=0,n1(klev+1)
c     do 40 j=0,n2(klev+1)
c       rc(i,j,0)=0.d0
c       rc(i,j,n3(klev+1))=0.d0
c       psic(i,j,0)=0.d0
c       psic(i,j,n3(klev+1))=0.d0
c  40 continue
c14aug06 discovery of a bug?
       do i=0,n1(klev+1)
        rc(i,1,0)=0.d0
        rc(i,1,n3(klev+1))=0.d0
        psic(i,1,0)=0.d0
        psic(i,1,n3(klev+1))=0.d0 
        psif(i,1,0)=0.d0
        psif(i,1,n3(klev+1))=0.d0 
       enddo
c14aug06 
c
      call trans1(rc,f,ibeg(klev+1),n1(klev+1),n2(klev+1),n3(klev+1))
c14aug06     call trans1(psic,v,ibeg(klev+1),n1(klev+1),n2(klev+1),n3(klev+1))
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
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,1,0:kp),psif(0:kr,1,0:kp)
      real*8 rc(0:kr,1,0:kp),rf(0:kr,1,0:kp),tempa(0:kr,1,0:kp)
      real*8 hr(6),hz(6),hp(6),hhr(6),hhz(6),hhp(6)
      integer ibeg(6),n1(6),n2(6),n3(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hz,hp,hhr,hhz,hhp
      common /ints/ibeg,n1,n2,n3
!$OMP THREADPRIVATE(/backup/,/arrays/)

c
      call trans2(v,psic,ibeg(klev),n1(klev),n2(klev),n3(klev))
c
c--------------------prolongation
c
      do 10 i=0,n1(klev-1)
!mpk2d      do 10 j=0,n2(klev-1)
      j=1
      do 10 k=0,n3(klev-1)
        i1=i/2
!mpk2d        j1=j/2
        j1=1
        k1=k/2
        i2=i1*2
!mpk2d        j2=j1*2
        k2=k1*2
        if(i2.eq.i)then
!mpk2d          if(j2.eq.j)then
            if(k2.eq.k)then
              psif(i,j,k)=psic(i1,j1,k1)
            else
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1,j1,k1+1))*0.5d0
            end if
!mpk2d          else
!mpk2d            if(k2.eq.k)then
!mpk2d              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1,j1+1,k1))*0.5d0
!mpk2d            else
!mpk2d              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1,j1+1,k1)
!mpk2d     +                    +psic(i1,j1,k1+1)+psic(i1,j1+1,k1+1))*0.25d0
!mpk2d            end if
!mpk2d          end if
        else
!mpk2d           if(j2.eq.j)then
            if(k2.eq.k)then
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1+1,j1,k1))*0.5d0
            else
              psif(i,j,k)=(psic(i1,j1,k1)+psic(i1+1,j1,k1)
     +                    +psic(i1,j1,k1+1)+psic(i1+1,j1,k1+1))*0.25d0
            end if
!mpk2d           else
!mpk2d             if(k2.eq.k)then
!mpk2d               psif(i,j,k)=(psic(i1,j1,k1)+psic(i1+1,j1,k1)
!mpk2d      +                    +psic(i1,j1+1,k1)+psic(i1+1,j1+1,k1))*0.25d0
!mpk2d             else
!mpk2d               psif(i,j,k)=(psic(i1,j1,k1)+psic(i1+1,j1,k1)
!mpk2d      +                     +psic(i1,j1+1,k1)+psic(i1+1,j1+1,k1)
!mpk2d      +                     +psic(i1,j1,k1+1)+psic(i1+1,j1,k1+1)
!mpk2d      +               +psic(i1,j1+1,k1+1)+psic(i1+1,j1+1,k1+1))*0.125d0
!mpk2d             end if
!mpk2d           end if
        end if
   10 continue
c
      call trans2(v,tempa,ibeg(klev-1),n1(klev-1),n2(klev-1),n3(klev-1))
c
      do 20 k=0,n3(klev-1)
      do 20 i=0,n1(klev-1)
!mpk2d      do 20 j=0,n2(klev-1)
      j=1
c8aug06      do 20 k=0,n3(klev-1)
        psif(i,j,k)=psif(i,j,k)+tempa(i,j,k)
        psic(i,j,k)=0.d0
   20 continue
c13aug06      call trans1(psic,v,ibeg(klev),n1(klev),n2(klev),n3(klev))
c14aug06      call trans1(psif,v,ibeg(klev-1),n1(klev-1),n2(klev-1),n3(klev-1))
      return
      end
c==================================================
      subroutine trans1(a3d,a1d,i1,n1,n2,n3)
      implicit real*8(a-h,o-z)
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 a3d(0:kr,1,0:kp),a1d(mgen)
c
      iimax=0
      do 10 k=0,n3
!mpk2d      do 10 j=0,n2
       j=1
      do 10 i=0,n1
!mpk2d        ii=k*(n1+1)*(n2+1) + j*(n1+1) + i + i1
             ii=k*(n1+1) + i + i1
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
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 a3d(0:kr,1,0:kp),a1d(mgen)
c
      do 10 k=0,n3
!mpk2d       do 10 j=0,n2
      j=1
      do 10 i=0,n1
!mpk2d         ii=k*(n1+1)*(n2+1) + j*(n1+1) + i + i1
         ii=k*(n1+1) + i + i1
        a3d(i,j,k)=a1d(ii)
   10 continue
      return
      end
c==================================================
      subroutine trans3(a3d1,a3d2,n1,n2,n3)
      implicit real*8(a-h,o-z)
      parameter(kr=290,kz=33,kp=1610,mgen=650000)
      real*8 a3d1(0:kr,1,0:kp),a3d2(0:kr,1,0:kp)
c
      do 10 k=0,n3
!mpk2d      do 10 j=0,n2
      j=1
      do 10 i=0,n1
        a3d2(i,j,k)=a3d1(i,j,k)
   10 continue
      return
      end
c==================================================
c      subroutine rufcon(a,nx,ny,nz)
c      parameter(kx=281,ky=33,kz=1210,mgen=500000)
c      implicit real*8(a-h,o-z)
c      real*8 a(0:kx,1,0:kz)
c      character*1 c(0:401)
c      character*1 ch(11)
c      ch(1)='0'
c      ch(2)='1'
c      ch(3)='2'
c      ch(4)='3'
c      ch(5)='4'
c      ch(6)='5'
c      ch(7)='6'
c      ch(8)='7'
c      ch(9)='8'
c      ch(10)='9'
c      ch(11)='t'
c
c--------Z-section at Z=Zmax/2
c
c      amax=biggest2(a,nx,ny,nz,ny/2)
c      amin=smallest2(a,nx,ny,nz,ny/2)
!mpk2d      print*,' '
!mpk2d      print*,'At z=zmax/2'
!mpk2d      print*,'arraymax,min=',amax,amin
!mpk2d      print*,'============'
!mpk2d      print*,'Abscissa = r-axis, Ordinate = phi-axis'
c
c      adiff=amax-amin
c      if(adiff.ne.0.)adiff=1./adiff
c
c      do 20 k=nz,0,-1
c        do 10 i=0,nx
!mpk2d          x=1.5+10.*(a(i,ny/2,k)-amin)*adiff
c          x=1.5+10.*(a(i,1,k)-amin)*adiff
c          n=x
c          c(i)=ch(n)
c  10    continue
c        write(6,100)(c(i),i=0,nx)
c  20  continue
c
c      print*,' '
c      print*,'arraymax,min=',amax,amin
c
c--------Z-section at r=rmax/2
c
!mpk2d      amax=biggest1(a,nx,ny,nz,nx/2)
!mpk2d      amin=smallest1(a,nx,ny,nz,nx/2)
!mpk2d      print*,' '
!mpk2d      print*,'At r=rmax/2'
!mpk2d      print*,'arraymax,min=',amax,amin
!mpk2d      print*,'============'
!mpk2d      print*,'Abscissa = z-axis, Ordinate = phi-axis'
c
!mpk2d      adiff=amax-amin
!mpk2d      if(adiff.ne.0.)adiff=1./adiff
c
!mpk2d      do 30 k=nz,0,-1
!mpk2d        do 40 j=0,ny
!mpk2d          x=1.5+10.*(a(nx/2,j,k)-amin)*adiff
!mpk2d          n=x
!mpk2d          c(j)=ch(n)
!mpk2d  40    continue
!mpk2d        write(6,100)(c(j),j=0,ny)
!mpk2d  30  continue
c
c
c--------Z-section at phi=phimax/2
c
!mpk2d      amax=biggest3(a,nx,ny,nz,nz/2)
!mpk2d      amin=smallest3(a,nx,ny,nz,nz/2)
!mpk2d      print*,' '
!mpk2d      print*,'At phi=phimax/2'
!mpk2d      print*,'arraymax,min=',amax,amin
!mpk2d      print*,'============'
!mpk2d      print*,'Abscissa = z-axis, Ordinate = r-axis'
c
!mpk2d      adiff=amax-amin
!mpk2d      if(adiff.ne.0.)adiff=1./adiff
c
!mpk2d      do 50 i=nx,0,-1
!mpk2d        do 60 j=0,ny
!mpk2d          x=1.5+10.*(a(i,j,nz/2)-amin)*adiff
!mpk2d          n=x
!mpk2d          c(j)=ch(n)
!mpk2d  60    continue
!mpk2d        write(6,100)(c(j),j=0,ny)
!mpk2d  50  continue
c
c 100  format(1x,200a1)
c      return
c      end
c==================================================
c      double precision function smallest1(a,nx,ny,nz,ii)
c      parameter(kx=281,ky=33,kz=1210,mgen=202590)
c      implicit real*8(a-h,o-z)
c      real*8 a(0:kx,1,0:kz)
c
!mpk2d      z=a(ii,0,0)
c      z=a(ii,1,0)
!mpk2d      do 10 j=0,ny
c       j=1
c      do 10 k=0,nz
c        if(z.gt.a(ii,j,k))z=a(ii,j,k)
c  10  continue
c      smallest1=z
c      return
c      end
c==================================================
c      double precision function biggest1(a,nx,ny,nz,ii)
c      parameter(kx=281,ky=33,kz=1210,mgen=202590)
c      implicit real*8(a-h,o-z)
c      real*8 a(0:kx,1,0:kz)
c
c      z=a(ii,0,0)
!mpk2d      do 10 j=0,ny
c      j=1
c      do 10 k=0,nz
c        if(z.lt.a(ii,j,k))z=a(ii,j,k)
c  10  continue
c      biggest1=z
c      return
c      end
c==================================================
c      double precision function smallest2(a,nx,ny,nz,jj)
c      parameter(kx=281,ky=33,kz=1210,mgen=202590)
c      implicit real*8(a-h,o-z)
c      real*8 a(0:kx,1,0:kz)
c
c      jj=1
c      z=a(0,jj,0)
c      do 10 i=0,nx
c      do 10 k=0,nz
c        if(z.gt.a(i,jj,k))z=a(i,jj,k)
c  10  continue
c      smallest2=z
c      return
c      end
c==================================================
c      double precision function biggest2(a,nx,ny,nz,jj)
c      parameter(kx=281,ky=33,kz=1210,mgen=202590)
c      implicit real*8(a-h,o-z)
c      real*8 a(0:kx,1,0:kz)
c
c      jj=1
c      z=a(0,jj,0)
c      do 10 i=0,nx
c      do 10 k=0,nz
c        if(z.lt.a(i,jj,k))z=a(i,jj,k)
c  10  continue
c      biggest2=z
c      return
c      end
c==================================================
c      double precision function smallest3(a,nx,ny,nz,kk)
c      parameter(kx=281,ky=33,kz=1210,mgen=202590)
c      implicit real*8(a-h,o-z)
c      real*8 a(0:kx,1,0:kz)
c
!mpk2d      z=a(0,0,kk)
c      zz=a(0,1,kk)
c      do 10 i=0,nx
!mpk2d      do 10 j=0,ny
c        j=1
c        if(z.gt.a(i,j,kk))z=a(i,j,kk)
c  10  continue
c      smallest3=z
c      return
c      end
c==================================================
c      double precision function biggest3(a,nx,ny,nz,kk)
c      parameter(kx=281,ky=33,kz=1210,mgen=202590)
c      implicit real*8(a-h,o-z)
c      real*8 a(0:kx,1,0:kz)
c
c      z=a(0,0,kk)
c      do 10 i=0,nx
!mpk2d      do 10 j=0,ny
c       j=1
c        if(z.lt.a(i,j,kk))z=a(i,j,kk)
c  10  continue
c      biggest3=z
c      return
c      end
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

