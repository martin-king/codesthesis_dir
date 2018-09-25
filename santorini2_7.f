c2345678901234567890123456789012345678901234567890123456789012345678901c
c     Program santorini2.f                                             c
c     ==============                                                   c
c     T.Lewis's rot7.f modified by martin king, 1 March 1999.          c                                               
c     Mathematical modelling of flow in a rotating annulus             c
c     in r-z plane.                                                    c
c     with heat flux applied in the axial direction(mpk,               c
c     boussinesq off, ie isothermal incompressible).                   c
c     Second order computation on a non-uniform grid.                  c
c----------------------------------------------------------------------c
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,mgen=100000,na=8)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax)!mpk,thetan(0:nrmax,0:nzmax)
!mpk      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      real*8 omegmax(0:mgen),omegmin(0:mgen)
!mpk      real*8 thetmax(0:mgen),thetmin(0:mgen)
      real*8 phimax(0:mgen),phimin(0:mgen),czz(nzmax,4,na)
      real*8 cr(nrmax,3,na),crr(nrmax,4,na),cz(nzmax,3,na)
      real*8 hr(nrmax,na),hz(nzmax,na),r(0:nrmax,na),z(0:nzmax,na)
      real*8 mean_vel(nrmax/2,3),mean_velsqr(nrmax/2,3),s(nrmax/2,3)
      real*8 mean_uel(nzmax/2,2),mean_uelsqr(nzmax/2,2),suel(nzmax/2,2)
      integer n1(na),n2(na),ibeg(na)
      character*20 fname1,fname2,fname3,fname4,fname5
      character*20 fname6,fname8,fname10,fname,ans
      common/val/Re,Ptl !,Ra
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
!mpk      common/sol2/thetan,theta,thetap
      common/num/hr,hz,n1,n2,dt,r,z
      common/coeff/cr,crr,cz,czz
      common/ints/ibeg
c
      print*,'grid and coefficient data filename: '
      read*,fname8
      print*,'filename for vorticity at the end of computation: '
      read*,fname1
      print*,'filename for stream function at the end of computation: '
      read*,fname2
      print*,'filename for tangential velocity at the end of ' 
      print*,'computation: '
      read*,fname3
!mpk      print*,'filename for temperature: '
!mpk      read*,fname4
      print*,'iterations between wrt file for time-averaged tangential ' 
      print*,'velocities & standard deviations?(0 for no): '
      read*,nff
      print*,'filename for max/min vorticity data' 
      print*,'at every timestep(filename or "no"): '
      read*,fname5
      print*,'filename for max/min stream function ' 
      print*,'at every timestep(filename or "no"): ' 
      read*,fname6
      print*,'time step: '
      read*,dt
      print*,'maximum number of time steps: '
      read*,nts
      print*,'iterations between writing to file(stream function): '
      read*,nf
      print*,'iterations to be completed before writing(stream '
      print*,'function): '
      read*,istart
      print*,'rotational Reynolds number: '
      read*,Re
      print*,'Prandtl number: '
      read*,Ptl
!mpk      print*,'T_h and T_c: '
!mpk      read*,th,tc
!mpk      if(th.lt.tc)stop
c
!mpk      Ra=(Re**2)*Ptl*(th-tc)/((th+tc)*0.5d0+273.d0)
!mpk      print*,Ra
          Ra=0.d0                   !to rm later. just in case.
c
c----------------------------------------------------------------------c
c     Read grid data file                                              c
c----------------------------------------------------------------------c
c
      open(7,file=fname8)
      read(7,*)rmax,r0,zmax,zmin
      read(7,*)nlev
      read(7,*)(ibeg(i),i=1,nlev)
      do 5 l=1,nlev
        read(7,*)n1(l),n2(l)
        read(7,*)(hr(i,l),i=1,n1(l))
        read(7,*)(hz(j,l),j=1,n2(l))
        read(7,*)((cr(i,k,l),k=1,3),i=1,n1(l)/2)
        read(7,*)((crr(i,k,l),k=1,4),i=1,n1(l)/2)
        read(7,*)((cr(i,k,l),k=1,3),i=n1(l)/2+1,n1(l)-1)
        read(7,*)((crr(i,k,l),k=1,4),i=n1(l)/2+1,n1(l)-1)
        read(7,*)((cz(j,k,l),k=1,3),j=1,n2(l)/2)
        read(7,*)((czz(j,k,l),k=1,4),j=1,n2(l)/2)
        read(7,*)((cz(j,k,l),k=1,3),j=n2(l)/2+1,n2(l)-1)
        read(7,*)((czz(j,k,l),k=1,4),j=n2(l)/2+1,n2(l)-1)
    5 continue
      close(7)
c
c----------------------------------------------------------------------c
c     Calculate constants and r,z values                               c
c----------------------------------------------------------------------c
c
      pi=dacos(-1.d0)
c
      do 3 l=1,nlev
        r(0,l)=r0
        z(0,l)=zmin
        do 1 i=1,n1(l)
          r(i,l)=r(i-1,l)+hr(i,l)
    1   continue
        do 2 j=1,n2(l)
          z(j,l)=z(j-1,l)+hz(j,l)
    2   continue
    3 continue
c
c----------------------------------------------------------------------c
c     Set initial conditions                                           c
c----------------------------------------------------------------------c
c
      do 10 i=0,n1(1)
      do 10 j=0,n2(1)
        omeg(i,j)=0.d0
        phi(i,j)=0.d0
c        vel(i,j)=-Ptl*Re*r(i,1)+0.01*Re*z(j,1) !var in z-dir
c        vel(n1(1),j)=-Ptl*Re*r(i,1)
        vel(i,j)=-Ptl*Re*r(i,1)   ! start with 0 rotation(rot frame)
        vel(i,0)=0.d0             !left disk
        vel(i,n2(1))=0.d0         !right disk
        vel(0,j)=0.d0             !hub
!        vel(n1(1),j)=0.d0         !shroud. this is a solid body rot.
!mpk  see subroutine update(line 660) for BC's
!mpk        theta(i,j)=0.d0
!mpk        theta(i,0)=1.d0
  10  continue
c
      print*,'initial conditions from file (yes/no):'
      read*,ans
      if(ans .ne. 'no')then
        print*,'filename for vorticity: '
        read*,fname10
        open(9,file=fname10,status='unknown')
        read(9,*)n1(1),n2(1)
        read(9,*)rmax,r0,zmax,zmin
        read(9,*)(r(i,1),i=0,n1(1))
        read(9,*)(z(j,1),j=0,n2(1))
        read(9,*)((omeg(i,j),i=0,n1(1)),j=0,n2(1))
        close(9)
        print*,'filename for strm function: '
        read*,fname10
        open(9,file=fname10,status='unknown')
        read(9,*)n1(1),n2(1)
        read(9,*)rmax,r0,zmax,zmin
        read(9,*)(r(i,1),i=0,n1(1))
        read(9,*)(z(j,1),j=0,n2(1))
        read(9,*)((phi(i,j),i=0,n1(1)),j=0,n2(1))
        close(9)
        print*,'filename for velocity: '
        read*,fname10
        open(9,file=fname10,status='unknown')
        read(9,*)n1(1),n2(1)
        read(9,*)rmax,r0,zmax,zmin
        read(9,*)(r(i,1),i=0,n1(1))
        read(9,*)(z(j,1),j=0,n2(1))
        read(9,*)((vel(i,j),i=0,n1(1)),j=0,n2(1))
        close(9)
!mpk        print*,'filename for temperature: '
!mpk        read*,fname10
!mpk        open(9,file=fname10,status='unknown')
!mpk        read(9,*)n1(1),n2(1)
!mpk        read(9,*)rmax,r0,zmax,zmin
!mpk        read(9,*)(r(i,1),i=0,n1(1))
!mpk        read(9,*)(z(j,1),j=0,n2(1))
!mpk        read(9,*)((theta(i,j),i=0,n1(1)),j=0,n2(1))
!mpk        close(9)
        call rufcon2(phi,omeg,n1(1),n2(1))!,amax,amin)
      endif
c
      do 50 i=0,n1(1)
      do 50 j=0,n2(1)
        omegp(i,j)=omeg(i,j)
        velp(i,j)=vel(i,j)
!mpk        thetap(i,j)=theta(i,j)
  50  continue
c
c----------------------------------------------------------------------c
c     Solve for vorticity, stream function, axial velocity and         c
c     (!not)temperature on internal points. buoyancy correction also   c
c     switched off.!mpk
c----------------------------------------------------------------------c
c
      ap=0
      do 20 k=1,nts
        print*,'---------------------------------'
        call vorticity
        call velocity
!mpk        call temperature
        call poisson(phi,omegn,nlev,vel)
        call update
c        amax=woppa(vel,n1(1),n2(1))
c        amin=titch(vel,n1(1),n2(1))
c        print*,amax,amin
!mpk        if(k/10.EQ.dreal(k)/10)then
!mpk        call rufcon2(phi,omeg,n1(1),n2(1))!,amax,amin)
!mpk        endif
c----------------------------------------------------------------------c
c     Convergence and others check                                     c
c----------------------------------------------------------------------c
c
        errmax1=0.d0
        errmax2=0.d0
        errmax3=0.d0
        gmax=-10000.d0
        gmin=10000.d0
        pimax=-10000.d0
        pimin=10000.d0
        thmax=-10000.d0
        thmin=10000.d0
c
        do 30 i=0,n1(1)
         do 30 j=0,n2(1)
          error1=dabs(omeg(i,j)-omegp(i,j))
          error2=dabs(vel(i,j)-velp(i,j))
!mpk          error3=dabs(theta(i,j)-thetap(i,j))
          if(error1 .gt. errmax1)errmax1=error1
          if(error2 .gt. errmax2)errmax2=error2
!mpk          if(error3 .gt. errmax3)errmax3=error3
          if(omeg(i,j) .gt. gmax)gmax=omeg(i,j)
          if(omeg(i,j) .lt. gmin)gmin=omeg(i,j)
          if(phi(i,j) .gt. pimax)pimax=phi(i,j)
          if(phi(i,j) .lt. pimin)pimin=phi(i,j)
!mpk          if(theta(i,j) .gt. thmax)thmax=theta(i,j)
!mpk          if(theta(i,j) .lt. thmin)thmin=theta(i,j)
!mpk  borrow thmax, thmin for checking solid body rotation criteria.          
c          if(ABS(vel(i,j)/(Ptl*Re*r(i,1))) .gt. thmax)thmax=ABS(
c&          vel(i,j)/(Ptl*Re*r(i,1)))
c          if(ABS(vel(i,j)/(Ptl*Re*r(i,1))) .lt. thmin)thmin=ABS(
c&          vel(i,j)/(Ptl*Re*r(i,1)))
  30    continue
c
        omegmax(k)=gmax
        omegmin(k)=gmin
        phimax(k)=pimax
        phimin(k)=pimin
!mpk        thetmax(k)=thmax
!mpk        thetmin(k)=thmin
c
        print*,errmax1,errmax2,k
c       print*,pimax,pimin,k
c
c      WRITE(10,33)k*dt,thmin,thmax
c  33  FORMAT(3(1x,e20.10E2))
        if(errmax1.lt.1d-6.and.errmax2.lt.1d-6)then
          print*,'converged'
          goto 40
        endif
!mpk    to calculate statistical infor at z/s=0.2, 0.5 & 0.8
        if(nff .ne. 0)then
         do 34 n=1,3
         do m1=1,n2(1)
          if(z(m1,1)/z(n2(1),1).gt.(0.2d0+(n-1)*0.3d0))then 
            j_z=m1
            GOTO 33
          endif
         enddo
c         j_z=n2(1)*(0.2d0+(n-1)*0.3d0)  
  33     do 35 m=0,n1(1)/2
           mean_vel(m,n)=(mean_vel(m,n)*(k-1)+vel(2*m,j_z))/k
           mean_velsqr(m,n)=(mean_velsqr(m,n)*(k-1)+vel(2*m,j_z)
&          **2)/k
           s(m,n)=(mean_velsqr(m,n)-(mean_vel(m,n)**2.d0))**0.5d0
  35  continue
  34  continue
        endif
!mpk    writing statistical at nx_3 & nx_5,at z/s=0.8,at evry tmstp
        nx_5=n1(1)/2*0.5
        nx_3=n1(1)/2*0.3
        PRINT*,mean_vel(nx_3,3),mean_vel(nx_5,3)
        WRITE(99,51)k*dt,mean_vel(nx_3,3),mean_vel(nx_5,3)
  51    FORMAT(3(1x,e15.5))
!mpk
        if(nff .ne. 0)then
        do 37 n=1,2
        do m1=1,n1(1)
         if(r(m1,1)/rmax.gt.(0.75d0+(n-1)*0.1d0))then
           j_r=m1
           GOTO 36
          endif
         enddo  
  36      do 38 m=0,n2(1)/2
           uel=dum8(phi,cz,j_r,2*m)
           mean_uel(m,n)=(mean_uel(m,n)*(k-1)+uel)/k
           mean_uelsqr(m,n)=(mean_uelsqr(m,n)*(k-1)+uel**2)/k
           suel(m,n)=(mean_uelsqr(m,n)-(mean_uel(m,n)**2.d0))**0.5d0
  38  continue
  37  continue
        endif             
        if(k.gt.istart)then
          if(k/nf.eq.dble(k)/dble(nf))then
            print*,'writing to file'
            numb=k/nf-istart/nf
c
            write(fname,400)'strm3',numb
  400       format(A5,I3.3)
            open(unit=400,file=fname)
            write(400,*)n1(1),n2(1)
            write(400,*)rmax,r0,zmax,0
            write(400,*)((phi(i,j),i=0,n1(1)),j=0,n2(1))
            close(400)
          endif
        endif
c        
!mpk wrtng to fl evry nff tmstp        
        if(nff .ne. 0)then
         if(k/nff .eq. dble(k)/dble(nff))then
          numb=k/nff 
c         
          write(fname,405)'stats',numb
  405     format(A5,I3.3)
          open(unit=405,file=fname)
          write(405,*),'        means','        s','        mean_velsqr'
          write(405,55)((in*1.d0,r(2*m,1),mean_vel(m,in),s(m,in),
&         mean_velsqr(m,in),m=0,n1(1)/2),in=1,3)
          write(405,55)((in*1.d0,z(2*m,1),mean_uel(m,in),suel(m,in),
&         mean_uelsqr(m,in),m=0,n2(1)/2),in=1,2)          
          close(405)
         endif
        endif
  55    format(5(1x,e20.10E4))
!mpk
      call cm(rmax,r0,zmax,zmin,k,nts,ap)
  20  continue
      WRITE(99,*)r(nx_3*2,1),r(nx_5*2,1)
c
      print*,'not converged---steady state criteria not reached'
c
  40  continue
c
c---------------------------------------------------------------------c
c     Write vorticity, stream function, velocity and (!not)           c
c     temperature data to files                                       c                                                         c
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
        write(13,*)errmax1,errmax2,nts
        close(13)
      endif
c
      open(7,file=fname1,status='unknown')
      write(7,*)n1(1),n2(1)
      write(7,*)rmax,r0,zmax,0
      write(7,*)(r(i,1),i=0,n1(1))
      write(7,*)(z(j,1),j=0,n2(1))
      write(7,*)((omeg(i,j),i=0,n1(1)),j=0,n2(1))
      close(7)
c
      open(8,file=fname2,status='unknown')
      write(8,*)n1(1),n2(1)
      write(8,*)rmax,r0,zmax,0
      write(8,*)(r(i,1),i=0,n1(1))
      write(8,*)(z(j,1),j=0,n2(1))
      write(8,*)((phi(i,j),i=0,n1(1)),j=0,n2(1))
      close(8)
c
      open(9,file=fname3,status='unknown')
      write(9,*)n1(1),n2(1)
      write(9,*)rmax,r0,zmax,0
      write(9,*)(r(i,1),i=0,n1(1))
      write(9,*)(z(j,1),j=0,n2(1))
      write(9,*)((vel(i,j),i=0,n1(1)),j=0,n2(1))
      close(9)
c
!      open(10,file=fname4,status='unknown')
!      write(10,*)n1(1),n2(1)
!      write(10,*)rmax,r0,zmax,0
!      write(10,*)(r(i,1),i=0,n1(1))
!      write(10,*)(z(j,1),j=0,n2(1))
!      write(10,*)((theta(i,j),i=0,n1(1)),j=0,n2(1))
!      close(10)
c
  100  format(3(1x,e16.8))
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
      parameter(nrmax=250,nzmax=250,mgen=100000,na=8)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax) !,thetan(0:nrmax,0:nzmax)
!      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      real*8 czz(nzmax,4,na)
      real*8 cr(nrmax,3,na),crr(nrmax,4,na),cz(nzmax,3,na)
      real*8 hr(nrmax,na),hz(nzmax,na),r(0:nrmax,na),z(0:nzmax,na)
      integer n1(na),n2(na)
      common/val/Re,Ptl !,Ra
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
!      common/sol2/thetan,theta,thetap
      common/num/hr,hz,n1,n2,dt,r,z
      common/coeff/cr,crr,cz,czz
c
      do 10 i=1,n1(1)-1
      do 10 j=1,n2(1)-1
        if(i.le.n1(1)/2)then
          omegrra=dum1a(omeg,crr,i,j)
          omegra=dum2a(omeg,cr,i,j)
          omegr=dum2(omeg,cr,i,j)        
          phir=dum2(phi,cr,i,j)
          cnr=crr(i,2,1)/2.d0+cr(i,2,1)/(2.d0*r(i,1))
c explicit          cnr=crr(i,2,1)+cr(i,2,1)/r(i,1)
        else
          omegrra=dum3a(omeg,crr,i,j)
          omegra=dum4a(omeg,cr,i,j)
          omegr=dum4(omeg,cr,i,j)
          phir=dum4(phi,cr,i,j)
          cnr=crr(i,3,1)/2.d0+cr(i,2,1)/(2.d0*r(i,1))
c explicit          cnr=crr(i,3,1)+cr(i,2,1)/r(i,1)
        endif
c	
        if(j.le.n2(1)/2)then
          omegzza=dum5a(omeg,czz,i,j)
          omegza=dum6a(omeg,cz,i,j)
          omegz=dum6(omeg,cz,i,j)
          phiz=dum6(phi,cz,i,j)
          velz=dum6(vel,cz,i,j)
!mpk          thetaz=dum6(theta,cz,i,j)
          cnz=czz(j,2,1)/2.d0
        else
          omegzza=dum7a(omeg,czz,i,j) 
          omegza=dum8a(omeg,cz,i,j) 
          omegz=dum8(omeg,cz,i,j)
          phiz=dum8(phi,cz,i,j) 
          velz=dum8(vel,cz,i,j) 
!mpk          thetaz=dum8(theta,cz,i,j)
          cnz=czz(j,3,1)/2.d0
        endif
c
        Pe=dabs(0.5d0*phiz*(hr(i,1)+hr(i+1,1))/Ptl)
        if(Pe.gt.2)then
        if(phiz.gt.0)then
          omegr=(omeg(i,j)-omeg(i-1,j))/hr(i,1)
          phir=(phi(i,j)-phi(i-1,j))/hr(i,1)
        elseif(phiz.lt.0)then
          omegr=(omeg(i+1,j)-omeg(i,j))/hr(i+1,1)
          phir=(phi(i+1,j)-phi(i,j))/hr(i+1,1)
        endif
        endif
c
        Pe=dabs(0.5d0*(-phir-phi(i,j)/r(i,1))*(hz(j,1)+hz(j+1,1))/Ptl)
        if(Pe.gt.2)then
        if(-phir-phi(i,j)/r(i,1).gt.0)then
          omegz=(omeg(i,j)-omeg(i,j-1))/hz(j,1)
          phiz=(phi(i,j)-phi(i,j-1))/hz(j,1)
        elseif(-phir-phi(i,j)/r(i,1).lt.0)then
          omegz=(omeg(i,j+1)-omeg(i,j))/hz(j+1,1)
          phiz=(phi(i,j+1)-phi(i,j))/hz(j+1,1)
        endif
        endif
c       
        RHS1=omegrra+omegra/r(i,1)+omegzza
     +     +omeg(i,j)/(dt*Ptl)
        ULT=((phiz*omegr-phir*omegz)
     +     -(phi(i,j)*omegz)/r(i,1))/Ptl
        RHS2=-2.d0*vel(i,j)*velz/(r(i,1)*Ptl)-2.d0*Re*velz
!mpk     +      +2.d0*Ra*(theta(i,j)*velz+thetaz*vel(i,j))/(Re*Ptl)
!mpk     +      +r(i,1)*Ra*thetaz
        cf1=1.d0/(Ptl*dt)+1.d0/(2.d0*r(i,1)**2)-cnr-cnz
     +     -phiz/(2.d0*r(i,1)*Ptl)
        cf2=-1.d0/(2.d0*r(i,1)**2)+cnr+cnz+phiz/(2.d0*r(i,1)*Ptl)
c
        omegn(i,j)=(cf2*omegp(i,j)-ULT-RHS2+RHS1)/cf1
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
      parameter(nrmax=250,nzmax=250,mgen=100000,na=8)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax) !mpk,thetan(0:nrmax,0:nzmax)
!mpk      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      real*8 czz(nzmax,4,na)
      real*8 cr(nrmax,3,na),crr(nrmax,4,na),cz(nzmax,3,na)
      real*8 hr(nrmax,na),hz(nzmax,na),r(0:nrmax,na),z(0:nzmax,na)
      integer n1(na),n2(na)
      common/val/Re,Ptl !,Ra
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
!mpk      common/sol2/thetan,theta,thetap
      common/num/hr,hz,n1,n2,dt,r,z
      common/coeff/cr,crr,cz,czz
c     
      do 10 i=1,n1(1)-1
      do 10 j=1,n2(1)-1
        if(i.le.n1(1)/2)then
          velrra=dum1a(vel,crr,i,j)
          velra=dum2a(vel,cr,i,j)
          velr=dum2(vel,cr,i,j)
          phir=dum2(phi,cr,i,j)
          cnr=crr(i,2,1)/2.d0+cr(i,2,1)/(2.d0*r(i,1))
c explicit          cnr=crr(i,2,1)+cr(i,2,1)/r(i,1)
        else
          velrra=dum3a(vel,crr,i,j)
          velra=dum4a(vel,cr,i,j)
          velr=dum4(vel,cr,i,j)
          phir=dum4(phi,cr,i,j)
          cnr=crr(i,3,1)/2.d0+cr(i,2,1)/(2.d0*r(i,1))
c explicit          cnr=crr(i,3,1)+cr(i,2,1)/r(i,1)
        endif
c
        if(j.le.n2(1)/2)then
          velzza=dum5a(vel,czz,i,j)
          velza=dum6a(vel,cz,i,j)
          velz=dum6(vel,cz,i,j)
          phiz=dum6(phi,cz,i,j)
          cnz=czz(j,2,1)/2.d0
        else
          velzza=dum7a(vel,czz,i,j)
          velza=dum8a(vel,cz,i,j)
          velz=dum8(vel,cz,i,j)
          phiz=dum8(phi,cz,i,j)
          cnz=czz(j,3,1)/2.d0
        endif
c
        Pe=dabs(0.5d0*phiz*(hr(i,1)+hr(i+1,1))/Ptl)
        if(Pe.gt.2)then
        if(phiz.gt.0)then
          velr=(vel(i,j)-vel(i-1,j))/hr(i,1)
          phir=(phi(i,j)-phi(i-1,j))/hr(i,1)
        elseif(phiz.lt.0)then
          velr=(vel(i+1,j)-vel(i,j))/hr(i+1,1)
          phir=(phi(i+1,j)-phi(i,j))/hr(i+1,1)
        endif
        endif
c
        Pe=dabs(0.5d0*(-phir-phi(i,j)/r(i,1))*(hz(j,1)+hz(j+1,1))/Ptl)
        if(Pe.gt.2)then
        if(-phir-phi(i,j)/r(i,1).gt.0)then
          velz=(vel(i,j)-vel(i,j-1))/hz(j,1)
          phiz=(phi(i,j)-phi(i,j-1))/hz(j,1)
        elseif(-phir-phi(i,j)/r(i,1).lt.0)then
          velz=(vel(i,j+1)-vel(i,j))/hz(j+1,1)
          phiz=(phi(i,j+1)-phi(i,j))/hz(j+1,1)
        endif
        endif
c
        RHS1=velrra+velra/r(i,1)+velzza+vel(i,j)/(Ptl*dt)
        ULT=((phiz*velr-phir*velz)
     +     -(phi(i,j)*velz)/r(i,1))/Ptl
        RHS2=2.d0*Re*phiz !mpk -2.d0*Ra*theta(i,j)*phiz/(Re*Ptl)
        cf1=1.d0/(Ptl*dt)+1.d0/(2.d0*r(i,1)**2)-cnr-cnz
     +     +phiz/(2.d0*r(i,1)*Ptl)
        cf2=-1.d0/(2.d0*r(i,1)**2)+cnr+cnz-phiz/(2.d0*r(i,1)*Ptl)
c
        veln(i,j)=(cf2*velp(i,j)-ULT-RHS2+RHS1)/cf1
c
  10  continue
c
        return
        end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     subroutine TEMPERATURE                                           c
c     ======================                                           c
c     Calculates temperature using the duFort-Frankel method           c
c----------------------------------------------------------------------c
c
!mpk      subroutine temperature
c
!mpk      implicit real*8(a-h,o-z)
!mpk      parameter(nrmax=128,nzmax=128,mgen=100000,na=8)
!mpk      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
!mpk      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
!mpk      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
!mpk      real*8 velp(0:nrmax,0:nzmax),thetan(0:nrmax,0:nzmax)
!mpk      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
!mpk      real*8 czz(nzmax,4,na)
!mpk      real*8 cr(nrmax,3,na),crr(nrmax,4,na),cz(nzmax,3,na)
!mpk      real*8 hr(nrmax,na),hz(nzmax,na),r(0:nrmax,na),z(0:nzmax,na)
!mpk      integer n1(na),n2(na)
!mpk      common/val/Ra,Re,Ptl
!      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
!      common/sol2/thetan,theta,thetap
!      common/num/hr,hz,n1,n2,dt,r,z
!      common/coeff/cr,crr,cz,czz
c
!      do 10 i=0,n1(1)
!      do 10 j=1,n2(1)-1
!        if(i.le.n1(1)/2.and.i.ne.0)then
!          thetarra=dum1a(theta,crr,i,j)
!          thetara=dum2a(theta,cr,i,j)
!          thetar=dum2(theta,cr,i,j)
!          phir=dum2(phi,cr,i,j)
!          cnr=crr(i,2,1)/2.d0+cr(i,2,1)/(2.d0*r(i,1))
!        elseif(i.eq.0)then
!          thetarra=2.d0*theta(i+1,j)/hr(1,1)**2
!          thetara=0.d0
!          thetar=0.d0
!          phir=0.d0
!          cnr=-1.d0/hr(1,1)**2
!        elseif(i.eq.n1(1))then
!          thetarra=2.d0*theta(i-1,j)/hr(i,1)**2
!          thetara=0.d0
!          thetar=0.d0
!          phir=0.d0
!          cnr=-1.d0/hr(i,1)**2
!        else
!          thetarra=dum3a(theta,crr,i,j)
!          thetara=dum4a(theta,cr,i,j)
!          thetar=dum4(theta,cr,i,j)
!          phir=dum4(phi,cr,i,j)
!          cnr=crr(i,3,1)/2.d0+cr(i,2,1)/(2.d0*r(i,1))
!        endif
c
!        if(j.le.n2(1)/2.and.j.ne.0)then
!          thetazza=dum5a(theta,czz,i,j)
!          thetaza=dum6a(theta,cz,i,j)
!          thetaz=dum6(theta,cz,i,j)
!          phiz=dum6(phi,cz,i,j)
!          cnz=czz(j,2,1)/2.d0
!        else
!          thetazza=dum7a(theta,czz,i,j)
!          thetaza=dum8a(theta,cz,i,j)
!          thetaz=dum8(theta,cz,i,j)
!          phiz=dum8(phi,cz,i,j)
!          cnz=czz(j,3,1)/2.d0
!        endif
c
!        Pe=dabs(0.5d0*phiz*(hr(i,1)+hr(i+1,1))/Ptl)
!        if(Pe.gt.2.d0.and.i.ge.1.and.i.le.n1(1)-1)then
!        if(phiz.gt.0)then
!          thetar=(theta(i,j)-theta(i-1,j))/hr(i,1)
!          phir=(phi(i,j)-phi(i-1,j))/hr(i,1)
!        elseif(phiz.lt.0)then
!          thetar=(theta(i+1,j)-theta(i,j))/hr(i+1,1)
!          phir=(phi(i+1,j)-phi(i,j))/hr(i+1,1)
!        endif
!        endif
c
!        Pe=dabs(0.5d0*(-phir-phi(i,j)/r(i,1))*(hz(j,1)+hz(j+1,1))/Ptl)
!        if(Pe.gt.2.d0)then
!        if(-phir-phi(i,j)/r(i,1).gt.0.d0)then
!          thetaz=(theta(i,j)-theta(i,j-1))/hz(j,1)
!          phiz=(phi(i,j)-phi(i,j-1))/hz(j,1)
!        elseif(-phir-phi(i,j)/r(i,1).lt.0.d0)then
!          thetaz=(theta(i,j+1)-theta(i,j))/hz(j+1,1)
!          phiz=(phi(i,j+1)-phi(i,j))/hz(j+1,1)
!        endif
!        endif
c
!        RHS1=thetarra+thetara/r(i,1)+thetazza+theta(i,j)/dt
!        ULT=phiz*thetar-phir*thetaz
!     +     -phi(i,j)*thetaz/r(i,1)
!        cf1=1.d0/(dt)-cnr-cnz
!        cf2=cnr+cnz
c
!        thetan(i,j)=(cf2*thetap(i,j)-ULT+RHS1)/cf1
c
!  10  continue
c
!      return
!mpk      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     subroutine UPDATE                                                c
c     =================                                                c
c     updates the values for vorticity, velocity and (not) temperature c
c     for the next iteration. Boundary values for the two parameters   c
c     are also calculated                                              c
c----------------------------------------------------------------------c
c
      subroutine update
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,mgen=100000,na=8)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax) !mpk ,thetan(0:nrmax,0:nzmax)
!mpk      real*8 theta(0:nrmax,0:nzmax),thetap(0:nrmax,0:nzmax)
      real*8 hr(nrmax,na),hz(nzmax,na),r(0:nrmax,na),z(0:nzmax,na)
      integer n1(na),n2(na)
      common/val/Re,Ptl !,Ra
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
!mpk      common/sol2/thetan,theta,thetap
      common/num/hr,hz,n1,n2,dt,r,z
c
      do 10 i=0,n1(1)
!        thetan(i,0)=1.d0
        omegn(i,0)=2.d0*phi(i,1)/hz(1,1)**2   !left disk
        veln(i,0)=0.d0                        !left disk
!        thetan(i,n2(1))=0.d0
        omegn(i,n2(1))=2.d0*phi(i,n2(1)-1)/hz(n2(1),1)**2    !right disk
        veln(i,n2(1))=0.d0                                   !right disk
  10  continue
c
      do 30 j=0,n2(1)
        omegn(0,j)=2.d0*phi(1,j)/hr(1,1)**2     !hub
        veln(0,j)=0.d0                          !hub
        omegn(n1(1),j)=2.d0*phi(n1(1)-1,j)/hr(n1(1),1)**2   !shroud
        veln(n1(1),j)=-Re*Ptl*r(n1(1),1)                    !shroud
!        veln(n1(1),j)=0.d0
  30  continue
c
      do 40 i=0,n1(1)
      do 40 j=0,n2(1)
        omegp(i,j)=omeg(i,j)
        omeg(i,j)=omegn(i,j)
        velp(i,j)=vel(i,j)
        vel(i,j)=veln(i,j)
 !       thetap(i,j)=theta(i,j)
 !       theta(i,j)=thetan(i,j)
  40  continue
c
      return
      end
c
c----------------------------------------------------------------------c
c     subroutine cm                                                    c
c     =============                                                    c  
c     calculates error in Cmi at every timestep and output to fort.96  c
c     print overall mean Cmi and e at end of computation.              c
c----------------------------------------------------------------------c 
      subroutine cm(rmax,r0,zmax,zmin,k,nts,ap)
c
      implicit real*8(a-h,o-z)
      parameter (nrmax=250,nzmax=250,na=8)
      parameter (pi=3.14159265359d0,vis=1.82d-5,denref=1.19d0)
      real*8 hr(nrmax,na),hz(nzmax,na),r(0:nrmax,na),z(0:nzmax,na)
      real*8 omegn(0:nrmax,0:nzmax),omeg(0:nrmax,0:nzmax)
      real*8 omegp(0:nrmax,0:nzmax),phi(0:nrmax,0:nzmax)
      real*8 veln(0:nrmax,0:nzmax),vel(0:nrmax,0:nzmax)
      real*8 velp(0:nrmax,0:nzmax) 
      real*8 vel_d(0:nrmax,0:nzmax)
      real*8 hr_d(nrmax),hz_d(nzmax),r_d(0:nrmax),z_d(0:nzmax) 
      real*8 amr(2),amz(2),cmr(2),cmz(2)
      real*8 cmr_mean(2),cmz_mean(2),I_mean
      integer n1(na),n2(na)
      common/val/Re,Ptl !,Ra
      common/sol1/phi,omegn,omeg,omegp,veln,vel,velp
!mpk      common/sol2/thetan,theta,thetap
      common/num/hr,hz,n1,n2,dt,r,z
      real*8 velave_d(1:nrmax,1:nzmax)
c      
      omref_d=Re*Ptl*2.1541d-5/(0.382d0**2)
      rmax_d=0.382d0*rmax
      r0_d=0.382d0*r0
      zmax_d=0.382d0*zmax
      zmin_d=0.382d0*zmin
      r_d(0)=0.382d0*r(0,1)
      z_d(0)=0.382d0*z(0,1)
c
      do 3 i=1,n1(1)
      hr_d(i)=0.382d0*hr(i,1)
      r_d(i)=0.382d0*r(i,1)
    3 continue
c
      do 5 j=1,n2(1)
      hz_d(j)=0.382d0*hz(j,1)
      z_d(j)=0.382d0*z(j,1)
    5 continue
c
      do 6 i=0,n1(1)
      do 7 j=0,n2(1)
      vel_d(i,j)=(vel(i,j)*(2.1541d-5)/0.382d0)+omref_d*r_d(i)
    7 continue
    6 continue   
c     
c-------------------------------------------------------------------
c     show starts:_-| |
c-------------------------------------------------------------------
c
      amz(1)=0.d0
      if(r0 .NE. 0)then                    !foolproof check
c      p=-hr_d(2)/(hr_d(1)*(hr_d(1)+hr_d(2)))
c      q=(hr_d(2)-hr_d(1))/hr_d(1)/hr_d(2)
c      r=hr_d(1)/(hr_d(2)*(hr_d(1)+hr_d(2)))
      do 10 j=1,n2(1)
      vn=(vel_d(2,j)/r_d(2))-(vel_d(0,j)/r_d(0))
      vp=(vel_d(1,j)/r_d(1))-(vel_d(0,j)/r_d(0))
      sswn=vis*r_d(0)*((-vn*hr_d(1)/(hr_d(1)+hr_d(2))
     1     +vp*(hr_d(1)+hr_d(2))/hr_d(1))/hr_d(2))
c      ssw=vis*r_d(0)*(vp/hr_d(1))
c      sswn=vis*r_d(0)*(p*vel_d(0,j)/r_d(0)+q*vel_d(1,j)/r_d(1)
c     1  +r*vel_d(2,j)/r_d(2))
      if(j.eq.1)then
      sswp=sswn
      endif
      ssw=(sswp+sswn)/2d0
      sswp=sswn
      s_temp=ssw*2*pi*r0_d**2*hz_d(j)*2/denref/omref_d**2/rmax_d**5
      print*,z_d(j),s_temp
      amz(1)=amz(1)+ssw*2*pi*r0_d**2*hz_d(j)
10    continue
      cmz(1)=amz(1)*2/denref/omref_d**2/rmax_d**5
      endif      
c
      amz(2)=0.d0
c      p=-hr_d(n1(1))/(hr_d(n1(1)-1)*(hr_d(n1(1)-1)+hr_d(n1(1))))
c      q=(hr_d(n1(1))-hr_d(n1(1)-1))/hr_d(n1(1)-1)/hr_d(n1(1))
c      r=hr_d(n1(1)-1)/(hr_d(n1(1))*(hr_d(n1(1)-1))+hr_d(n1(1))))    
      do 20 j=1,n2(1)
      vs=(vel_d(n1(1)-2,j)/r_d(n1(1)-2))-(vel_d(n1(1),j)/r_d(n1(1)))
      vp=(vel_d(n1(1)-1,j)/r_d(n1(1)-1))-(vel_d(n1(1),j)/r_d(n1(1)))
      sswn=vis*r_d(n1(1))*((-vs*hr_d(n1(1))/(hr_d(n1(1))+hr_d(n1(1)-1))
     1     +vp*(hr_d(n1(1))+hr_d(n1(1)-1))/hr_d(n1(1)))
     2      /hr_d(n1(1)-1))           
c      ssw=vis*r_d(n1(1))*(vp/hr_d(n1(1)))
c      sswn=-vis*r_d(n1(1))*(p*vel_d(n1(1)-2,j)/r_d(n1(1)-2)
c     1  +q*vel_d(n1(1)-1,j)/r_d(n1(1)-1)
c     2  +r*vel_d(n1(1),j)/r_d(n1(1))
      if(j.eq.1)then
      sswp=sswn
      endif
      ssw=(sswp+sswn)/2d0
      sswp=sswn
      s_temp=ssw*2*pi*rmax_d**2*hz_d(j)*2/denref/omref_d**2/rmax_d**5
      print*,z_d(j),s_temp
      amz(2)=amz(2)+ssw*2*pi*rmax_d**2*hz_d(j)
20    continue
      cmz(2)=amz(2)*2/denref/omref_d**2/rmax_d**5
c
      amr(1)=0.d0
c      p=-hz_d(2)/(hz_d(1)*(hz_d(1)+hz_d(2)))
c      q=(hz_d(2)-hz_d(1))/hz_d(1)/hz_d(2)
c      r=hz_d(1)/(hz_d(2)*(hz_d(1)+hz_d(2)))
      do 30 i=1,n1(1)
      ve=vel_d(i,2)-vel_d(i,0)
      vp=vel_d(i,1)-vel_d(i,0)
      sswn=vis*(-ve*hz_d(1)/(hz_d(1)+hz_d(2))
     1     +vp*(hz_d(1)+hz_d(2))/hz_d(1))/hz_d(2)*r_d(i)**2
c      ssw=vis*(vp/hz_d(1))
c      sswn=vis*(p*vel_d(i,0)+q*vel_d(i,1)+r*vel_d(i,2))
      if(i.eq.1)then
      sswp=sswn
      endif
      ssw=(sswp+sswn)/2d0
      sswp=sswn
      print*,r_d(i),ssw*2*pi*hr_d(i)*2/denref/omref_d**2/rmax_d**5
      amr(1)=amr(1)+ssw*2*pi*hr_d(i)
30    continue
      cmr(1)=amr(1)*2/denref/omref_d**2/rmax_d**5
c
      amr(2)=0.d0
c      p=-hz_d(n2(1))/(hz_d(n2(1)-1)*(hz_d(n2(1)-1)+hz_d(n2(1))))
c      q=(hz_d(n2(1))-hz_d(n2(1)-1))/hz_d(n2(1)-1)/hz_d(n2(1))
c      r=hz_d(n2(1)-1)/(hz_d(n2(1))*(hz_d(n2(1)-1))+hz_d(n2(1))))  
      do 40 i=1,n1(1)
      vw=vel_d(i,n2(1)-2)-vel_d(i,n2(1))
      vp=vel_d(i,n2(1)-1)-vel_d(i,n2(1))
      sswn=vis*(-vw*hz_d(n2(1))/(hz_d(n2(1))+hz_d(n2(1)-1))
     1     +vp*(hz_d(n2(1))+hz_d(n2(1)-1))/hz_d(n2(1)))
     2     /hz_d(n2(1)-1)*r_d(i)**2
c      ssw=vis*(vp/hz_d(n2(1)))
c      sswn=-vis*(p*vel_d(i,n2(1)-2)+q*vel_d(i,n2(1)-1)
c     1    +r*vel_d(i,n2(1)))
      if(i.eq.1)then
      sswp=sswn
      endif
      ssw=(sswp+sswn)/2d0
      sswp=sswn
      print*,r_d(i),ssw*2*pi*hr_d(i)*2/denref/omref_d**2/rmax_d**5
      amr(2)=amr(2)+ssw*2*pi*hr_d(i)
40    continue
      cmr(2)=amr(2)*2/denref/omref_d**2/rmax_d**5
c      
c      et=(cmz(1)+cmz(2)+cmr(1)+cmr(2))
c      aet=(ABS(cmz(1))+ABS(cmz(2))+ABS(cmr(1))+
c     1   ABS(cmr(2)))
c      WRITE(98,500)k*dt,e,aet
c  500 FORMAT(3(1x,e20.10E2))
c
c      REWIND 11      
c---------------------------------------------------------------
c      calculate the rate of change of total momentum in fluid(see Owen)
c      ideally its normalised form equals Scmi at every timestep. 
c      it is zero over a long time
c---------------------------------------------------------------
c     METHOD 1
      do i=1,n1(1)
      do j=1,n2(1)
      velave_d(i,j)=0.25d0*(vel_d(i-1,j-1)+vel_d(i-1,j)+
     1 vel_d(i,j)+vel_d(i,j-1))
      enddo
      enddo
c      
      a=0   
      do i=1,n1(1)
      do j=1,n2(1)
      a=a+(velave_d(i,j)*hr_d(i)*hz_d(j)*(r_d(i)-hr_d(i)/2)**2)
      enddo
      enddo
c     METHOD 2
c      a=0
c     do i=1,n1(1)-1
c      do j=1,n2(1)-1
c      a=a+(vel_d(i,j)*(hr_d(i)+hr_d(i+1))*(hz_d(j)+hz_d(j+1))/4
c     1  *r_d(i)**2)
c      enddo
c      enddo
c      
      a=a*2*pi
!by all means, avoid numerical differentiation. when cal the mean.
      b=2*(a-ap)/(dt*6774.24446d0)/omref_d**2/rmax_d**5   
!mpk
c     e=(cmz(1)+cmr(1)+cmr(2)+cmz(2))/cmz(2)
      cmz_mean(1)=(cmz_mean(1)*(k-1)+cmz(1))/k
      cmz_mean(2)=(cmz_mean(2)*(k-1)+cmz(2))/k
      cmr_mean(1)=(cmr_mean(1)*(k-1)+cmr(1))/k
      cmr_mean(2)=(cmr_mean(2)*(k-1)+cmr(2))/k
c     
      if(k.eq.1)then
      a_1st=a
      WRITE(95,*)n1,n2
      WRITE(95,*)omref_d,vis,denref
      WRITE(95,*)rmax_d,r0_d,zmax_d,zmin_d
      WRITE(95,*)r_d(n1(1)),r_d(0),z_d(n2(1)),z_d(0)
      endif    
c      
      if(k.eq.nts)then 
      e=(cmz_mean(1)+cmz_mean(2)+cmr_mean(1)+cmr_mean(2))/
     1  (ABS(cmz_mean(1))+ABS(cmz_mean(2))+ABS(cmr_mean(1))+
     2   ABS(cmr_mean(2)))
      I_mean=(a-a_1st)/(nts*dt)*2/6774.24446d0/omref_d**2/rmax_d**5 
      WRITE(97,*)'          INNER AXIAL WALL'
      WRITE(97,1002)(cmz_mean(1))
      WRITE(97,*)'          OUTER AXIAL WALL'
      WRITE(97,1002)(cmz_mean(2))
      WRITE(97,*)'          LEFT HAND WALL'
      WRITE(97,1002)(cmr_mean(1))
      WRITE(97,*)'          RIGHT HAND WALL'
      WRITE(97,1002)(cmr_mean(2))
      WRITE(97,*)'          MEAN Im'
      WRITE(97,*)I_mean
      WRITE(97,*)e
      endif
 1002 FORMAT('   Mean CM              ',10(1X,E12.4))
c
      c=cmz(1)+cmz(2)+cmr(1)+cmr(2)
      print*,b,-c
      ap=a
      WRITE(96,400)k*dt,-cmz(1),-cmz(2),-cmr(1),-cmr(2),b,a
  400 FORMAT(7(1x,e20.10E2))
c      
      RETURN
      END
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      real*8 function dum1a(v,crr,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),crr(nrmax,4,na)
c
      dum1a=crr(i,1,1)*v(i-1,j)+crr(i,3,1)*v(i+1,j)+crr(i,4,1)*v(i+2,j)
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum2a(v,cr,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),cr(nrmax,3,na)
c
      dum2a=cr(i,1,1)*v(i-1,j)+cr(i,3,1)*v(i+1,j)
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum2(v,cr,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),cr(nrmax,3,na)
c 
      dum2=cr(i,1,1)*v(i-1,j)+cr(i,2,1)*v(i,j)+cr(i,3,1)*v(i+1,j)
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum3a(v,crr,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),crr(nrmax,4,na)
c
      dum3a=crr(i,1,1)*v(i-2,j)+crr(i,2,1)*v(i-1,j)+crr(i,4,1)*v(i+1,j)
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum4a(v,cr,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),cr(nrmax,3,na)
c
      dum4a=cr(i,1,1)*v(i-1,j)+cr(i,3,1)*v(i+1,j)
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum4(v,cr,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),cr(nrmax,3,na) 
c 
      dum4=cr(i,1,1)*v(i-1,j)+cr(i,2,1)*v(i,j)+cr(i,3,1)*v(i+1,j) 
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum5a(v,czz,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),czz(nzmax,4,na)
c
      dum5a=czz(j,1,1)*v(i,j-1)+czz(j,3,1)*v(i,j+1)+czz(j,4,1)*v(i,j+2)
c
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum6a(v,cz,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),cz(nzmax,3,na)
c 
      dum6a=cz(j,1,1)*v(i,j-1)+cz(j,3,1)*v(i,j+1) 
c
      return
      end
c----------------------------------------------------------------------c
      real*8 function dum6(v,cz,i,j)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),cz(nzmax,3,na) 
c  
      dum6=cz(j,1,1)*v(i,j-1)+cz(j,2,1)*v(i,j)+cz(j,3,1)*v(i,j+1)  
c
      return
      end
c----------------------------------------------------------------------c 
      real*8 function dum7a(v,czz,i,j)
c 
      implicit real*8(a-h,o-z) 
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),czz(nzmax,4,na)  
c   
      dum7a=czz(j,1,1)*v(i,j-2)+czz(j,2,1)*v(i,j-1)+czz(j,4,1)*v(i,j+1)   
c
      return 
      end 
c----------------------------------------------------------------------c
      real*8 function dum8a(v,cz,i,j) 
c  
      implicit real*8(a-h,o-z)  
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),cz(nzmax,3,na)   
c    
      dum8a=cz(j,1,1)*v(i,j-1)+cz(j,3,1)*v(i,j+1)    
c
      return
      end
c----------------------------------------------------------------------c 
      real*8 function dum8(v,cz,i,j)  
c   
      implicit real*8(a-h,o-z)   
      parameter(nrmax=250,nzmax=250,na=8)
      real*8 v(0:nrmax,0:nzmax),cz(nzmax,3,na)    
c     
      dum8=cz(j,1,1)*v(i,j-1)+cz(j,2,1)*v(i,j)+cz(j,3,1)*v(i,j+1)     
c
      return 
      end 
c----------------------------------------------------------------------c
c======================================================================c
      subroutine poisson(psi,rhs,nlev,vel) 
c  
c     Multigrid subroutine for the Poisson equation, 4 levels, G-S relax
c    ------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
      parameter(kr=250,kz=250,mgen=100000,na=8)
      real*8 psi(0:kr,0:kz),rhs(0:kr,0:kz)
      real*8 v(mgen),f(mgen),vel(0:kr,0:kz)
      real*8 psic(0:kr,0:kz),psif(0:kr,0:kz)
      real*8 rc(0:kr,0:kz),rf(0:kr,0:kz),tempa(0:kr,0:kz)
      real*8 hr(kr,na),hz(kz,na),hhr(kr,na),hhz(kz,na)
      real*8 cr(kr,3,na),crr(kr,4,na),cz(kz,3,na),czz(kz,4,na)
      real*8 r(0:kr,na),z(0:kz,na)
      integer ibeg(na),n1(na),n2(na)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /num/hr,hz,n1,n2,dt,r,z
      common /step/hhr,hhz
      common /ints/ibeg
      common /coeff/cr,crr,cz,czz
c
c
      do 3 l=1,nlev
        do 1 i=1,n1(l)
          hhr(i,l)=hr(i,l)**2
    1   continue
        do 2 j=1,n2(l)
          hhz(j,l)=hz(j,l)**2
    2   continue
    3 continue
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
c       print*,'Sweep ',i
        if(i.eq.1.or.nlev.eq.1)then
          call relax(nrel1,1,nconv,nlev)
        end if
        if(nlev.gt.1)then
          do 51 ilev=1,nlev-1
            call corser(ilev)
            call relax(nrel1,ilev+1,nconv,nlev)
   51     continue
          do 52 ilev=nlev,2,-1
            call finer(ilev)
            call relax(nrel2,ilev-1,nconv,nlev)
   52     continue
        end if
        if(nconv.eq.1)then
          print*,'Sweeps converged',i,'sweeps'
          do 60 ii=0,n1(1)
          do 60 jj=0,n2(1)
            psi(ii,jj)=psif(ii,jj)
   60     continue
c         call rufcon2(psi,vel,n1(1),n2(1),amax,amin)
          return
        end if
   50 continue
      return
      end
c==================================================
      subroutine relax(nrel,klev,nconv,nlev)
      implicit real*8(a-h,o-z)
      parameter(kr=250,kz=250,mgen=100000,na=8)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz),psif(0:kr,0:kz)
      real*8 rc(0:kr,0:kz),rf(0:kr,0:kz),tempa(0:kr,0:kz)
      real*8 hr(kr,na),hz(kz,na),hhr(kr,na),hhz(kz,na)
      real*8 cr(kr,3,na),crr(kr,4,na),cz(kz,3,na),czz(kz,4,na)
      real*8 r(0:kr,na),z(0:kz,na)
      real*8 pdmaa(kr),pdmab(kr),pdmac(kr),pdmar(kr),pdmas(kr)
      real*8 pdmad(kr),pdmae(kr)
      integer ibeg(na),n1(na),n2(na)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /num/hr,hz,n1,n2,dt,r,z
      common /step/hhr,hhz
      common /ints/ibeg
      common /coeff/cr,crr,cz,czz
c
      tolcon=1d-6
      nconv=0
c
      call trans2(f,rf,ibeg(klev),n1(klev),n2(klev))
      call trans2(v,psif,ibeg(klev),n1(klev),n2(klev))
c
c------------------------------relaxing with some G & S.
c
      do 200 irel=1,nrel
c                                 line relax in r-dir
        do 10 j=1,n2(klev)-1
        do 11 i=1,n1(klev)-1
          if(i.le.n1(klev)/2)then
            pdmaa(i)=0.d0
            pdmab(i)=crr(i,1,klev)+cr(i,1,klev)/r(i,klev)
            pdmac(i)=crr(i,2,klev)+cr(i,2,klev)/r(i,klev)
     +              -1.d0/r(i,klev)**2
            pdmad(i)=crr(i,3,klev)+cr(i,3,klev)/r(i,klev)
            pdmae(i)=crr(i,4,klev)
          else
            pdmaa(i)=crr(i,1,klev)
            pdmab(i)=crr(i,2,klev)+cr(i,1,klev)/r(i,klev)
            pdmac(i)=crr(i,3,klev)+cr(i,2,klev)/r(i,klev)
     +              -1.d0/r(i,klev)**2
            pdmad(i)=crr(i,4,klev)+cr(i,3,klev)/r(i,klev)
            pdmae(i)=0.d0
          endif
c
          if(j.le.n2(klev)/2)then
            pdmac(i)=pdmac(i)+czz(j,2,klev)
            pdmar(i)=-(czz(j,1,klev)*psif(i,j-1)+czz(j,3,klev)
     +           *psif(i,j+1)+czz(j,4,klev)*psif(i,j+2))+rf(i,j)
          else
            pdmac(i)=pdmac(i)+czz(j,3,klev)
            pdmar(i)=-(czz(j,1,klev)*psif(i,j-2)+czz(j,2,klev)
     +           *psif(i,j-1)+czz(j,4,klev)*psif(i,j+1))+rf(i,j)
          endif
c
   11   continue
        call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n1(klev)-1)
        do 12 i=1,n1(klev)-1
          psif(i,j)=pdmas(i)
   12   continue
   10   continue
c                                 line relax in z-dir
        do 13 i=1,n1(klev)-1
        do 14 j=1,n2(klev)-1
          if(j.le.n2(klev)/2)then
            pdmaa(j)=0.d0
            pdmab(j)=czz(j,1,klev)
            pdmac(j)=czz(j,2,klev)-1.d0/r(i,klev)**2
            pdmad(j)=czz(j,3,klev)
            pdmae(j)=czz(j,4,klev)
          else
            pdmaa(j)=czz(j,1,klev)
            pdmab(j)=czz(j,2,klev)
            pdmac(j)=czz(j,3,klev)-1.d0/r(i,klev)**2
            pdmad(j)=czz(j,4,klev)
            pdmae(j)=0.d0
          endif
c
          if(i.le.n1(klev)/2)then
            pdmac(j)=pdmac(j)+crr(i,2,klev)+cr(i,2,klev)/r(i,klev)
            pdmar(j)=-(crr(i,1,klev)*psif(i-1,j)+crr(i,3,klev)
     +              *psif(i+1,j)+crr(i,4,klev)*psif(i+2,j)+(cr(i,1,klev)
     +              *psif(i-1,j)+cr(i,3,klev)*psif(i+1,j))/r(i,klev))
     +              +rf(i,j)
          else
            pdmac(j)=pdmac(j)+crr(i,3,klev)+cr(i,2,klev)/r(i,klev)
            pdmar(j)=-(crr(i,1,klev)*psif(i-2,j)+crr(i,2,klev)
     +              *psif(i-1,j)+crr(i,4,klev)*psif(i+1,j)+(cr(i,1,klev)
     +              *psif(i-1,j)+cr(i,3,klev)*psif(i+1,j))/r(i,klev))
     +              +rf(i,j)
          endif
c
   14   continue
        call pdma(pdmaa,pdmab,pdmac,pdmad,pdmae,pdmar,pdmas,n2(klev)-1)
        do 15 j=1,n2(klev)-1
          psif(i,j)=pdmas(j)
   15   continue
   13   continue
c
c                                  Calculating residual.
c
        resmax=0.d0
        do 20 i=1,n1(klev)-1
        do 20 j=1,n2(klev)-1
          if(i.le.n1(klev)/2)then
            temp=crr(i,1,klev)*psif(i-1,j)+crr(i,2,klev)*psif(i,j)
     +          +crr(i,3,klev)*psif(i+1,j)+crr(i,4,klev)*psif(i+2,j)
     +          +(cr(i,1,klev)*psif(i-1,j)+cr(i,2,klev)*psif(i,j)
     +          +cr(i,3,klev)*psif(i+1,j))/r(i,klev)-psif(i,j)
     +          /r(i,klev)**2
          else
            temp=crr(i,1,klev)*psif(i-2,j)+crr(i,2,klev)*psif(i-1,j)
     +          +crr(i,3,klev)*psif(i,j)+crr(i,4,klev)*psif(i+1,j)
     +          +(cr(i,1,klev)*psif(i-1,j)+cr(i,2,klev)*psif(i,j)
     +          +cr(i,3,klev)*psif(i+1,j))/r(i,klev)-psif(i,j)
     +          /r(i,klev)**2
          endif
c
          if(j.le.n2(klev)/2)then
            temp=temp+czz(j,1,klev)*psif(i,j-1)+czz(j,2,klev)*psif(i,j)
     +          +czz(j,3,klev)*psif(i,j+1)+czz(j,4,klev)*psif(i,j+2)
          else
           temp=temp+czz(j,1,klev)*psif(i,j-2)+czz(j,2,klev)*psif(i,j-1)
     +          +czz(j,3,klev)*psif(i,j)+czz(j,4,klev)*psif(i,j+1)
          endif
c
          temp=rf(i,j)-temp
          rc(i,j)=temp
          temp=dabs(temp)
          if(temp.gt.resmax)resmax=temp
   20   continue
c
c       if(klev.eq.1)write(6,100)klev,resmax
  100   format(1x,'Level',i1,'  resmax=',f16.10)
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
      parameter(kr=250,kz=250,mgen=100000,na=8)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz),psif(0:kr,0:kz)
      real*8 rc(0:kr,0:kz),rf(0:kr,0:kz),tempa(0:kr,0:kz)
      real*8 hr(kr,na),hz(kz,na),hhr(kr,na),hhz(kz,na)
      real*8 r(0:kr,na),z(0:kz,na)
      integer ibeg(na),n1(na),n2(na)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /num/hr,hz,n1,n2,dt,r,z
      common /step/hhr,hhz
      common /ints/ibeg
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
      parameter(kr=250,kz=250,mgen=100000,na=8)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kz),psif(0:kr,0:kz)
      real*8 rc(0:kr,0:kz),rf(0:kr,0:kz),tempa(0:kr,0:kz)
      real*8 hr(kr,na),hz(kz,na),hhr(kr,na),hhz(kz,na)
      real*8 r(0:kr,na),z(0:kz,na)
      integer ibeg(na),n1(na),n2(na)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /num/hr,hz,n1,n2,dt,r,z
      common /step/hhr,hhz
      common /ints/ibeg
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
      parameter(kr=250,kz=250,mgen=100000)
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
      parameter(kr=250,kz=250,mgen=100000)
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
      parameter(kr=250,kz=250,mgen=100000)
      real*8 a2d1(0:kr,0:kz),a2d2(0:kr,0:kz)
c
      do 10 j=0,n2
      do 10 i=0,n1
        a2d2(i,j)=a2d1(i,j)
   10 continue
      return
      end
c==================================================
      subroutine rufcon2(a,b,nx,ny) !,amax,amin)
      parameter(mx=250,my=250,mgen=100000)
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
      print*,'phimax,min=',amax,amin
      print*,'omegmax,min=',bmax,bmin
c      adiff=amax-amin
!mpk      bdiff=bmax-bmin
c      if(adiff.ne.0.)adiff=1./adiff
!mpk      if(bdiff.ne.0.)bdiff=1./bdiff
c      do 30 j=0,ny,2
c      jj=ny-j
c      do 10 i=0,nx,2
c      xa=1.5+10.*(a(i,jj)-amin)*adiff
!mpk      xb=1.5+10.*(b(i,jj)-bmin)*bdiff
c      na=xa
!mpk   nb=xb
c      c(i)=ch(na)
!mpk      c(i+nx+2)=ch(nb)
c  10  continue
c      m2=nx+1
c      m1=0
c      if(nx.gt.130)m1=nx*2+2-130
c      write(6,20)(c(k),k=m1,m2,2)
c  20  format(1x,200a1)
c  30  continue
      return
      end
c==================================================
c====to search for minimum entry in an array=======
      double precision function titch(a,nr,nz)
      parameter(kr=250,kz=250,mgen=100000)
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
c=====to search for maximum entry in an array======
      double precision function woppa(a,nr,nz)
      parameter(kr=250,kz=250,mgen=100000)
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
c==================================================
