c    from disk13/enptwl/proj9a.dir
      program research

      implicit real*8(a-h,o-z)
      parameter(kr=256,kp=256,g=9.81,mgen=100000)
      real*8 thp(0:kr,0:kp),thm(0:kr,0:kp),th(0:kr,0:kp)
      real*8 psi(0:kr,0:kp),psip(0:kr,0:kp),psim(0:kr,0:kp)
      real*8 vorp(0:kr,0:kp),vor(0:kr,0:kp),vorm(0:kr,0:kp)
      real*8 thold(0:kr,0:kp),vorold(0:kr,0:kp),Qconv(100000)
      real*8 ssumc(10,mgen),ssums(10,mgen)
      character*20 fname,fname1,fname2,fname3,fname4,fname5,ans
      common /points/nr,nphi,dphi,dr,rmin,rmax,dt,nt
      common /vortic/vor,vorm,vorp,w,Pr,Ra,Re
      common /temper/th,thm,thp
      common /converg/thold,vorold
      common /stream/psi
      common /filenames/fname1,fname2,fname3,fname4
      common /iiii/it
      common /four/ssumc,ssums
     
      m=1                   
      print*,'Enter Num steps in radial direction(128 max):'
      read*,nr
      print*,'Enter No step in phi direction(192 max):'
      read*,nphi
      print*,'Enter radius of inner ring:'
      read*,rmin
      print*,'Enter radius of outer ring:'
      read*,rmax
      print*,'Enter the time step and number of time steps:'
      read*,dt,nt
      print*,'iterations to start: '
      read*,istart
      print*,'write interval: '
      read*,nf
      print*,'Enter Raleigh number,Reynolds no,:'
      read*,Ra,Re
c      print*,'Enter the Prandtl number(approx 0.7 prob):'
c      read*,Pr
      Pr=0.7d0

****** Specify Filenames to Write Results to ***********************

      print*,'Enter Filename Vorticity'
      read*,fname1
      print*,'Enter Filename Temperature'
      read*,fname2
      print*,'Enter Filename Stream Function'
      read*,fname3
      print*,'Enter Filename Nusselt Thru time'
      read*,fname4
c
      pi=dacos(-1.d0)
      dphi=2.d0*pi/nphi
      dr=(rmax-rmin)/nr
c      
***   Set initial values of vorticity,stream function ***
***   and temperatures at t(n),t(n-1) to 0,1 on 0     ***
      pi=dacos(-1.d0)
      do 10 i=0,nr
      do 10 j=0,nphi
        r=rmin+dble(i)*dr
        vorm(i,j)=0.d0
        vor(i,j)=0.d0
        vorold(i,j)=0.d0
        psi(i,j)=0.d0 
        psip(i,j)=0.d0
        psim(i,j)=0.d0
        thm(i,j)=1.d0-dlog(r)/dlog(rmin)
        th(i,j)=1.d0-dlog(r)/dlog(rmin)
        thp(i,j)=1.d0-dlog(r)/dlog(rmin)
        thold(i,j)=1.d0-dlog(r)/dlog(rmin)
   10 continue
c
   11 continue  
      print*,'n,amp'
      read*,nval,amp
      if(nval.eq.0)goto 13
      do 12 i=0,nr
      do 12 j=0,nphi
        r=rmin+dble(i)*dr
        phi=8.d0*datan(1.d0)*dble(j)/dble(nphi)
        thm(i,j)=thm(i,j)+amp*(rmax-r)*(r-rmin)*dcos(dble(nval)*phi)
        th(i,j)=thm(i,j)
        thp(i,j)=thm(i,j)
        thold(i,j)=thm(i,j)
   12 continue
      goto 11
   13 continue
c     call rufcon2(psi,thp,nr,nphi)
c
!mpk      th(nr/4,nphi/4)=th(nr/4,nphi/4)-0.4d0*th(nr/4,nphi/4)
!mpk      thm(nr/4,nphi/4)=thm(nr/4,nphi/4)-0.4d0*thm(nr/4,nphi/4)
      th(nr/4,nphi/2)=th(nr/4,nphi/2)-0.02d0
      thm(nr/4,nphi/2)=thm(nr/4,nphi/2)-0.02d0
!mpk      th(nr/4,nphi/4*3)=th(nr/4,nphi/4*3)-0.4d0*th(nr/4,nphi/4*3)
!mpk      thm(nr/4,nphi/4*3)=thm(nr/4,nphi/4*3)-0.4d0*thm(nr/4,nphi/4*3)
!mpk      th(nr/4,nphi/2)=th(nr/4,nphi/2)-0.4d0*th(nr/4,nphi/2)
!mpk      thm(nr/4,nphi/2)=thm(nr/4,nphi/2)-0.4d0*thm(nr/4,nphi/2)
!mpk      th(nr/4,0)=th(nr/4,0)-0.4d0*th(nr/4,0)
!mpk      thm(nr/4,0)=thm(nr/4,0)-0.4d0*thm(nr/4,0)
!mpk      call rufcon(th,nr,nphi)
c
      print*,'read initial conditions from file (y/n)?'
      read*,ans
      if(ans.eq.'yes'.or.ans.eq.'y')then
        print*,'Enter Filename Vorticity'
        read*,fname5
        open(11,file=fname5)
        read(11,*)nr,nphi
        read(11,*)rmax,rmin,phimax,phimin
        read(11,*)((vor(i,j),i=0,nr),j=0,nphi)
        close(11)
        print*,'Enter Filename Temperature'
        read*,fname5
        open(11,file=fname5)
        read(11,*)nr,nphi
        read(11,*)rmax,rmin,phimax,phimin
        read(11,*)((th(i,j),i=0,nr),j=0,nphi)
        close(11)
        print*,'Enter Filename Stream Function'
        read*,fname5 
        open(11,file=fname5)
        read(11,*)nr,nphi
        read(11,*)rmax,rmin,phimax,phimin
        read(11,*)((psi(i,j),i=0,nr),j=0,nphi)
        close(11)
        do 95 i=0,nr
        do 95 j=0,nphi
          vorm(i,j)=vor(i,j)
          thm(i,j)=th(i,j)
   95   continue
      endif
c
  100   continue
c
********************************************************************
***** Start Of Main Time Loop Follows ******************************
********************************************************************

      do 5 it=1,nt
       
         print*,'time =',dble(it)*dt

***   Set values of vibration in x,y direction respectively ***       
c        
c
***   Call subroutines to calculate temp,vort,stream  ***
c
         call tempr
         call vorticity
         call aaaaaaaa(it)
         call poisson(psi,vorp,nr,nphi,rmax,rmin)
c      
***   Update values of vor,vorm,th,thm                ***
c
	 do 30 j=0,nphi
           vorp(0,j)=2.d0*psi(1,j)/dr**2
           vorp(nr,j)=2.d0*psi(nr-1,j)/dr**2
           do 30 i=0,nr
             vorm(i,j)=vor(i,j)
             vor(i,j)=vorp(i,j)
             thm(i,j)=th(i,j)
             th(i,j)=thp(i,j)
   30    continue

********* Call Subroutines To Show To Screen A Representation ***
*********  Of Stream Function, Temperature and Vorticity ********

c         call rufcon(psi,nr,nphi)
          call rufcon(th,nr,nphi)
c         call rufcon(vor,nr,nphi)
c         call rufcon2(psi,thp,nr,nphi)

***   Call Nusselt Subroutine To Calculate Nusselt Number     ***

      call Nusselt(th,Qconv,it)


***   At The End of the timestep                      ***
***   Call convergence subroutine to check for        ***
***   convergence, if converged write data to file    ***
c
        if(k.gt.istart)then
          if(k/nf.eq.dble(k)/dble(nf))then
            print*,'writing to file'
            numb=(k-istart)/nf
c
            write(fname,400)'5e6temp',numb
  400       format(A7,I3.3)
            open(unit=400,file=fname)
            write(400,*)nr,nphi
            write(400,*)rmax,r0,2.d0*pi,0
            write(400,*)((th(i,j),i=0,nr),j=0,nphi)
            close(400)
          endif
        endif
c
        call conver(Qconv,it)

   5  continue    
      stop
      end
*******************************************************
      subroutine tempr
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592654,kr=256,kp=256)
      real*8 thp(0:kr,0:kp),thm(0:kr,0:kp),th(0:kr,0:kp)
      real*8 psi(0:kr,0:kp)
      real*8 LHS
      common /points/nr,nphi,dphi,dr,rmin,rmax,dt,nt
      common /temper/th,thm,thp
      common /stream/psi

***  Start of Loops to Calculate Temperature At each Grid Point **

      do 60 i=1,nr-1
      do 60 j=0,nphi
        im=i-1
        ip=i+1
        jm=j-1
        jp=j+1
        if(j .eq. 0)jm=nphi-1
        if(j .eq. nphi)jp=1
        ri=rmin+dble(i)*dr

        arak=arakawa(psi,th,i,j)/ri

        LHS=1.d0/(2.d0*dt)+1.d0/dr**2+1.d0/(ri*dphi)**2

        RHS=(1.d0/(2.d0*dt)-1.d0/dr**2-1.d0/(ri*dphi)**2)*thm(i,j)
     +     +(th(ip,j)+th(im,j))/dr**2
     +     +(th(ip,j)-th(im,j))/(2.d0*dr*ri)
     +     +(th(i,jp)+th(i,jm))/(ri*dphi)**2
     +     +arak

        thp(i,j)=RHS/LHS

   60 continue
   
      return
      end

********************************************************
      subroutine vorticity
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592654,kr=256,kp=256)
      real*8 vorp(0:kr,0:kp),vor(0:kr,0:kp),vorm(0:kr,0:kp)
      real*8 psi(0:kr,0:kp),th(0:kr,0:kp)
      real*8 thm(0:kr,0:kp),thp(0:kr,0:kp)
      real*8 LHS
      common /points/nr,nphi,dphi,dr,rmin,rmax,dt,nt
      common /vortic/vor,vorm,vorp,w,Pr,Ra,Re
      common /temper/th,thm,thp
      common /stream/psi
      
***  Start Of Loops To Calculate Vorticity At Each Grid Point **

      do 70 i=1,nr-1
      do 70 j=0,nphi
        ri=rmin+dble(i)*dr
        phi=dphi*dble(j)
        jm=j-1
        jp=j+1
        im=i-1
        ip=i+1
        if(j .eq. 0)jm=nphi-1
        if(j .eq. nphi)jp=1
        arak1=arakawa(psi,vor,i,j)/(ri*Pr)
        arak2=arakawa(psi,th,i,j)
        arak2=2.d0*arak2*Ra/(Re*Pr*ri)

        LHS=1.d0/(2.d0*dt*Pr)+1.d0/dr**2+1.d0/(ri*dphi)**2       

        buoy=-Ra*(th(i,jp)-th(i,jm))/(2.d0*dphi)

        Rest=(1.d0/(2.d0*dt*Pr)-1.d0/dr**2-1.d0/(ri*dphi)**2)
     +      *vorm(i,j)
     +      +(vor(ip,j)+vor(im,j))/dr**2
     +      +(vor(ip,j)-vor(im,j))/(2.d0*dr*ri)
     +      +(vor(i,jp)+vor(i,jm))/(ri*dphi)**2
     +      +arak1+arak2
  
        RHS=buoy+Rest
        vorp(i,j)=RHS/LHS
  70  continue
c
c     call rufcon(vorp,nr,nphi)
c     
      return
      end

**********************************************************
      subroutine conver(Qconv,it)
      implicit real*8(a-h,o-z)
      parameter(kr=256,kp=256,pi=3.141592654,mgen=100000)
      real*8 vor(0:kr,0:kp),vorold(0:kr,0:kp)
      real*8 th(0:kr,0:kp),thold(0:kr,0:kp)
      real*8 psi(0:kr,0:kp),vorm(0:kr,0:kp),vorp(0:kr,0:kp)
      real*8 thm(0:kr,0:kp),thp(0:kr,0:kp),Qconv(100000)
      real*8 ssumc(10,mgen),ssums(10,mgen)
      character*20 fname1,fname2,fname3,fname4
      common /converg/thold,vorold
      common /temper/th,thm,thp
      common /vortic/vor,vorm,vorp,w,Pr,Ra,Re
      common /points/nr,nphi,dphi,dr,rmin,rmax,dt,nt
      common /stream/psi
      common /filenames/fname1,fname2,fname3,fname4
      common /four/ssumc,ssums
      

      convor=0.d0
      conth=0.d0
      do 90 i=0,nr
      do 90 j=0,nphi
        vergvor=dabs(vor(i,j)-vorold(i,j))
        if(vergvor .gt. convor)then
          convor=vergvor   
        endif
        vorold(i,j)=vor(i,j)
      
        vergth=dabs(th(i,j)-thold(i,j))
        if(i .gt. 0 .and. i .lt. nr)then
          if(vergth .gt. conth) then
            conth=vergth
          endif
          thold(i,j)=th(i,j) 
        endif
   90 continue

*** Print to Screen Convergence Values             ****

      print*,'*******************************************************'
      print*,convor,' Convor',conth,' Conth',it
      print*,'*******************************************************'
      
      iend=0
   
      if(conth .le. 1.d-6.and.it.gt.50 )then
        print*,'You have convergence'
        iend=1
      end if
      if(it.eq.nt)then
        print*,'Thats it, mate! Out of time.'
        iend=1
      end if
c
c     if(it.ge.5)then
c     if(dabs(Qconv(it)-Qconv(it-1)).lt.1d-8
c    +    .and.dabs(Qconv(it)-Qconv(it-2)).lt.1d-8)then
c         print*,'converged'
c         iend=1
c     endif
c     endif
c
      if(iend.eq.1)then
c
      open(12,file=fname1)
      open(13,file=fname2)
      open(14,file=fname3)
      open(15,file=fname4)
c
        write(12,*)nr,nphi
        write(12,*)rmax,rmin,2*pi,0.d0
        write(13,*)nr,nphi
        write(13,*)rmax,rmin,2*pi,0.d0
        write(14,*)nr,nphi
        write(14,*)rmax,rmin,2*pi,0.d0

***    Once Convergence Found Then Results Written to File **

          write(12,*)((vor(i,j),i=0,nr),j=0,nphi)
          write(13,*)((th(i,j),i=0,nr),j=0,nphi)
          write(14,*)((psi(i,j),i=0,nr),j=0,nphi)
c
!mpk        write(15,200)(l*dt,Qconv(l),ssumc(1,l),ssums(1,l),
!mpk     +  ssumc(2,l),ssums(2,l),ssumc(3,l),ssums(3,l),ssumc(4,l),
!mpk     +  ssums(4,l),ssumc(5,l),ssums(5,l),ssumc(6,l),ssums(6,l),
!mpk     +  ssumc(7,l),ssums(7,l),ssumc(8,l),ssums(8,l),ssumc(9,l),
!mpk     +  ssums(9,l),ssumc(10,l),ssums(10,l),l=1,it,10)
!mpkc200  format(22(1x,e16.8))
        write(15,200)(l*dt,Qconv(l),l=1,it,10)
 200  format(2(1x,e16.8))
c
        close(12)
        close(13)
        close(14)
        close(15)
c
        stop
      endif
c
      return
      end

***********************************************************************
      subroutine Nusselt(th,Qconv,it)
      implicit real*8(a-h,o-z)
      parameter(kr=256,kp=256)
      real*8 Qconv(100000),th(0:kr,0:kp)
      common /points/nr,nphi,dphi,dr,rmin,rmax,dt,nt

      Anul=0.d0
      ht=0.d0
      pi=dacos(-1.d0)

***  Calculate Local, and Consequently Global Nusselt Number **
c
      do 10 j=1,nphi
        ht=ht+dphi*(-25.d0*th(0,j)+48.d0*th(1,j)-36.d0
     +       *th(2,j)+16.d0*th(3,j)-3.d0*th(4,j))/(12.d0*dr)
  10  continue

      Qconv(it)=-rmin*dlog(rmin)*ht/(2.d0*pi)
      print*,'Nusselt no. = ',Qconv(it)

      return
      end
***********************************************************************
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     subroutine INT                                                   c
c     ==============                                                   c
c     Calculates the volume integrals                                  c
c----------------------------------------------------------------------c
c
      subroutine aaaaaaaa(kt)
c
      implicit real*8(a-h,o-z)
      parameter(kr=256,kp=256,mgen=100000)
      real*8 Qconv(100000)
      real*8 ssumc(10,mgen),ssums(10,mgen)
      real*8 thp(0:kr,0:kp),thm(0:kr,0:kp),th(0:kr,0:kp)
      common /temper/th,thm,thp
      common /points/nr,nphi,dphi,dr,rmin,rmax,dt,nt
      common /four/ssumc,ssums

c
      do 10 n=1,10
        ssumc(n,kt)=0.d0
        ssums(n,kt)=0.d0
        sumc=0.d0
        sums=0.d0
        do 20 k=1,nphi
          ang=dble(k)*dphi
          do 30 i=0,nr+1
            r=rmin+dble(i)*dr
            fact=1.d0
            if(i.eq.0.or.i.eq.nr+1)fact=fact/2.d0
            fc=r*th(i,k)*dcos(dble(n)*ang)
            fs=r*th(i,k)*dsin(dble(n)*ang)
            sumc=sumc+fc*fact
            sums=sums+fs*fact
   30     continue
          sumc=sumc*dr
          sums=sums*dr
c
          ssumc(n,kt)=ssumc(n,kt)+sumc*dphi
          ssums(n,kt)=ssums(n,kt)+sums*dphi
   20   continue
   10 continue
c
      return
      end
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
      double precision function arakawa(psi,fu,i,j)
      implicit real*8(a-h,o-z)
      parameter(kr=256,kp=256)
      real*8 psi(0:kr,0:kp),fu(0:kr,0:kp)
      common /points/nr,nphi,dphi,dr,rmin,rmax,dt,nt
   
      ip=i+1
      im=i-1
      jp=j+1
      jm=j-1
      if(j .eq. 0)jm=nphi-1
      if(j .eq. nphi)jp=1
      arak=(psi(ip,j)-psi(im,j))*(fu(i,jp)-fu(i,jm))
     +    -(psi(i,jp)-psi(i,jm))*(fu(ip,j)-fu(im,j))
     +    +psi(ip,j)*(fu(ip,jp)-fu(ip,jm))
     +    -psi(im,j)*(fu(im,jp)-fu(im,jm))
     +    -psi(i,jp)*(fu(ip,jp)-fu(im,jp))
     +    +psi(i,jm)*(fu(ip,jm)-fu(im,jm))
     +    +fu(i,jp)*(psi(ip,jp)-psi(im,jp))
     +    -fu(i,jm)*(psi(ip,jm)-psi(im,jm))
     +    -fu(ip,j)*(psi(ip,jp)-psi(ip,jm))
     +    +fu(im,j)*(psi(im,jp)-psi(im,jm))
          arak=arak/(12.d0*dr*dphi)
          arakawa=arak

       return
       end


c======================================================================
      subroutine poisson(psi,rhs,nr,nth,rmax,rmin)
c
c     Multigrid subroutine for the Poisson equation, up to 6 levels, G-S relax
c     ------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      parameter(kr=256,kth=256,mgen=50000)
      real*8 psi(0:kr,0:kth),rhs(0:kr,0:kth)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kth),psif(0:kr,0:kth)
      real*8 rc(0:kr,0:kth),rf(0:kr,0:kth),tempa(0:kr,0:kth)
      real*8 hr(6),hth(6),hhr(6),hhth(6)
      real*8 rval(0:kr)
      integer ibeg(6),n1(6),n2(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hth,hhr,hhth,rval
      common /ints/ibeg,n1,n2
c
c---------------Determining number of multigrid levels
c
      nlev=0
      itest=nr
      itest2=nr
      jtest=nth
      jtest2=nth
c      print*,itest,itest2,jtest,jtest2
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
        if(jtest.eq.1)then
          nlev=nlev-1
          goto 2
        end if
    1 continue
    2 continue
      if(nlev.gt.6)nlev=6
c      print*,'Number of levels=',nlev
c
      dr=(rmax-rmin)/dble(nr)
      thmax=8.d0*datan(1.d0)
      dth=thmax/dble(nth)
      nvs=10000
      nrel1=2
      nrel2=2
c
      do 10 i=1,mgen
        v(i)=0.d0
        f(i)=0.d0
   10 continue
c
      ii=1
      do 20 i=1,nlev
        n1(i)=nr*2/(2**i)
        n2(i)=nth*2/(2**i)
        hr(i)=dr*dble(2**i)/2.d0
        hth(i)=dth*dble(2**i)/2.d0
        hhr(i)=hr(i)**2
        hhth(i)=hth(i)**2
        ibeg(i)=ii
        ii=ii+(n1(i)+1)*(n2(i)+1)+1
c        print*,i,n1(i),n2(i),hr(i),hth(i),ibeg(i)
   20 continue
c
      do 30 i=0,nr
        rval(i)=dble(i)*dr+rmin
        do 30 j=0,nth
          psic(i,j)=0.d0
          psif(i,j)=psi(i,j)
          rc(i,j)=0.d0
          rf(i,j)=0.d0
   30 continue
c
      do 40 i=0,nr
      do 40 j=0,nth
        rf(i,j)=rhs(i,j)
   40 continue
c
      call trans1(rf,f,ibeg(1),n1(1),n2(1))
      call trans1(psif,v,ibeg(1),n1(1),n2(1))
c
      do 50 i=1,nvs
c        print*,'Sweep ',i
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
        if(nconv.eq.1.or.i.eq.nvs)then
c          print*,'Converged'
          do 60 ii=0,nr
          do 60 jj=0,nth
            psi(ii,jj)=psif(ii,jj)
   60     continue
          print*,'Sweep ',i          
c         call rufcon(psi,nr,nth)
          return
        end if
   50 continue
      return
      end
c==================================================
      subroutine relax(nrel,klev,nconv,nlev)
      implicit real*8(a-h,o-z)
      parameter(kr=256,kth=256,mgen=50000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kth),psif(0:kr,0:kth)
      real*8 rc(0:kr,0:kth),rf(0:kr,0:kth),tempa(0:kr,0:kth)
      real*8 hr(6),hth(6),hhr(6),hhth(6)
      real*8 rval(0:kr)
c
      real*8 tdmaa(kr),tdmab(kr),tdmac(kr),tdmar(kr),tdmas(kr)
c
      integer ibeg(6),n1(6),n2(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hth,hhr,hhth,rval
      common /ints/ibeg,n1,n2
      common /iiii/it
c
      tolcon=1.d-6
c     if(it.lt.40)tolcon=1.d-4
c     if(it.lt.30)tolcon=5.d-3
c     if(it.lt.100)tolcon=1.d-3
 
      nconv=0
c
      call trans2(f,rf,ibeg(klev),n1(klev),n2(klev))
      call trans2(v,psif,ibeg(klev),n1(klev),n2(klev))
c
c------------------------------relaxing with some G & S.
c
      levfact=(2**klev)/2
c
      do 200 irel=1,nrel
c                                 line relax in r-dir
        do 10 j=1,n2(klev)
          jp=j+1
          jm=j-1
          if(j.eq.0)jm=n2(klev)-1
          if(j.eq.n2(klev))jp=1
          do 11 i=1,n1(klev)-1
            ii=i*levfact
            tdmaa(i)=1.d0/hhr(klev)-0.5d0/hr(klev)/rval(ii)
            tdmac(i)=1.d0/hhr(klev)+0.5d0/hr(klev)/rval(ii)
            tdmab(i)=-2.d0/hhr(klev)-2.d0/hhth(klev)/rval(ii)**2
            tdmar(i)=-(psif(i,jp)+psif(i,jm))/hhth(klev)/rval(ii)**2
     +              +rf(i,j)
c123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
   11     continue
          call tdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n1(klev)-1)
          do 12 i=1,n1(klev)-1
            psif(i,j)=tdmas(i)
            if(j.eq.n2(klev))psif(i,0)=psif(i,j)
   12     continue
   10   continue
c                                 line relax in th-dir
        do 13 i=1,n1(klev)-1
          ii=i*levfact
          ip=i+1
          im=i-1
          do 14 j=1,n2(klev)
            tdmaa(j)=1.d0/hhth(klev)/rval(ii)**2
            tdmac(j)=1.d0/hhth(klev)/rval(ii)**2
            tdmab(j)=-2.d0/hhr(klev)-2.d0/hhth(klev)/rval(ii)**2
            tdmar(j)=-(psif(ip,j)+psif(im,j))/hhr(klev)+rf(i,j)
     +              -(psif(ip,j)-psif(im,j))*0.5d0/rval(ii)/hr(klev)
   14     continue
          call ptdma(tdmaa,tdmab,tdmac,tdmas,tdmar,n2(klev))
          do 15 j=0,n2(klev)
            if(j.eq.0)then
              psif(i,j)=tdmas(n2(klev))
            else
              psif(i,j)=tdmas(j)
            end if
   15     continue
   13   continue
c
c                                  Calculating residual.
c
        resmax=0.d0
        do 20 i=1,n1(klev)-1
        do 20 j=0,n2(klev)
          ii=i*levfact
          jp=j+1
          jm=j-1
          if(j.eq.0)jm=n2(klev)-1
          if(j.eq.n2(klev))jp=1
          temp=(psif(i+1,j)-2.d0*psif(i,j)+psif(i-1,j))/hhr(klev)
     +        +(psif(i+1,j)-psif(i-1,j))*0.5d0/hr(klev)/rval(ii)
     +        +(psif(i,jp)-2.d0*psif(i,j)+psif(i,jm))/hhth(klev)
     +                                                      /rval(ii)**2
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
          temp=rf(i,j)-temp
          rc(i,j)=temp
          temp=dabs(temp)
          if(temp.gt.resmax)resmax=temp
   20   continue
c
c       if(klev.eq.1)write(6,100)klev,resmax
c       write(6,100)klev,resmax
  100   format(1x,'Level',i1,'  resmax=',f14.10)
c
  200 continue
c
  210 continue
      call trans1(psif,v,ibeg(klev),n1(klev),n2(klev))
      if(klev.eq.1.and.resmax.lt.tolcon)then
        nconv=1
      elseif(klev.eq.1.and.(resmax.lt.resmaxprev+1.d-7.and.
     +       resmax.gt.resmaxprev-1.d-7))then
        nconv=1
      endif
      if(klev.eq.1)resmaxprev=resmax
         
      return
      end
c==================================================
      subroutine corser(klev)
      implicit real*8(a-h,o-z)
      parameter(kr=256,kth=256,mgen=50000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kth),psif(0:kr,0:kth)
      real*8 rc(0:kr,0:kth),rf(0:kr,0:kth),tempa(0:kr,0:kth)
      real*8 hr(6),hth(6),hhr(6),hhth(6)
      real*8 rval(0:kr)
      integer ibeg(6),n1(6),n2(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hth,hhr,hhth,rval
      common /ints/ibeg,n1,n2
c
c------------------Start off with restriction
c
      call trans3(rc,rf,n1(klev),n2(klev))
c
      do 10 i=1,n1(klev+1)-1
      do 10 j=0,n2(klev+1)
        ii=i+i
        jj=j+j
        jjp=jj+1
        jjm=jj-1
        if(j.eq.0)jjm=n2(klev+1)-1
        if(j.eq.n2(klev+1))jjp=1
        temp=rf(ii,jj)
        temp=temp*2.d0+rf(ii+1,jj)+rf(ii-1,jj)
        temp=temp+rf(ii,jjp)+rf(ii,jjm)
        temp=temp*2.d0+rf(ii+1,jjp)+rf(ii-1,jjm)
        temp=temp+rf(ii+1,jjm)+rf(ii-1,jjp)
        temp=temp/16.d0
        rc(i,j)=temp
        psic(i,j)=0.d0
   10 continue
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
      parameter(kr=256,kth=256,mgen=50000)
      real*8 v(mgen),f(mgen)
      real*8 psic(0:kr,0:kth),psif(0:kr,0:kth)
      real*8 rc(0:kr,0:kth),rf(0:kr,0:kth),tempa(0:kr,0:kth)
      real*8 hr(6),hth(6),hhr(6),hhth(6)
      real*8 rval(0:kr)
      integer ibeg(6),n1(6),n2(6)
      common /backup/v,f
      common /arrays/psic,psif,rc,rf,tempa
      common /step/hr,hth,hhr,hhth,rval
      common /ints/ibeg,n1,n2
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
      parameter(kr=256,kth=256,mgen=50000)
      real*8 a2d(0:kr,0:kth),a1d(mgen)
c
      do 10 j=0,n2
      do 10 i=0,n1
        ii=j*(n1+1) + i + i1
        if(ii.gt.50000)then
          print*,'FATAL ERROR',ii
c          stop 
c          end
        endif
        a1d(ii)=a2d(i,j)
   10 continue
      return
      end
c==================================================
      subroutine trans2(a1d,a2d,i1,n1,n2)
      implicit real*8(a-h,o-z)
      parameter(kr=256,kth=256,mgen=50000)
      real*8 a2d(0:kr,0:kth),a1d(mgen)
c
      do 10 j=0,n2
      do 10 i=0,n1
        ii=j*(n1+1) + i + i1
        if(ii.gt.50000)then
          print*,'FATAL ERROR',ii
c          stop 
c          end
        endif
        a2d(i,j)=a1d(ii)
   10 continue
      return
      end
c==================================================
      subroutine trans3(a2d1,a2d2,n1,n2)
      implicit real*8(a-h,o-z)
      parameter(kr=256,kth=256,mgen=50000)
      real*8 a2d1(0:kr,0:kth),a2d2(0:kr,0:kth)
c
      do 10 j=0,n2
      do 10 i=0,n1
        a2d2(i,j)=a2d1(i,j)
   10 continue
      return
      end
c==================================================
      subroutine rufcon(a,nr,nth)
      parameter(kr=256,kth=256,mgen=24000)
      implicit real*8(a-h,o-z)
      real*8 a(0:kr,0:kth)
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
      amax=woppa(a,nr,nth)
      amin=titch(a,nr,nth)
      print*,'arraymax,min=',amax,amin
      adiff=amax-amin
      if(adiff.ne.0.)adiff=1./adiff
      do 30 j=0,nth
      jj=nth-j
      do 10 i=0,nr
      x=1.5+10.*(a(i,jj)-amin)*adiff
      n=x
      c(i)=ch(n)
  10  continue
      m1=0
      m2=nr
      ifact=2
      if(m2.gt.130)ifact=2
      write(6,20)(c(k),k=m1,m2,ifact)
  20  format(1x,200a1)
  30  continue
      return
      end
c==================================================
      subroutine rufcon2(a,b,nr,nth)
      parameter(kr=256,kth=256,mgen=24000)
      implicit real*8(a-h,o-z)
      real*8 a(0:kr,0:kth),b(0:kr,0:kth)
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
      amax=woppa(a,nr,nth)
      amin=titch(a,nr,nth)
      bmax=woppa(b,nr,nth)
      bmin=titch(b,nr,nth)
      print*,'arraymax,min=',amax,amin
      print*,'arraymax,min=',bmax,bmin
      adiff=amax-amin
      bdiff=bmax-bmin
      if(adiff.ne.0.)adiff=1./adiff
      if(bdiff.ne.0.)bdiff=1./bdiff
      do 30 j=0,nth
      jj=nth-j 
      do 10 i=0,nr
      x=1.5+10.*(a(i,jj)-amin)*adiff
      y=1.5+10.*(b(i,jj)-bmin)*bdiff
      na=x
      nb=y
      c(i)=ch(na) 
      c(i+nr+2)=ch(nb)
  10  continue 
      c(nr+1)=' '
      m1=0
      m2=nr*2+2
      ifact=2
      if(m2.gt.130)ifact=1
      write(6,20)(c(k),k=m1,m2,ifact)
  20  format(1x,200a1)
  30  continue
      return
      end
c==================================================
      double precision function titch(a,nr,nth)
      parameter(kr=256,kth=256,mgen=24000)
      implicit real*8(a-h,o-z)
      real*8 a(0:kr,0:kth)
c
      z=a(0,0)
      do 10 i=0,nr
      do 10 j=0,nth
      if(z.gt.a(i,j))z=a(i,j)
  10  continue
      titch=z
      return
      end
c==================================================
      double precision function woppa(a,nr,nth)
      parameter(kr=256,kth=256,mgen=24000)
      implicit real*8(a-h,o-z)
      real*8 a(0:kr,0:kth)
c
      z=a(0,0)
      do 10 i=0,nr
      do 10 j=0,nth
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
c 
c 
c 
c 
c     _ _ _ _ _ _ _----|      | ----|   \ `  ' /
c     ]-I-I-I-I-[  ----|      |     |    |. ` |
c      \ `   '_/       |     / \    |    | /^\|
c       []  `__|       ^    / ^ \   ^    | |*||
c       |__   ,|      / \  / ^ ^`\ / \   | ===|
c    ___| ___ ,|__   / ^  /=_=_=_=\ ^ \  |, `_|
c    I_I__I_I__I_I  (====(_________)_^___|____|____
c    \-\--|-|--/-/  |     I  [ ]__I I_I__|____I_I_| 
c     |[] `    '|_  |_   _|`__  ._[  _-\--|-|--/-/ 
c    / \  [] ` .| |-| |-| |_| |_| |_| | []   [] |
c   <===>      .|-=-=-=-=-=-=-=-=-=-=-|        / \ 
c   ] []|` ` [] | .   _________   .   |-      <===>
c   <===>  `  ' ||||  |       |  |||  |  []   <===>
c    \_/     -- ||||  |       |  |||  | .  '   \_/
c   ./|' . . . .|||||/|_______|\|||| /|. . . . .|\_
c 
