c given a field on a mesh, interpolate to a finer (in both r and p-dir) mesh 
c just a simple 2-points linear interpolation
c Martin King, ICTP, 11-4-05
      
      implicit real*8(a-h,o-z)
      parameter(kr=400,kp=4000)
      real*8 A(0:kr,1,0:kp),B(0:kr,1,0:kp),C(0:kr,1,0:kp)
      real*8 uel(0:kr,1,0:kp),vel(0:kr,1,0:kp),wel(0:kr,1,0:kp)
      real*8 th(0:kr,1,0:kp)
      real*8 An(0:kr,1,0:kp),Bn(0:kr,1,0:kp),Cn(0:kr,1,0:kp)
      real*8 ueln(0:kr,1,0:kp),veln(0:kr,1,0:kp),weln(0:kr,1,0:kp)
      real*8 thn(0:kr,1,0:kp)
      real*8 rold(0:kr),rnew(0:kr)
      real*8 pold(0:kp),pnew(0:kp)
      character*30 fname9
c      common/val/A,B,C,uel,vel,wel,th
c      common/valn/An,Bn,Cn,ueln,veln,weln,thn
c      common/num/nr,np,nrnew
c      common/radius/dr,drnew,r0,rmax,rold,rnew

c---------- some necessary inputs

      print*,'number of intervals in r-direction: '
      read*,nr
      nz=1
      print*,'number of intervals in phi-direction: '
      read*,np
      print*,'r0,rmax'
      read*,r0,rmax

      print*,'number of new intervals in r-direction: '
      read*,nrnew

      print*,'number of new intervals in p-direction: '
      read*,npnew

      dr=(rmax-r0)/dble(nr)
      drnew=(rmax-r0)/dble(nrnew)
      dp=1./dble(np)
      dpnew=1./dble(npnew)

c--------------- reading files
        print*,'filename for vorticity A: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)(((A(i,j,k),i=0,nr+1),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for vorticity B: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)(((B(i,j,k),i=0,nr),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for vorticity C: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)(((C(i,j,k),i=0,nr),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for velocity u: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)(((uel(i,j,k),i=0,nr),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for velocity v: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)(((vel(i,j,k),i=0,nr+1),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for velocity w: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)(((wel(i,j,k),i=0,nr+1),j=1,1),k=0,np)
        close(11)
c
        print*,'filename for temperature: '
        read*,fname9
        open(11,file=fname9,form='unformatted')
        read(11)(((th(i,j,k),i=0,nr+1),j=1,1),k=0,np+1)
        close(11)        

c reading end

c------------ defining coordinate of old points
        do  i=0,nr+1
         if(i.eq.0)then
           rold(i)=r0
         elseif(i.eq.nr+1)then
           rold(i)=rmax
         else
           rold(i)=r0+(dble(i)-0.5d0)*dr
         endif
        enddo

c defining coordinate of new points
        do  i=0,nrnew+1
         if(i.eq.0)then
           rnew(i)=r0
         elseif(i.eq.nrnew+1)then
           rnew(i)=rmax
         else
           rnew(i)=r0+(dble(i)-0.5d0)*drnew
         endif
        enddo

c interpolating for A, vel, wel, th in radial direction first  
        do k=0,np    
        do in=0,nrnew+1
          if(in.eq.0)then 
            An(in,1,k)= A(0,1,k)
            veln(in,1,k)= vel(0,1,k)
            weln(in,1,k)= wel(0,1,k)
            thn(in,1,k)= th(0,1,k)
          elseif(in.eq.nrnew+1)then
            An(in,1,k)= A(nr+1,1,k)
            veln(in,1,k)= vel(nr+1,1,k)
            weln(in,1,k)= wel(nr+1,1,k)
            thn(in,1,k)= th(nr+1,1,k)
          else
            do i=0,nr
             if (rnew(in).ge.rold(i) .and. rnew(in).lt.rold(i+1))then 
              An(in,1,k)=(A(i,1,k)+A(i+1,1,k))/2.d0
              veln(in,1,k)=(vel(i,1,k)+vel(i+1,1,k))/2.d0
              weln(in,1,k)=(wel(i,1,k)+wel(i+1,1,k))/2.d0
              thn(in,1,k)=(th(i,1,k)+th(i+1,1,k))/2.d0
cmpk 30 mar 2009 added for velvor2dhorconvtank code
              thn(in,1,np+1)=(th(i,1,np+1)+th(i+1,1,np+1))/2.d0
c
             endif
            enddo
           endif
c           if(k.eq.240)print*,thn(in,1,k)
         enddo     
         enddo

c saving back to A, vel, wel and th
         do i=0,nrnew+1
          do k=0,np
            A(i,1,k)=An(i,1,k)
            vel(i,1,k)=veln(i,1,k)
            wel(i,1,k)=weln(i,1,k)
            th(i,1,k)=thn(i,1,k)
cmpk 30 mar 2009
            th(i,1,np+1)=thn(i,1,np+1)
          enddo
         enddo
         
         pold(0)=0.d0
         do k=1,np
           pold(k)=pold(k-1)+dp
         enddo

         pnew(0)=0.d0
         do k=1,npnew
            pnew(k)=pnew(k-1)+dpnew
         enddo

c interpolating in tangential direction
         do i=0,nrnew+1
          do k=0,npnew
            if(k.eq.0)then
              An(i,1,k)=A(i,1,0)
              veln(i,1,k)=vel(i,1,0)
              weln(i,1,k)=wel(i,1,0)
              thn(i,1,k)=th(i,1,0)
            elseif(k.eq.npnew)then
              An(i,1,k)=A(i,1,np)
              veln(i,1,k)=vel(i,1,np)
              weln(i,1,k)=wel(i,1,np)
              thn(i,1,k)=th(i,1,np)
            else
             do kk=0,np
              if(pnew(k).ge.pold(kk) .and. pnew(k).lt.pold(kk+1))then
               An(i,1,k)=(A(i,1,kk)+A(i,1,kk+1))*0.5d0
               veln(i,1,k)=(vel(i,1,kk)+vel(i,1,kk+1))*0.5d0
               weln(i,1,k)=(wel(i,1,kk)+wel(i,1,kk+1))*0.5d0
               thn(i,1,k)=(th(i,1,kk)+th(i,1,kk+1))*0.5d0
              endif
             enddo
            endif
           enddo
          enddo

c--------------------- interpolating for B, C, uel,
c------------ defining coordinate of old points
        do  i=0,nr
         if(i.eq.0)then
           rold(i)=r0
         elseif(i.eq.nr)then
           rold(i)=rmax
         else
           rold(i)=r0+(dble(i))*dr
         endif
        enddo
c defining coordinate of new points
        do  i=0,nrnew
         if(i.eq.0)then
           rnew(i)=r0
         elseif(i.eq.nrnew)then
           rnew(i)=rmax
         else
           rnew(i)=r0+dble(i)*drnew
         endif
        enddo   

c interpolating in radial direction
        do k=0,np    
        do in=0,nrnew
          if(in.eq.0)then 
            Bn(in,1,k)= B(0,1,k)
            Cn(in,1,k)= C(0,1,k)
            ueln(in,1,k)= uel(0,1,k)
          elseif(in.eq.nrnew)then
            Bn(in,1,k)= B(nr,1,k)
            Cn(in,1,k)= C(nr,1,k)
            ueln(in,1,k)= uel(nr,1,k)
          else
            do i=0,nr-1
             if (rnew(in).ge.rold(i) .and. rnew(in).lt.rold(i+1))then 
              Bn(in,1,k)=(B(i,1,k)+B(i+1,1,k))/2.d0
              Cn(in,1,k)=(C(i,1,k)+C(i+1,1,k))/2.d0
              ueln(in,1,k)=(uel(i,1,k)+uel(i+1,1,k))/2.d0
             endif
            enddo
          endif
         enddo     
         enddo

c saving back to
         do i=0,nrnew
          do k=0,np
            B(i,1,k)=Bn(i,1,k)
            C(i,1,k)=Cn(i,1,k)
            uel(i,1,k)=ueln(i,1,k)
          enddo
         enddo

c interpolating in tangential direction
         do i=0,nrnew
          do k=0,npnew
            if(k.eq.0)then
              Bn(i,1,k)=B(i,1,0)
              Cn(i,1,k)=C(i,1,0)
              ueln(i,1,k)=uel(i,1,0)
            elseif(k.eq.npnew)then
              Bn(i,1,k)=B(i,1,np)
              Cn(i,1,k)=C(i,1,np)
              ueln(i,1,k)=uel(i,1,np)
            else
             do kk=0,np
              if(pnew(k).ge.pold(kk) .and. pnew(k).lt.pold(kk+1))then
               Bn(i,1,k)=(B(i,1,kk)+B(i,1,kk+1))*0.5d0
               Cn(i,1,k)=(C(i,1,kk)+C(i,1,kk+1))*0.5d0
               ueln(i,1,k)=(uel(i,1,kk)+uel(i,1,kk+1))*0.5d0
              endif
             enddo
            endif
           enddo
          enddo

         
c interpolating ends   

c writing 
       
       fname9='avortresfiner.dat'
       open(11,file=fname9,form='unformatted')
       write(11)(((An(i,j,k),i=0,nrnew+1),j=1,1),k=0,npnew)
       close(11)

       fname9='bvortresfiner.dat'
       open(11,file=fname9,form='unformatted')
       write(11)(((Bn(i,j,k),i=0,nrnew),j=1,1),k=0,npnew)
       close(11)

       fname9='cvortresfiner.dat'
       open(11,file=fname9,form='unformatted')
       write(11)(((Cn(i,j,k),i=0,nrnew),j=1,1),k=0,npnew)
       close(11)

       fname9='uelresfiner.dat'
       open(11,file=fname9,form='unformatted')
       write(11)(((ueln(i,j,k),i=0,nrnew),j=1,1),k=0,npnew)
       close(11)

       fname9='velresfiner.dat'
       open(11,file=fname9,form='unformatted')
       write(11)(((veln(i,j,k),i=0,nrnew+1),j=1,1),k=0,npnew)
       close(11)

       fname9='welresfiner.dat'
       open(11,file=fname9,form='unformatted')
       write(11)(((weln(i,j,k),i=0,nrnew+1),j=1,1),k=0,npnew)
       close(11)

       fname9='thresfiner.dat'
       open(11,file=fname9,form='unformatted')
       write(11)(((thn(i,j,k),i=0,nrnew+1),j=1,1),k=0,npnew+1)
       close(11)

  201  format(1(1x,e16.10))           
        stop
        end
 
          




