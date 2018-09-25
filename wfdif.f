c23456789012345678901234567890123456789012345678901234567890123456789012
c     Program wfdif.f
c     Martin P. King, University of Bath, 4 Feb 2000.
c----------------------------------------------------------------------c   
c     Using Langrange formula to fit a curve at ANY x. With 
c     Taylor's series to obtain recursion relations for the weights.  
c     Users provide nx gridpoints
c     program returns C_i,k, weights for each y_i(of nx) for each 
c     k(from 0 to nx-1) derivative of y.     
c     Ref: Bengt Fornberg,1998,Calculation of weights in FD formulas,
c          SIAM REV,Vol.40 No.3,pp.685-691.
c     obtainable fr bath IP address:
c       http://epubs.siam.org/sam-bin/dbq/toc/SIREV/40/3
c     Possible development:   
c     Oder of accuracy. and higher order errors.
c     Calculation of derivatives can be performed with y_i's supplied.  
c----------------------------------------------------------------------c
      parameter(mx=49)
c     you could increase mx if necessary
      implicit real*8(a-h,o-z)
      real*8 x(0:mx),C(0:mx,0:mx)
c
    
      print*,'This program gives you C_i, the coefficients for y_i,'  
      print*,'in a finite difference formula'
      print*,'---------------------------------------------------------'
c
      do 1 ii=1,100
      print*,'Number of points(max. 50), x of reference point(begins 0)'
      read*,nx,xx
      print*,'Constant step(1)(x begins at 0), user defined(2), GEF(3)?'
      read*,nans
c
      if(nans.eq.1)then
        do 5 i=0,nx-1
          x(i)=dble(i)
   5    continue
      end if
      if(nans.eq.2)then
        print*,'values of x'
        read*,(x(i),i=0,nx-1)
      end if
      if(nans.eq.3)then
        h=1.d0
        x(0)=0.d0
        print*,'GEF='
        read*,gef
        do 7 i=1,nx-1
          x(i)=x(i-1)+h
          h=h*gef
   7    continue
      end if
c    initialise all coefficients to zero
      do 10 k=0,mx
       do 10 i=0,mx
         C(i,k)=0.0d0
  10  continue
c
      C(0,0)=1.0d0
      c1=1.0d0
      c4=x(0)-xx
c    coefficient by coefficient
      do 50 i=1,nx-1    
        c2=1.0d0
        c5=c4
        c4=x(i)-xx
c    going up point by point until j=i-1
        do 40 j=0,i-1
          c3=x(i)-x(j)
          c2=c2*c3
c    checking for i=j
          if(j.eq.i-1)then
            do 20 k=i,1,-1
            C(i,k)=c1*(k*C(i-1,k-1)-c5*C(i-1,k))/c2
  20  continue
            C(i,0)=-c1*c5*C(i-1,0)/c2
          endif
c    going down derivative by derivative until k=1, ie. k=i,1,-1
          do 30 k=i,1,-1
            C(j,k)=(c4*C(j,k)-k*C(j,k-1))/c3
  30  continue
          C(j,0)=c4*C(j,0)/c3
  40  continue
      c1=c2
  50  continue
c
      print*,' '
      write(6,51)'x_n',(x(i),i=0,nx-1)
      write(6,54)('=',i=1,13*nx+37)
      write(6,53)nx,nx+1,('C_',i,i=0,nx-1)
      write(6,54)('=',i=1,13*nx+37)
      do 60 j=0,nx-1
        if(j.le.9)then
        write(6,55)0,0,0,'y^',j,
     +              (C(i,j),i=0,nx-1)
        else
        write(6,56)0,0,0,'y^',j,
     +              (C(i,j),i=0,nx-1)
        end if
   60 continue
   51 format(32x,'| 'a3,1h=,100(1x,f12.8))
   55 format(2(1x,f12.8)1x,i3,'  | ',a2,i1,1h=,100(1x,f12.8))
   56 format(2(1x,f12.8)1x,i3,' | ',a2,i2,1h=,100(1x,f12.8))
   53 format(7x,2hy^,i2,9x,2hy^,i2,4x,'Ord | ',11x,100(a2,i2,9x))
   54 format(1x,1000a1)   
    1 continue
c
      stop
      end
