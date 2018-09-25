      parameter(mx=50,mx2=52)
      implicit real*8(a-h,o-z)
      real*8 x(mx),a(mx,mx),ainv(mx2,mx),b(mx2,mx)
      integer order(mx)
c
      print*,'Number of points, reference point'
      read*,nx,nref
      print*,'Constant step(1), user defined(0), GEF(2)?'
      read*,nans
      if(nans.eq.1)then
        do 5 i=1,nx
          x(i)=dble(i)
   5    continue
      end if
      if(nans.eq.0)then
        print*,'values of x'
        read*,(x(i),i=1,nx)
      end if
      if(nans.eq.2)then
        h=1.d0
        x(1)=0.d0
        print*,'GEF='
        read*,gef
        do 7 i=2,nx
          x(i)=x(i-1)+h
          h=h*gef
   7    continue
      end if
c
      do 8 i=1,nx
        order(i)=nx-i+1
   8  continue
c
      do 10 j=1,nx
        temp=1.d0
        h=x(j)-x(nref)
        a(1,j)=1.d0
        ainv(1,j)=0.d0
        b(1,j)=0.d0
        do 10 i=2,nx+2
          if(j.eq.nref)then
            temp=0.d0
          else
            temp=temp*h/dble(i-1)
          end if
          if(i.le.nx)a(i,j)=temp
          ainv(i,j)=0.d0
          if(i.le.nx)then
            b(i,j)=0.d0
          else
            b(i,j)=-temp
          end if
   10 continue
c
      do 20 j=1,nx
        ainv(j,j)=1.d0
        b(j,j)=1.d0
   20 continue
c
c     print*,'----------------------------'
c     do 25 j=1,nx
c       write(6,40)(a(i,j),i=1,nx+2)
c  25 continue
c     print*,'----------------------------'
c     do 26 j=1,nx
c       write(6,40)(b(i,j),i=1,nx+2)
c  26 continue
c     print*,'----------------------------'
c
      call gausse(a,b,ainv,nx)
c
      do 28 j=1,nx
        if(dabs(ainv(nx+1,j)).lt.1.d-6)then
          if(dabs(ainv(nx+2,j)).lt.1.d-6)then
            order(j)=0
          else
            order(j)=order(j)+1
          end if
        end if
   28 continue
c
      print*,' '
      write(6,38)'x_n',(x(i),i=1,nx)
      write(6,43)('=',i=1,13*nx+37)
      write(6,42)nx,nx+1,('x_',i,i=1,nx)
      write(6,43)('=',i=1,13*nx+37)
      do 30 j=1,nx
        if(j.le.10)then
        write(6,39)ainv(nx+1,j),ainv(nx+2,j),order(j),'y^',j-1,
     +              (ainv(i,j),i=1,nx)
        else
        write(6,40)ainv(nx+1,j),ainv(nx+2,j),order(j),'y^',j-1,
     +              (ainv(i,j),i=1,nx)
        end if
   30 continue
   38 format(32x,'| 'a3,1h=,100(1x,f12.8))
   39 format(2(1x,f12.8),1x,i3,'  | ',a2,i1,1h=,100(1x,f12.8))
   40 format(2(1x,f12.8),1x,i3,' | ',a2,i2,1h=,100(1x,f12.8))
   42 format(7x,2hy^,i2,9x,2hy^,i2,4x,'Ord | ',11x,100(a2,i2,9x))
   43 format(1x,1000a1)
c
      stop
      end
c--------------------------------------------------
      subroutine gausse(a,b,x,n)
      parameter(mx=50,mx2=52)
      implicit real*8(a-h,o-z)
      real*8 a(mx,mx),b(mx2,mx),x(mx2,mx)
c
c     Nb This code assumes that
c
c                 a11 a21 a31 ....
c                 a12 a22 a32 ....
c          A =    a13 a23 a33 ....
c                  :   :   :
c                  :   :   :
c
      if(n.eq.1)then
        x(1,1)=b(1,1)/a(1,1)
        x(2,1)=b(2,1)/a(1,1)
        x(3,1)=b(3,1)/a(1,1)
        return
      end if
c
      do 10 k=1,(n-1)
        kpvt=k
        kp1=k+1
        do 35 i=kp1,n
          if (abs(a(k,kpvt)).lt.abs(a(k,i))) kpvt=i
 35     continue
        if (abs(a(k,kpvt)).eq.0.d0)then
          print*,'Sorry, the matrix is singular, i.e. det=0.'
          print*,'Gaussian Elimination does not work.'
          stop
        end if
        if (kpvt.eq.k) goto 25
        do 45 j=k,n
          save=a(j,k)
          a(j,k)=a(j,kpvt)
          a(j,kpvt)=save
 45     continue
        do 46 j=1,n+2
          save=b(j,k)
          b(j,k)=b(j,kpvt)
          b(j,kpvt)=save
 46     continue
c
 25     continue
c
        do 15 i=k+1,n
          q=-a(k,i)/a(k,k)
          a(k,i)=0.d0
          do 20 j=k+1,n
            a(j,i)=a(j,i)+q*a(j,k)
 20       continue
          do 21 j=1,n+2
            b(j,i)=b(j,i)+q*b(j,k)
 21       continue
 15     continue
 10   continue
      t=1.d0
      do 70 i=1,n
        t=t*a(i,i)
   70 continue
      write(6,203)t
  203 format('determinant of matrix of coefficients =',e11.4)
c     print*,'----------------------------'
c     do 26 j=1,n
c       write(6,27)(a(i,j),i=1,n)
c  26 continue
c  27 format(100(1x,f10.6))
c     print*,'----------------------------'
      do 71 j=1,n+2
        x(j,n)=b(j,n)/a(n,n)
  71  continue
      do 75 i=n-1,1,-1
      do 75 k=1,n+2
        s=0.d0
        do 80 j=i+1,n
          s=s+a(j,i)*x(k,j)
 80     continue
        x(k,i)=(b(k,i)-s)/a(i,i)
 75   continue
      return
      end
