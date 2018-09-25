c
c-----simple nag interactive contouring program
c     -----------------------------------------
c     modified to output staggered grid profile into file "fort.7"
c     and calcuate finite difference coefficients.
c
c-----program: contour.f
c
c
c-----nag routines used:
c
c         j06waf - initialise nag graphics system
c         j06wcf - set viewport
c         j06wbf - set data (axis) range
c         j06xff - set character set (hardware/software)
c         j06baf - plot points with symbols and/or straight lines
c         j06ccf - plot points using curves
c         j06aef - plot axes
c         j06ajf - add an axis title
c         j06ahf - add a plot title
c         j06wzf - close nag graphics system
c         j06ycf - draw to a point
c         j06yaf - move to a point
c
c----integer variables for plot routine j06baf
c
c         ifail     - nag error flag (set to zero)
c         istyl     - plot style: 0 - plot symbols only
c                                 1 - plot lines only
c                                 2 - plot symbols and lines
c         isym      - nag symbol number (1-8)
c
c----main program variables
c
c         label     - character string for axis labels
c         lablen    - character string length
c         margin    - leave margin for axis annotation (0 - no,1 - yes)
c
c         wx1,wy1   - nag screen position for window origin        
c         wx2,wy2   - nag screen position for window top right corner
c
c         x         - x data array dimensioned to nx
c         y         - y data array dimensioned to ny
c         arr       - function to be contoured; dimensions nx,ny
c
c         xmin,xmax - x axis limits
c         ymin,ymax - y axis limits
c
c
c-----nag implementation is double precision
c     --------------------------------------
c
      program contour
      implicit double precision(a-h,o-z)
c
c-----plot data arrays and other variables
c     ------------------------------------
c
      parameter(nx=200,ny=200,ng=2020,na=50,na2=52,nb=8)
c
      dimension x(0:nx,nb),y(0:ny,nb),arr(0:nx,0:ny)
      real*8 hx(nx,nb),hy(ny,nb)
      real*8 dx(nx),dy(ny)
      real*8 ainv(na2,na),cr(nx,3,nb),crr(nx,4,nb)
      real*8 cz(ny,3,nb),czz(ny,4,nb)
      real*8 r(0:nx),z(0:ny)
      integer order(na)
      integer ibeg(nb),n1(nb),n2(nb)
      character*50 label,fname
      character*3 ans
c
      call grid(dx,dy,mx,my,xmax,xmin,ymax,ymin)
c
      n1(1)=mx
      n2(1)=my
      do 5 i=1,mx
        hx(i,1)=dx(i)
    5 continue
c
      do 10 j=1,my
        hy(j,1)=dy(j)
   10 continue
c

      call level(hx,hy,x,y,n1,n2,ibeg,nlev,xmin,ymin)
c
      do 15 i=0,mx
      do 15 j=0,my
        arr(i,j)=0
   15 continue
c
      r(0)=xmin
      z(0)=ymin
      do 31 i=1,mx
        r(i)=r(i-1)+dx(i)
   31 continue
      do 32 j=1,my
        z(j)=z(j-1)+dy(j)
   32 continue
c
c-----initialise graphics system (screen or plotter)
c
      call xxxxxx
c
c-----initialise nag graphical supplement
c
      call j06waf
c

c------define viewport on the plotting surface (display position)
c
c     wx1 = 0.18d0
c     wx2 = 1.00d0
c     wy1 = 0.10d0
c     wy2 = 0.92d0
      yrange=ymax-ymin
      xrange=xmax-xmin
      if(xrange.gt.yrange)then
        wx1=0.18d0
        wx2=1.00d0
        factor=(1.d0-yrange/xrange)/2.d0
        wy1=0.1d0+0.82d0*factor
        wy2=0.92d0-factor*0.82d0
       else 
        factor=(1.d0-xrange/yrange)/2.d0
        wx1=0.18d0+0.82d0*factor
        wx2=1.d0-0.82d0*factor
        wy1=0.10d0
        wy2=0.92d0
      end if
c
      call j06wcf(wx1,wx2,wy1,wy2)
c
c-----map data region to viewport (set axis extrema)
c
      margin = 1
      call j06wbf(xmin,xmax,ymin,ymax,margin)
c
c-----use high quality software characters
c
      call j06xff(2)
c
c-----add axes 
c
c     call j06aef
c
c-----add axis labels (function "len" returns string length)
c
c     print*,'x-axis label'
c     label='R'
c     lablen = len(label)
c     call j06ajf(1,label,lablen)
c     print*,'y-axis label'
c     label='Z'
c     lablen = len(label)
c     call j06ajf(2,label,lablen)
c
c-----add plot title
c
c     print*,'Title of plot'
      label='GRID'
      call j06ahf(label)
c
c-----call contouring routine

      call contor(r,z,arr,mx,my)
c
c-----terminate nag
c
      call j06wzf
c
      print*,'Create coefficient file (y/n): '
      read*,ans
      if(ans.eq.'n'.or.ans.eq.'no') stop
      print*,'data filename'
      read*,fname
c
      do 35 l=1,nlev
        do 40 i=1,n1(l)/2
          call fdif(3,2,0,i,l,x,ainv,order)
          cr(i,1,l)=ainv(1,2)
          cr(i,2,l)=ainv(2,2)
          cr(i,3,l)=ainv(3,2)
          call fdif(4,2,0,i,l,x,ainv,order)
          crr(i,1,l)=ainv(1,3)
          crr(i,2,l)=ainv(2,3)
          crr(i,3,l)=ainv(3,3)
          crr(i,4,l)=ainv(4,3)
   40   continue
c
        do 45 i=n1(l)/2+1,n1(l)-1
          call fdif(3,2,0,i,l,x,ainv,order) 
          cr(i,1,l)=ainv(1,2) 
          cr(i,2,l)=ainv(2,2) 
          cr(i,3,l)=ainv(3,2) 
          call fdif(4,3,0,i,l,x,ainv,order) 
          crr(i,1,l)=ainv(1,3) 
          crr(i,2,l)=ainv(2,3) 
          crr(i,3,l)=ainv(3,3) 
          crr(i,4,l)=ainv(4,3)
   45   continue
c
        do 50 j=1,n2(l)/2
          call fdif(3,2,0,j,l,y,ainv,order) 
          cz(j,1,l)=ainv(1,2) 
          cz(j,2,l)=ainv(2,2) 
          cz(j,3,l)=ainv(3,2) 
          call fdif(4,2,0,j,l,y,ainv,order) 
          czz(j,1,l)=ainv(1,3) 
          czz(j,2,l)=ainv(2,3) 
          czz(j,3,l)=ainv(3,3) 
          czz(j,4,l)=ainv(4,3)
   50   continue
c 
        do 55 j=n2(l)/2+1,n2(l)-1
          call fdif(3,2,0,j,l,y,ainv,order) 
          cz(j,1,l)=ainv(1,2) 
          cz(j,2,l)=ainv(2,2) 
          cz(j,3,l)=ainv(3,2)  
          call fdif(4,3,0,j,l,y,ainv,order)  
          czz(j,1,l)=ainv(1,3)  
          czz(j,2,l)=ainv(2,3)  
          czz(j,3,l)=ainv(3,3)  
          czz(j,4,l)=ainv(4,3) 
   55   continue 
   35 continue
c
      call fdif(3,1,0,0,1,x,ainv,order)
      cr(0,1,1)=ainv(1,2)
      cr(0,2,1)=ainv(2,2)
      cr(0,3,1)=ainv(3,2)
      call fdif(3,3,0,n1(1),1,x,ainv,order)
      cr(n1(1),1,1)=ainv(1,2)
      cr(n1(1),2,1)=ainv(2,2)
      cr(n1(1),3,1)=ainv(3,2)
      call fdif(3,1,0,0,1,y,ainv,order)
      cz(0,1,1)=ainv(1,2)
      cz(0,2,1)=ainv(2,2)
      cz(0,3,1)=ainv(3,2)
      call fdif(3,3,0,n2(1),1,y,ainv,order)
      cz(n2(1),1,1)=ainv(1,2)
      cz(n2(1),2,1)=ainv(2,2)
      cz(n2(1),3,1)=ainv(3,2)
c
      open(7,file=fname)
      write(7,*)xmax,xmin,ymax,ymin
      write(7,*)nlev
      write(7,*)(ibeg(i),i=1,nlev)
      do 60 l=1,nlev
        write(7,*)n1(l),n2(l)
        write(7,*)(hx(i,l),i=1,n1(l))
        write(7,*)(hy(j,l),j=1,n2(l))
        write(7,*)((cr(i,k,l),k=1,3),i=1,n1(l)/2)
        write(7,*)((crr(i,k,l),k=1,4),i=1,n1(l)/2)
        write(7,*)((cr(i,k,l),k=1,3),i=n1(l)/2+1,n1(l)-1)
        write(7,*)((crr(i,k,l),k=1,4),i=n1(l)/2+1,n1(l)-1)
        write(7,*)((cz(j,k,l),k=1,3),j=1,n2(l)/2)
        write(7,*)((czz(j,k,l),k=1,4),j=1,n2(l)/2)
        write(7,*)((cz(j,k,l),k=1,3),j=n2(l)/2+1,n2(l)-1)
        write(7,*)((czz(j,k,l),k=1,4),j=n2(l)/2+1,n2(l)-1)
   60 continue
        write(7,*)(cr(0,k,1),k=1,3)
        write(7,*)(cr(n1(1),k,1),k=1,3)
        write(7,*)(cz(0,k,1),k=1,3)
        write(7,*)(cz(n2(1),k,1),k=1,3)
      close(12)
c
      stop
      end
c
c------------------------------------------------------------
c============================================================
      subroutine contor(x,y,f,mx,my)
c
c
c      *--------------------------------*
c      *                                *
c      *    Contour Plotting Routine    *
c      *                                *
c      *    c  D.A.S.Rees.   1982       *
c      *                                *
c      *                                *
c      * You are welcome to use this    *
c      * subroutine freely but I would  *
c      * require all users to acknowledge *
c      * its authorship in all reports  *
c      * papers and pulications of any  *
c      * kind
c      *--------------------------------*
c
c
      parameter(nx=200,ny=200,ng=2020)
      implicit real*8(a-h,o-z)
      real*8 f(0:nx,0:ny),x(0:nx),y(0:ny)
      real*8 gx(ng),gy(ng),cv(100)
      integer iw1(ng),iw2(ng),iord(ng),itest(4),jtest(4)
      logical l1,l2
      character*20 answer
c
c      ng=hoped for maximum number of points on a contour
c      n.b. ng>10*max(nx,ny) just to be on the safe side
c      should ng be too small the code will inform the user
c
      mxy=mx*my
      fmax=arrmax(f,nx,ny,mx,my)
      fmin=arrmin(f,nx,ny,mx,my)
      nmaxi=ng
c----                         Plot boundary
      call gramov(x(0),y(0))
      call gralin(x(0),y(my))
      call gralin(x(mx),y(my))
      call gralin(x(mx),y(0))
      call gralin(x(0),y(0))
c
  50  continue
c      print*,'Max=',fmax
c      print*,'Min=',fmin
c      print*,'Ok? (y/n) '
c      read*,answer
c      if(answer.eq.'n'.or.answer.eq.'N')then
c          print*,'Max,min= '
c          read*,fmax,fmin
c      end if
   51 continue
c      print*,'#/contours='
c      read*,nc
      nc=1
      if(nc.gt.100)then
          print*,'nc>100'
          goto 51
      end if
      if(nc.le.0)return
c      print*,'Automatic level selection (y/n)? '
c      read*,answer
      answer='y'
      if(answer.eq.'y'.or.answer.eq.'Y')then
          cint=(fmax-fmin)/dble(nc+1)
          do 60 i=1,nc
  60      cv(i)=fmin+dble(i)*cint
      else
          print*,'Type in the values='
          read(5,*)(cv(i),i=1,nc)
      end if
c
c--------------------Main loop.(For each contour value)
c
      do 210 k=1,nc
          write(6,70)k,cv(k)
   70     format(1x,'Value(',i2,')=',f12.7,3x,$)
          if((k/3)*3.eq.k)print*,' '
          nd=0
c                                            Searching each square.
          do 90 i=0,mx-1
          do 90 j=0,my-1
              na=0
              nb=0
              c1=f(i,j)-cv(k)
              c2=f(i+1,j)-cv(k)
              c3=f(i+1,j+1)-cv(k)
              c4=f(i,j+1)-cv(k)
              ic1=ising(c1)
              ic2=ising(c2)
              ic3=ising(c3)
              ic4=ising(c4)
c                                                 Along bottom.
              l1=ic1.gt.0.and.ic2.lt.0
              l2=ic1.lt.0.and.ic2.gt.0
              if(l1.or.l2)then
                  na=na+1
                  itest(na)=nmaxi*i+j
              end if
              if(ic1.eq.0)then
                  nb=nb+1
                  jtest(nb)=nmaxi*i+j
              end if
c                                                  Up right side.
              l1=ic2.gt.0.and.ic3.lt.0
              l2=ic2.lt.0.and.ic3.gt.0
              if(l1.or.l2)then
                  na=na+1
                  itest(na)=-(nmaxi*i+nmaxi+j)
              end if
              if(ic2.eq.0)then
                  nb=nb+1
                  jtest(nb)=nmaxi*(i+1)+j
              end if
c                                                 Along top.
              l1=ic3.gt.0.and.ic4.lt.0
              l2=ic3.lt.0.and.ic4.gt.0
              if(l1.or.l2)then
                  na=na+1
                  itest(na)=nmaxi*i+j+1
              end if
              if(ic3.eq.0)then
                  nb=nb+1
                  jtest(nb)=nmaxi*(i+1)+j+1
              end if
c                                                 Down left side.
              l1=ic4.gt.0.and.ic1.lt.0
              l2=ic4.lt.0.and.ic1.gt.0
              if(l1.or.l2)then
                  na=na+1
                  itest(na)=-(nmaxi*i+j)
              end if
              if(ic4.eq.0)then
                  nb=nb+1
                  jtest(nb)=nmaxi*i+j+1
              end if
c
c-----------------Test section.
c                                                          case 13
              if(nb.eq.4)then
                  call gramov(x(i),y(j))
                  call gralin(x(i+1),y(j))
                  call gralin(x(i+1),y(j+1))
                  call gralin(x(i),y(j+1))
                  call gralin(x(i),y(j))
                  goto 80
              end if
c                                                          case 12
              if(nb.eq.3)then
                  ii1=jtest(1)/nmaxi
                  jj1=jtest(1)-ii1*nmaxi
                  ii2=jtest(2)/nmaxi
                  jj2=jtest(2)-ii2*nmaxi
                  ii3=jtest(3)/nmaxi
                  jj3=jtest(3)-ii3*nmaxi
                  if(iabs(ii1-ii3).eq.1.and.iabs(jj1-jj3).eq.1)then
                      nd=nd+1
                      iw1(nd)=jtest(1)
                      iw2(nd)=jtest(2)
                      nd=nd+1
                      iw1(nd)=jtest(2)
                      iw2(nd)=jtest(3)
                  else
                      nd=nd+1
                      iw1(nd)=jtest(3)
                      iw2(nd)=jtest(1)
                      nd=nd+1
                      iw1(nd)=jtest(1)
                      iw2(nd)=jtest(2)
                  end if
                  goto 80
              end if
c                                                          cases 1 & 5
              if((na+nb).lt.2)goto 80
c                                                          cases 2 & 3
              if(na.eq.2.and.nb.eq.0)then
                  nd=nd+1
                  iw1(nd)=itest(1)
                  iw2(nd)=itest(2)
                  goto 80
              end if
c                                                          case 6
              if(na.eq.1.and.nb.eq.1)then
                  nd=nd+1
                  iw1(nd)=itest(1)
                  iw2(nd)=jtest(1)
                  goto 80
              end if
c                                                          cases 8,10,11
              if(na.eq.0.and.nb.eq.2)then
                  isum=ic1+ic2+ic3+ic4
                  isum=isum*isum
                  jsum=(ic2+ic4)*(ic3+ic1)
c                                                           case 10
                  if(isum.eq.4.and.jsum.eq.0)goto 80
c                                                           cases 8 & 11
                  nd=nd+1
                  iw1(nd)=jtest(1)
                  iw2(nd)=jtest(2)
                  goto 80
              end if
c                                                          case 7
              if(na.eq.2.and.nb.eq.1)then
                  nd=nd+1
                  iw1(nd)=jtest(1)
                  iw2(nd)=itest(1)
                  nd=nd+1
                  iw1(nd)=itest(1)
                  iw2(nd)=itest(2)
                  nd=nd+1
                  iw1(nd)=itest(2)
                  iw2(nd)=jtest(1)
                  goto 80
              end if
c                                                          case 9
              if(na.eq.1.and.nb.eq.2)then
                  nd=nd+1
                  iw1(nd)=jtest(1)
                  iw2(nd)=jtest(2)
                  nd=nd+1
                  iw1(nd)=jtest(2)
                  iw2(nd)=itest(1)
                  nd=nd+1
                  iw1(nd)=itest(1)
                  iw2(nd)=jtest(1)
                  goto 80
              end if
c                                                          case 4
              if(na.eq.4)then
                  nd=nd+1
                  iw1(nd)=itest(1)
                  iw2(nd)=itest(2)
                  nd=nd+1
                  iw1(nd)=itest(2)
                  iw2(nd)=itest(3)
                  nd=nd+1
                  iw1(nd)=itest(3)
                  iw2(nd)=itest(4)
                  nd=nd+1
                  iw1(nd)=itest(4)
                  iw2(nd)=itest(1)
                  goto 80
              end if
              print*,'Error in program logic.'
              print*,'i,j=',i,j
              print*,'k,cv(k)=',k,cv(k)
              print*,f(i,j+1),f(i+1,j+1)
              print*,f(i,j),f(i+1,j)
              print*,c4,c3
              print*,c1,c2
              print*,ic4,ic3
              print*,ic1,ic2
              stop
  80      continue
          if(nd.gt.ng)then
              write(6,999)nd,ng
  999         format(1x,//,'   *** ERROR ***',//,' na=',i4,'  ng=',i4)
              print*,'Increase the value of ng.'
              print*,'Setting nd=ng & continuing.... '
              print*,'Contour will NOT be complete'
              goto 91
          end if
  90      continue
  91      continue
          if(nd.eq.0)goto 210
c
c----------------------Ordering section.
c
  100     continue
          ntest=0
          na=2
          do 110 i=2,nd
              l1=iw1(1).eq.iw1(i)
              l2=iw1(1).eq.iw2(i)
              if(l1.or.l2)then
                 ntest=i
                  goto 120
              end if
  110     continue
  120     continue
          if(ntest.eq.0)then
              iord(1)=iw1(1)
              iord(2)=iw2(1)
          else
              iord(1)=iw2(1)
              iord(2)=iw1(1)
          end if
          iw1(1)=0
          iw2(1)=0
  130     continue
          do 140 i=1,nd
              l1=iord(na).eq.iw1(i)
              l2=iord(na).eq.iw2(i)
              if(l1)then
                  na=na+1
                  iord(na)=iw2(i)
                  iw1(i)=0
                  iw2(i)=0
                  goto 150
              end if
              if(l2)then
                  na=na+1
                  iord(na)=iw1(i)
                  iw1(i)=0
                  iw2(i)=0
                  goto 150
              end if
  140     continue
          goto 160
c........................goto 160 ==> no more points
  150     continue
          goto 130
  160     continue
          nagain=1
          nd=nd-na+1
          if(nd.eq.0)nagain=0
c........................i.e. not again
c
c--------------------Conversion to real coords. & plotting.
c
          do 170 i=1,na
              ii=iabs(iord(i)/nmaxi)
              jj=iabs(iord(i))-nmaxi*ii
              if(iord(i).lt.0)then
c                                              Upright crossings.
                  c1=f(ii,jj)-cv(k)
                  c2=f(ii,jj+1)-cv(k)
                  prop=c1/(c1-c2)
                  gx(i)=x(ii)
                  gy(i)=y(jj)+prop*(y(jj+1)-y(jj))
              else
c                                         Horiz. or corner crossings.
                  if(ii.eq.mx)then
c                                               Top right corner.
                      gx(i)=x(mx)
                      gy(i)=y(jj)
                  else
                      c1=f(ii,jj)-cv(k)
                      c2=f(ii+1,jj)-cv(k)
                      if(c1.eq.0.d0.or.c1-c2.eq.0.d0)then
c                                                         corner pts.
                          gx(i)=x(ii)
                          gy(i)=y(jj)
                      else
c                                                  horiz crossing.
                          prop=c1/(c1-c2)
                          gx(i)=x(ii)+prop*(x(ii+1)-x(ii))
                          gy(i)=y(jj)
                      end if
                  end if
              end if
  170     continue
          call grapol(gx,gy,na)
          if(nagain.eq.0)goto 210
          ntest=0
c
c---------------------Packing 'iord' for next bit of contour.
c
          do 190 i=1,ng
              if(iw1(i).ne.0)then
                  ntest=ntest+1
                  iw1(ntest)=iw1(i)
                  iw2(ntest)=iw2(i)
                  if(ntest.eq.nd)goto 200
              end if
  190     continue
  200     continue
c                                  Compaction completed.
          goto 100
  210 continue
c      print*,'More contours? (y/n)'
c      read*,answer
      answer='n'
      if(answer.eq.'y')goto 50
      print*,'Leaving Contour'
      return
      end
c---------------------------------------------------------------
c---------------------------------------------------------------
      function ising(a)
      implicit real*8(a-h,o-z)
      ising=0
      if(a.gt.0.)ising=1
      if(a.lt.0.)ising=-1
      end
c---------------------------------------------------------------
      double precision function arrmin(a,nx,ny,mx,my)
      implicit real*8(a-h,o-z)
      real*8 a(0:nx,0:ny)
      z=a(1,1)
      do 10 i=0,mx
      do 10 j=0,my
        if(z.gt.a(i,j))z=a(i,j)
   10 continue
      arrmin=z
      return
      end
c---------------------------------------------------------------
      double precision function arrmax(a,nx,ny,mx,my)
      implicit real*8(a-h,o-z)
      real*8 a(0:nx,0:ny)
      z=a(1,1)
      do 10 i=0,mx
      do 10 j=0,my
        if(z.lt.a(i,j))z=a(i,j)
   10 continue
      arrmax=z
      return
      end
c=================================================================
      subroutine grapol(x,y,n)
      implicit real*8(a-h,o-z)
      real*8 x(n),y(n)
c
      ifail=0
      call J06BAF(x,y,n,1,0,ifail)
      if(ifail.ne.0)then
        write(6,100)ifail
  100   format(1x,'IFAIL=',i3,' in GRAPOL')
      end if
      return
      end
c================================================================
      subroutine gracur(x,y,n)
      implicit real*8(a-h,o-z)
      real*8 x(n),y(n)
c
      ifail=10
      call J06CCF(x,y,n,1,ifail)
      if(ifail.ne.0)then
        write(6,100)ifail
  100   format(1x,'IFAIL=',i3,' in GRACUR')
      end if
      return
      end
c================================================================
      subroutine gramov(x,y)
      implicit real*8(a-h,o-z)
c
      call J06YAF(x,y)
      return
      end
c================================================================
      subroutine gralin(x,y)
      implicit real*8(a-h,o-z)
c
      call J06YCF(x,y)
      return
      end
c================================================================
      subroutine grid(dr,dz,nr,nz,rmax,r0,zmax,zmin)
c
      implicit real*8(a-h,o-z)
      parameter(nrmax=200,nzmax=200)
      real*8 dr(nrmax),dz(nzmax)
      character*4 ans
c
      print*,'number of intervals in r-direction: '
      read*,nr
      print*,'number of intervals in z-direction: '
      read*,nz
      print*,'r0,rmax,zmax,zmin'
      read*,r0,rmax,zmax,zmin
c
      print*,'staggered grid in r-direction (y/n): '
      read*,ans
      if(ans.eq.'y')then
        print*,'r distance for GEF application: '
        read*,rdist
        print*,'GEF in r direction: '
        read*,rgef
        print*,'first dr: '
        read*,dr(1)
c
        drt=dr(1)
      
        do 10 i=2,nr/2
          dr(i)=dr(i-1)*rgef
          drt=drt+dr(i)
          if(drt.gt.rdist)then
            imax=i
            drt=drt-dr(i)
            goto 20
          endif
   10   continue
        print*,'too many grid points in GEF region'
        stop
   20   continue
c
        if(2*imax-2.gt.nr)then
          print*,'grid points exceeded'
          stop
        endif
c
        rover=rmax-r0-2.d0*drt
c
        do 30 i=imax,nr-imax+1
          dr(i)=rover/dble(nr-2*imax+2)
   30   continue
c
        do 40 i=nr-imax+1,nr
          dr(i)=dr(nr-i+1)
   40   continue
c
      else
        imax=0
        do 50 i=1,nr
          dr(i)=(rmax-r0)/dble(nr)
   50   continue
c
      endif
c
      print*,'staggered grid in z-direction (y/n): '
      read*,ans
      if(ans.eq.'y')then
        print*,'z distance for GEF application: '
        read*,zdist
        print*,'GEF in z direction: '
        read*,zgef
        print*,'first dz: '
        read*,dz(1)
c
        dzt=dz(1)
      
        do 60 j=2,nz/2
          dz(j)=dz(j-1)*zgef
          dzt=dzt+dz(j)
          if(dzt.gt.zdist)then
            jmax=j
            dzt=dzt-dz(j)
            goto 70
          endif
   60   continue
        print*,'too many grid points in GEF region'
        stop
   70   continue
c
        if(2*jmax-2.gt.nz)then 
          print*,'grid points exceeded' 
          stop 
        endif
c
        zover=zmax-zmin-2.d0*dzt 
c 
        do 80 j=jmax,nz-jmax+1 
          dz(j)=zover/dble(nz-2*jmax+2) 
   80   continue 
c 
        do 90 j=nz-jmax+1,nz 
          dz(j)=dz(nz-j+1) 
   90   continue 
c 
      else 
        jmax=0
        do 100 j=1,nz 
          dz(j)=(zmax-zmin)/dble(nz) 
  100   continue 
c 
      endif 
c
      if(nr.gt.nz)then
        nmax=nr
      else
        nmax=nz
      endif
c
      do 110 i=1,nmax
        print*,dr(i),dz(i)
  110 continue
c
      rsum=0.d0
      do 120 i=1,nr
        rsum=rsum+dr(i)
  120 continue
c
      zsum=0.d0
      do 130 j=1,nz
        zsum=zsum+dz(j)
  130 continue
c
      print*,'check rgap = ',rsum
      print*,'check zgap = ',zsum
c 
      return
      end
c------------------------------------------------------------
      subroutine level(hr,hz,r,z,n1,n2,ibeg,nlev,r0,zmin)
      implicit real*8(a-h,o-z)
      parameter(kr=200,kz=200,nb=8)
      character*12 ans
      real*8 hr(kr,nb),hz(kz,nb),r(0:kr,nb),z(0:kz,nb)
      integer ibeg(nb),n1(nb),n2(nb)
c
      nlev=0
      itest=n1(1)
      itest2=n1(1)
      jtest=n2(1)
      jtest2=n2(1)
      print*,itest,itest2,jtest,jtest2
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
      if(nlev.gt.6)nlev=6
      print*,'Number of levels=',nlev
      print*,'Change? (y/n)'
      read*,ans
       if(ans.eq.'y'.or.ans.eq.'yes')then
         print*,'nlev'
         read*,nlev
       endif
c
      nvs=10000
      nrel1=2
      nrel2=2
c
      do 15 l=1,nlev
        r(0,l)=r0
        n1(l)=n1(1)*2/(2**l)
      do 15 i=1,n1(l)
        if(l.ne.1)then
          hr(i,l)=hr(2*i,l-1)+hr(2*i-1,l-1)
        endif
        r(i,l)=r(i-1,l)+hr(i,l)
   15 continue
c
      do 20 l=1,nlev
        n2(l)=n2(1)*2/(2**l)
        z(0,l)=zmin
      do 20 j=1,n2(l)
        if(l.ne.1)then
          hz(j,l)=hz(2*j,l-1)+hz(2*j-1,l-1)
        endif
        z(j,l)=z(j-1,l)+hz(j,l)
   20 continue
c
      ii=1
      do 25 l=1,nlev
        ibeg(l)=ii
        ii=ii+(n1(l)+1)*(n2(l)+1)+1
   25 continue
c
      return
      end
c============================================================
      subroutine fdif(nx,nref,nans,k,l,r,ainv,order)
      parameter(mx=50,mx2=52,na=200,nb=8)
      implicit real*8(a-h,o-z)
      real*8 x(mx),r(0:na,nb),a(mx,mx),ainv(mx2,mx),b(mx2,mx)
      integer order(mx)
c
c     print*,'Number of points, reference point'
c     read*,nx,nref
c     print*,'Constant step(1), user defined(0), GEF(2)?'
c     read*,nans
      if(nans.eq.1)then
        do 5 i=1,nx
          x(i)=dble(i)
   5    continue
      end if
      if(nans.eq.0)then
c       print*,'values of x'
c       read*,(x(i),i=1,nx)
        do 6 i=1,nx
          x(i)=r(k-nref+i,l)
   6    continue
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
c     print*,' '
c     write(6,38)'x_n',(x(i),i=1,nx)
c     write(6,43)('=',i=1,13*nx+37)
c     write(6,42)nx,nx+1,('x_',i,i=1,nx)
c     write(6,43)('=',i=1,13*nx+37)
c     do 30 j=1,nx
c       if(j.le.10)then
c       write(6,39)ainv(nx+1,j),ainv(nx+2,j),order(j),'y^',j-1,
c    +              (ainv(i,j),i=1,nx)
c       else
c       write(6,40)ainv(nx+1,j),ainv(nx+2,j),order(j),'y^',j-1,
c    +              (ainv(i,j),i=1,nx)
c       end if
c  30 continue
c  38 format(32x,'| 'a3,1h=,100(1x,f12.8))
   39 format(2(1x,f12.8),1x,i3,'  | ',a2,i1,1h=,100(1x,f12.8))
   40 format(2(1x,f12.8),1x,i3,' | ',a2,i2,1h=,100(1x,f12.8))
c  42 format(7x,2hy^,i2,9x,2hy^,i2,4x,'Ord | ',11x,100(a2,i2,9x))
c  43 format(1x,1000a1)
c
      return
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
c     write(6,203)t
c 203 format('determinant of matrix of coefficients =',e11.4)
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

