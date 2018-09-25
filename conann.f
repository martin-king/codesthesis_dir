c
c-----simple nag interactive contouring program
c     -----------------------------------------
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
      parameter(nx=201,ny=201,ng=2020)
c
      dimension r(nx),th(ny),arr(nx,ny)
      dimension x(nx,ny),y(nx,ny)
      character*50 label,fname
c
      print*,'data filename'
      read*,fname
      open(12,file=fname)
      read(12,*)mx,my
      mx=mx+1
      my=my+1
      if(mx.gt.nx.or.my.gt.ny)then
        print*,'too many points'
        close(12)
        stop
      end if
      read(12,*)rmax,rmin,thmax,thmin
      read(12,*)((arr(i,j),i=1,mx),j=1,my)
      close(12)
c
      print*,'Read data file; file closed'
c
      dr=(rmax-rmin)/dble(mx-1)
      dth=(thmax-thmin)/dble(my-1)
      do 1 i=1,mx
        r(i)=rmin+dble(i-1)*dr
    1 continue
      do 2 j=1,my
        th(j)=thmin+dble(j-1)*dth
    2 continue
      xmax=r(1)*dcos(th(1))
      xmin=xmax
      ymax=r(1)*dsin(th(1))
      ymin=ymax
      do 3 i=1,mx
      do 3 j=1,my
        x(i,j)=r(i)*dcos(th(j))
        y(i,j)=r(i)*dsin(th(j))
        if(x(i,j).gt.xmax)xmax=x(i,j)
        if(y(i,j).gt.ymax)ymax=y(i,j)
        if(x(i,j).lt.xmin)xmin=x(i,j)
        if(y(i,j).lt.ymin)ymin=y(i,j)
    3 continue
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
c     read*,label
c     lablen = len(label)
c     call j06ajf(1,label,lablen)
c     print*,'y-axis label'
c     read*,label
c     lablen = len(label)
c     call j06ajf(2,label,lablen)
c
c-----add plot title
c
      print*,'Title of plot'
      read*,label
      call j06ahf(label)
c
c-----call contouring routine

      call contor(x,y,arr,mx,my,0)
c
c-----terminate nag
c
      call j06wzf
c
      stop
      end
c------------------------------------------------------------
c============================================================
      subroutine contor(x,y,f,mx,my)
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
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
      parameter(nx=201,ny=201,ng=2020)
      implicit real*8(a-h,o-z)
      real*8 f(nx,ny),x(nx,ny),y(nx,ny)
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
c                           assumes circular boundary
      call gramov(x(1,1),y(1,1))
      do 48 j=2,my
        call gralin(x(1,j),y(1,j))
   48 continue
      call gramov(x(mx,1),y(mx,1))
      do 49 j=2,my
        call gralin(x(mx,j),y(mx,j))
   49 continue
      if(dabs(x(mx,1)-x(mx,my)).gt.1.d-8)then
        call gramov(x(mx,1),y(mx,1))
        call gralin(x(1,1),y(1,1))
        call gramov(x(mx,my),y(mx,my))
        call gralin(x(1,my),y(1,my))
      end if
c
  50  continue
      print*,'Max=',fmax
      print*,'Min=',fmin
      print*,'Ok? (y/n) '
      read*,answer
      if(answer.eq.'n'.or.answer.eq.'N')then
          print*,'Max,min= '
          read*,fmax,fmin
      end if
   51 continue
      print*,'#/contours='
      read*,nc
      if(nc.gt.100)then
          print*,'nc>100'
          goto 51
      end if
      if(nc.le.0)return
      print*,'Automatic level selection (y/n)? '
      read*,answer
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
          do 90 i=1,mx-1
          do 90 j=1,my-1
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
                  call gramov(x(i,j),y(i,j))
                  call gralin(x(i+1,j),y(i+1,j))
                  call gralin(x(i+1,j+1),y(i+1,j+1))
                  call gralin(x(i,j+1),y(i,j+1))
                  call gralin(x(i,j),y(i,j))
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
                  gx(i)=x(ii,jj)+prop*(x(ii,jj+1)-x(ii,jj))
                  gy(i)=y(ii,jj)+prop*(y(ii,jj+1)-y(ii,jj))
              else
c                                         Horiz. or corner crossings.
                  if(ii.eq.mx)then
c                                               Top right corner.
                      gx(i)=x(mx,jj)
                      gy(i)=y(mx,jj)
                  else
                      c1=f(ii,jj)-cv(k)
                      c2=f(ii+1,jj)-cv(k)
                      if(c1.eq.0.d0.or.c1-c2.eq.0.d0)then
c                                                         corner pts.
                          gx(i)=x(ii,jj)
                          gy(i)=y(ii,jj)
                      else
c                                                  horiz crossing.
                          prop=c1/(c1-c2)
                          gx(i)=x(ii,jj)+prop*(x(ii+1,jj)-x(ii,jj))
                          gy(i)=y(ii,jj)+prop*(y(ii+1,jj)-y(ii,jj))
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
      print*,'More contours? (y/n)'
      read*,answer
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
      real*8 a(nx,ny)
      z=a(1,1)
      do 10 i=1,mx
      do 10 j=1,my
        if(z.gt.a(i,j))z=a(i,j)
   10 continue
      arrmin=z
      return
      end
c---------------------------------------------------------------
      double precision function arrmax(a,nx,ny,mx,my)
      implicit real*8(a-h,o-z)
      real*8 a(nx,ny)
      z=a(1,1)
      do 10 i=1,mx
      do 10 j=1,my
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
