c
c-----simple nag graphics demonstration program
c     -----------------------------------------
c
c-----program: graffit.f
c
c     to compile and run on ss1, mungo etc. type
c                           ===  ====
c
c     f77 -o graffit.run graffit.f  -lnaggl04 -lnagaps -lnagfls
c
c     or copy f77g from my filespace and type
c
c     f77g graffit
c     
c
c    **** You may need to type chmod u+x f77 to get execute access ***
c                              =================
c
c-----nag routines used:
c
c         j06waf - initialise nag graphics system
c         j06wcf - set viewport
c         j06wbf - set data (axis) range
c         j06xff - set character set (hardware/software)
c         j06baf - plot points with symbols and/or straight lines
c         j06aef - plot axes
c         j06ajf - add an axis title
c         j06ahf - add a plot title
c         j06wzf - close nag graphics system
c
c----integer variables for plot routine j06baf
c
c         ifail     - nag error flag (set to zero)
c         istyl     - plot style: 0 - plot symbols only
c                                 1 - plot lines only
c                                 2 - plot symbols and lines
c         isym      - nag symbol number (1-8)
c         nrowmax      - dimension of arrays x,y
c         ncolmax      - second dimension of data2d
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
c         x,y       - x,y data arrays dimensioned to ndim
c
c         xmin,xmax - x axis limits
c         ymin,ymay - y ayis limits
c
c
c-----nag implementation is double precision
c     --------------------------------------
c
      program graffit
      implicit double precision(a-h,o-z)
c
c-----plot data arrays and other variables
c     ------------------------------------
c
      parameter(nrowmax=100000,ncolmax=50)
c
      dimension x(nrowmax),y(nrowmax),data2d(nrowmax,ncolmax)
      dimension test(nrowmax*ncolmax)
      character*30 label,fname,answer
      common /limits/xmax,xmin,ymax,ymin

c
c-----Setting up the graphics first.
c
c=============================================================
c=============================================================
c
c-----Open graphics file
c
      print*,'Filename for graphics output (use a .ps suffix)'
      read*,fname
      open(7,file=fname)
c
c-----initialise graphics system (screen or plotter)
c
      call xxxxxx
c
c-----initialise nag graphical supplement
c
      call j06waf
c
c-----define viewport on the plotting surface (display position)
c
      wx1 = 0.15d0
      wx2 = 1.00d0
      wy1 = 0.10d0
      wy2 = 0.92d0
c
      call j06wcf(wx1,wx2,wy1,wy2)
c
c-----map data region to viewport (set axis extrema)
c
      print*,'xmax,xmin,ymax,ymin for plotting'
      read*,xmax,xmin,ymax,ymin
c
      margin = 1
c
      call j06wbf(xmin,xmax,ymin,ymax,margin)
c
c----- size  of characters 1=normal
c
      scalex=1.d0
      scaley=1.d0
      print*,'Scale of text characters (x,y)'
      read*,scalex, scaley
      call j06xgf(scalex,scaley)
c
c-----use high quality software characters
c
c     call j06xff(2)
c
c-----mark points with symbols and join with straight lines
c       for further details see nag documentation
c
c     istyl = 2 ==> lines and symbols
c     istyl = 1 ==> lines
c     istyl = 0 ==> symbols
c      isymb takes values between 0 and 9 inclusive
      istyl=1
      isymb=0
      ifail = 0
c
c-----add axes
c
      print*,'Axes? (y/n)'
      read*,answer
      if(answer.eq.'y')call j06aef
c
c-----add axis labels (function "len" returns string length)
c
c
c     print*,'x-axis label....'
c     read*,label
c     label = 'x - axis      '
c     call j06ajf(1,label,lablen)
c
c     print*,'y-axis label....'
c     read*,label
c     label = 'y  axis      '
c     call j06ajf(2,label,lablen)
c
c-----add plot title
c
      print*,'plot title....'
      read*,label
      call j06ahf(label)
c
c=============================================================

c-----reading data
c     ------------
   10 continue

      print*,'Enter data filename (or `no`)....'
      read*,fname
      if(fname.eq.'no')goto 999
      open(12,file=fname)
c
c-----Check for number of rows in datafile
c
      do i=1,nrowmax+1
        npts=i-1
        read(12,*,err=888)dum
      end do
      print*,'Too many points'
      stop

  888 continue
      print*,'Number of rows=',npts
      rewind(12)
c
c-----Check for number of columns in datafile
c
      do j=1,ncolmax+1
        ncols=j
        read(12,*,err=887)(test(i),i=1,j*npts)
        read(12,*,err=889)dum
        rewind(12)
      end do
      print*,'Too many columns'
      print*,ncols,' columns'
      stop

  887 continue
      print*,'One column is no good, matey'
      print*,char(7)
      stop

  889 continue
      print*,'Number of columns=',ncols
      rewind(12)
c
c-----Finally read the datafile
c
      do i=1,npts
        read(12,*)(data2d(i,j),j=1,ncols)
      end do

      close(12)

      print*,'read data file'
c
c====================================================
c
c-----Plotting now
c
c     print*,'Line (1) or dashed (2)?'
c     read*,line
   20 continue
      line=1
      print*,'column for `y`-variable= '
      read*,icoly
      if(icoly.eq.0)goto 10
      print*,'column for `x`-variable='
      read*,icolx

      print*,'line(0)  dashed(-1)  symbol(1-9)'
      read*,iline
      ifail=0
      if(iline.eq.0)then
        line=1
        call j06yrf(line)
        istyl=1
        isymb=0
        call j06yrf(line)
        istyl=1
        isymb=0
      end if
      if(iline.gt.0)then
        istyl=0
        isymb=iline
      end if
      do 30 i=1,npts
        x(i)=data2d(i,icolx)
        y(i)=data2d(i,icoly)
   30 continue
c     call j06baf(x,y,npts,istyl,isymb,ifail)
      call myline(x,y,npts,istyl,isymb,ifail)
      goto 20
c====================================================

  999 continue
      call j06yrf(line)
c
c------------Axes to plot
c
      print*,'plot axes (x/y/b/n)'
      read*,fname
      if(fname.eq.'x'.or.fname.eq.'b')then
        x(1)=xmin
        x(2)=xmax
        y(1)=0.d0
        y(2)=0.d0
        call myline(x,y,2,istyl,isymb,ifail)
      end if
      if(fname.eq.'y'.or.fname.eq.'b')then
        y(1)=ymin
        y(2)=ymax
        x(1)=0.d0
        x(2)=0.d0
        call myline(x,y,2,istyl,isymb,ifail)
      end if
c
c-----terminate nag
c
    
      call j06wzf
      close(7)
      stop
      end
c============================================================
c============================================================
      subroutine myline(x,y,npts,istyl,isymb,ifail)
      implicit real*8(a-h,o-z)
c
c     Subroutine joins data points with straight lines but
c     performs data-clipping at the computational boundaries
c     ------------------------------------------------------
c
      parameter(nrowmax=100000,ncolmax=50)
c
      dimension x(nrowmax),y(nrowmax)
      dimension xp(2),yp(2)
      dimension xpp(2),ypp(2)
      common /limits/xmax,xmin,ymax,ymin
      
c
c     Clipping graphs. l1 and l2 are codes for where the
c     datapoints are relative to the plotting area
c
c              02 |  12  | 22
c              --------------
c              01 |  11  | 21
c              --------------
c              00 |  10  | 20
c
      do 1 i=2,npts
        xp(1)=x(i-1)
        xp(2)=x(i)
        yp(1)=y(i-1)
        yp(2)=y(i)
        lx1=0
        lx2=0
        ly1=0
        ly2=0
        if(xp(1).ge.xmin.and.xp(1).le.xmax)lx1=1
        if(xp(2).ge.xmin.and.xp(2).le.xmax)lx2=1
        if(yp(1).ge.ymin.and.yp(1).le.ymax)ly1=1
        if(yp(2).ge.ymin.and.yp(2).le.ymax)ly2=1
c
c------- Both points in plotting area
c
        if(lx1+lx2+ly1+ly2.eq.4)then
          call j06baf(xp,yp,2,istyl,isymb,ifail)
          goto 999
        end if
c
c------- Neither point in plotting area but line outside
c
        if(xp(1).lt.xmin.and.xp(2).lt.xmin)then
          goto 999
        end if
        if(xp(1).gt.xmax.and.xp(2).gt.xmax)then
          goto 999
        end if
        if(yp(1).lt.ymin.and.yp(2).lt.ymin)then
          goto 999
        end if
        if(yp(1).gt.ymax.and.yp(2).gt.ymax)then
          goto 999
        end if
c
c-------- One point in cases
c
        if(lx1*ly1+lx2*ly2.eq.0)goto 99

        if(lx2+ly2.eq.2)then
          temp=xp(1)
          xp(1)=xp(2)
          xp(2)=temp
          temp=yp(1)
          yp(1)=yp(2)
          yp(2)=temp
          ltemp=lx1
          lx1=lx2
          lx2=ltemp
          ltemp=ly1
          ly1=ly2
          ly2=ltemp
        end if

        if(xp(2)-xp(1).ne.0.d0)then
          eps=(xmax-xp(1))/(xp(2)-xp(1))
          if(eps.gt.1.0.or.eps.lt.0.0)goto 10
          ynew=yp(1)+eps*(yp(2)-yp(1))
          if(ynew.gt.ymax.or.ynew.lt.ymin)goto 10
          yp(2)=ynew
          xp(2)=xmax
          call j06baf(xp,yp,2,istyl,isymb,ifail)
          goto 999

   10     continue

          eps=(xmin-xp(1))/(xp(2)-xp(1))
          if(eps.gt.1.0.or.eps.lt.0.0)goto 20
          ynew=yp(1)+eps*(yp(2)-yp(1))
          if(ynew.gt.ymax.or.ynew.lt.ymin)goto 20
          yp(2)=ynew
          xp(2)=xmin
          call j06baf(xp,yp,2,istyl,isymb,ifail) 
          goto 999
        end if

   20   continue

        if(yp(2)-yp(1).ne.0.d0)then
          eps=(ymax-yp(1))/(yp(2)-yp(1))
          if(eps.gt.1.0.or.eps.lt.0.0)goto 30
          xnew=xp(1)+eps*(xp(2)-xp(1))
          if(xnew.gt.xmax.or.xnew.lt.xmin)goto 30
          yp(2)=ymax
          xp(2)=xnew
          call j06baf(xp,yp,2,istyl,isymb,ifail) 
          goto 999
 
   30     continue 

          eps=(ymin-yp(1))/(yp(2)-yp(1))
          if(eps.gt.1.0.or.eps.lt.0.0)goto 40
          xnew=xp(1)+eps*(xp(2)-xp(1))
          if(xnew.gt.xmax.or.xnew.lt.xmin)goto 40
          yp(2)=ymin
          xp(2)=xnew
          call j06baf(xp,yp,2,istyl,isymb,ifail) 
          goto 999
        end if

   40   continue
        print*,'ERROR in logic'
        stop
         
c
c-------No points in plotting area
c
   99   continue

        ipts=0

        if(xp(2)-xp(1).ne.0.d0)then
          eps=(xmax-xp(1))/(xp(2)-xp(1))
          if(eps.gt.1.0.or.eps.lt.0.0)goto 60
          ynew=yp(1)+eps*(yp(2)-yp(1))
          if(ynew.le.ymax.and.ynew.ge.ymin)then
            ipts=1
            ypp(1)=ynew
            xpp(1)=xmax
          end if

   60     continue

          eps=(xmin-xp(1))/(xp(2)-xp(1))
          if(eps.gt.1.0.or.eps.lt.0.0)goto 70
          ynew=yp(1)+eps*(yp(2)-yp(1))
          if(ynew.le.ymax.and.ynew.ge.ymin)then
            ipts=ipts+1
            ypp(ipts)=ynew
            xpp(ipts)=xmin
          end if
          if(ipts.eq.2)then
            call j06baf(xpp,ypp,2,istyl,isymb,ifail) 
            goto 999
          end if
        end if

   70   continue

        if(yp(2)-yp(1).ne.0.d0)then
          eps=(ymax-yp(1))/(yp(2)-yp(1))
          if(eps.gt.1.0.or.eps.lt.0.0)goto 80
          xnew=xp(1)+eps*(xp(2)-xp(1))
          if(xnew.lt.xmax.and.xnew.gt.xmin)then
            ipts=ipts+1
            ypp(ipts)=ymax
            xpp(ipts)=xnew
          end if
          if(ipts.eq.2)then
            call j06baf(xpp,ypp,2,istyl,isymb,ifail) 
            goto 999
          end if
 
   80     continue 

          eps=(ymin-yp(1))/(yp(2)-yp(1))
          if(eps.gt.1.0.or.eps.lt.0.0)goto 90
          xnew=xp(1)+eps*(xp(2)-xp(1))
          if(xnew.lt.xmax.and.xnew.gt.xmin)then
            ipts=ipts+1 
            ypp(ipts)=ymin
            xpp(ipts)=xnew 
          end if
          if(ipts.eq.2)then
            call j06baf(xpp,ypp,2,istyl,isymb,ifail) 
            goto 999
          end if
        end if

   90   continue

        if(ipts.eq.0)then
          goto 999
        end if

        if(ipts.eq.1)then
          goto 999
        end if

  999   continue
  777   format(1x,a7,8(1x,f9.6))

    1 continue
      return
      end

