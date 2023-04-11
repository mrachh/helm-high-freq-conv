
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the code below contains several user-callable routines
c       for the construction of polygons and polygons with
c       smoothed corners. they are as follows:
c
c       corners_pack - 
c
c       chunkpolysmooth - given a set of vertices on a closed or
c           open curve, constructs a smoothed version of the geometry
c           using convolutional and adaptive discretization methods
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
c
c
        subroutine chunkpolysmooth(ier, eps, widths, ibell, p1, p2, 
     1      i1, i2, nverts, verts, ifclosed, nover, k, nch, chunks,
     2      adjs, ders, ders2, hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        real *8 widths(1), verts(2,1), chunks(2,k,1), ders(2,k,1),
     1      ders2(2,k,1), hs(1), xnodes(1000), whts(1000),
     2      pars(10000), xymid(100), cmass(100)
        real *8, allocatable :: evec(:,:)
        external fgauss
c
        allocate( evec(2,nverts+10) )
c
c       construct a polygon with rounded corners
c
c       input:
c
c         eps - the precision to which the curve should be resolved
c         widths - the size to cut out at each corner, should not be larger
c           than half of any adjacent line segment, should contain nverts
c           values, with widths(1)=widths(nverts) if a closed curve.
c           one should set widths(1)=widths(nverts)=0 if the curve is open.
c         ibell - the type of rounding to use
c             ibell=1    Gaussian
c             ibell=2    not implemented 
c             ibell=3    not implemented
c         p1,p2,i1,i2 - parameters for use in the bell
c             ibell=1    no parameters needed
c             ibell=2    i1=korder, the order of the bell
c         nverts - the number of vertices, if curve is closed then
c             we expect that verts(1) = verts(nverts)
c         verts - (2,nverts) array containing vertices in the plane, assumed
c             assumed to be oriented counterclockwise
c         ifclosed - set to 1 if the curve is closed, 0 if open
c         nover - number of times to oversample the chunks upon completion
c         k - number of gaussian nodes to put on each panel
c
c       output:
c
c         nch - number of chunks created
c         chunks - array of node locations on the chunks
c         adjs - adjancency information, adjs(1,i) is the chunk to the
c           chunk to the left of chunk i, adjs(2,i) to the right
c         ders - derivatives, scaled w.r.t. hs
c         ders2 - 2nd derivatives, scaled w.r.t. hs
c         hs - scaling parameters such that for t in [-1,1]
c
c       NOTE: in order to calculate dsdt relative to t in [-1,1], the
c         the scaling factors must be used:
c
c                   dsdt = sqrt(ders(1)**2+ders(2)**2)*hs
c
c
        done=1
        pi=4*atan(done)

c
c       check that widths is of proper size
c
        ier=0
        do 1400 i=1,nverts-1
c
        dx1=verts(1,i+1)-verts(1,i)
        dy1=verts(2,i+1)-verts(2,i)
        dx2=verts(1,i)-verts(1,i+1)
        dy2=verts(2,i)-verts(2,i+1)
c
        r1=sqrt(dx1**2+dy1**2)
        r2=sqrt(dx2**2+dy2**2)
c
        if (widths(i+1) .gt. r1/2) ier=2
        if (widths(i+1) .gt. r2/2) ier=2
c
 1400 continue
c
        if (ier .ne. 0) then
            call prinf('widths is too large, ier=*',ier,1)
            return
        endif

c
c       construct unit vectors which point along the edges,
c       will be useful
c
        do i = 1,nverts-1
          dx1=verts(1,i+1)-verts(1,i)
          dy1=verts(2,i+1)-verts(2,i)
          r1=sqrt(dx1**2+dy1**2)
          evec(1,i)=dx1/r1
          evec(2,i)=dy1/r1
        end do


        call prin2('in chkpsmooth, evec = *', evec, 2*(nverts-1))

c
c       construct all chunks and points, then scan them and
c       divide if necessary - start with the first segment, and
c       proceed - not the first corner
c
        ifwhts=1
        call legerts(ifwhts,k,xnodes,whts)
c
        if (ifclosed .eq. 0) then
          widths(1)=0
          widths(nverts)=0
        end if
c
        nch=0
        nedges=nverts-1
c

        ifstop = 0
        
        do 5000 iseg=1,nedges
c
c       construct the single chunk on the middle part of the
c       edge
c
        x1=verts(1,iseg)+widths(iseg)*evec(1,iseg)
        y1=verts(2,iseg)+widths(iseg)*evec(2,iseg)
c
        x2=verts(1,iseg+1)-widths(iseg+1)*evec(1,iseg)
        y2=verts(2,iseg+1)-widths(iseg+1)*evec(2,iseg)
c
        nch=nch+1
c
        slopex=(x2-x1)
        slopey=(y2-y1)
        r=sqrt((x2-x1)**2+(y2-y1)**2)
        hs(nch)=r/2
c
        do j = 1,k
          chunks(1,j,nch)=x1+(xnodes(j)+1)/2*slopex
          chunks(2,j,nch)=y1+(xnodes(j)+1)/2*slopey
          ders(1,j,nch)=slopex/r
          ders(2,j,nch)=slopey/r
          ders2(1,j,nch)=0
          ders2(2,j,nch)=0
        end do
c
        if (iseg .eq. 1) then
          adjs(1,nch)=-1
          adjs(2,nch)=-1
          goto 3100
        endif

ccc        call prin2('chunks = *', chunks, 2*nch)
ccc        ifstop = ifstop + 1
ccc        if (ifstop .eq. 2) stop
        
c
c       find the last segment created on a rounded corner
c       and update adjacencies
c
        do 2600 i=1,nch-1
          if (adjs(2,i) .lt. 0) then
            iright=i
            goto 2700
          endif
 2600   continue
 2700   continue
c
        adjs(2,iright)=nch
        adjs(1,nch)=iright
        adjs(2,nch)=-1
 3100   continue

        if ((ifclosed .ne. 1) .and. (iseg .eq. nedges)) goto 5100


c
c       . . . now construct a corner chunk
c
        x=verts(1,iseg+1)
        y=verts(2,iseg+1)
        w=widths(iseg+1)
c
        if (iseg .ne. nedges) then
          x3=verts(1,iseg+1)+w*evec(1,iseg+1)
          y3=verts(2,iseg+1)+w*evec(2,iseg+1)
        endif
c
        if (iseg .eq. nedges) then
            x3=verts(1,1)+w*evec(1,1)
            y3=verts(2,1)+w*evec(2,1)
        endif
c
        r23=sqrt((x3-x2)**2+(y3-y2)**2)
        b=sqrt(w*w-(r23/2)**2)
        a=-b/(r23/2)
        xmid=(x2+x3)/2
        ymid=(y2+y3)/2


c
c       calculate the bandwidth parameter h depending on which bell
c       specified by the user
c

        if (ibell .eq. 1) then
c
c       override the scale parameter so eps is 5.0 x 10^{-15}
c
          h = abs(b/a)/8
          goto 1235

cccc            call prin2('h = *', h, 1)
cccc            stop
c
            wd = r23/2
            thresh=eps/50
c
            h=wd/5
            do 1234 i=1,30
            h=sqrt(wd*wd/(-2*log(sqrt(2*pi)*thresh*h)))
            f1=1/sqrt(2*pi)/h*exp(-wd**2/2/h**2)
            if (abs(f1) .le. eps/10) goto 1235
 1234 continue
c
            call prinf('bell width did not converge!!! last h=*',h,1)
            call prin2('eps=*',eps,1)
            call prin2('f1=*',f1,1)
            stop
c
 1235 continue
c
        endif

c
c
c       create the chunk using either a gaussian bell or a finite bell
c
        ta=b/a
        tb=-b/a
        chsmall=1000
        nover5=1
        ifc=0
c
        
        if (ibell .eq. 1) then
c
          pars(1)=a
          pars(2)=b
cccc          pars(3)=h/2
          pars(3)=h

          nch5 = 0
c
          call chunkfunc(eps,ifc,chsmall,ta,tb,fgauss,
     2        pars,nover5,k,nch5,chunks(1,1,nch+1),adjs(1,nch+1),
     3        ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
c
        end if

        nch999 = nch+nch5
cc        call prinf('nch = *', nch, 1)
cc        call prinf('nch5 = *', nch5, 1)
cc        call prinf('nch999 = *', nch999, 1)
cc        call prin2('after chunkfunc, chunks = *', chunks, 2*nch999*k)
cc        stop

        
c

c
c       if left hand turn, reverse the orientation
c
        u1=x-x2
        u2=y-y2
        v1=x3-x
        v2=y3-y
c
        z=u1*v2-u2*v1
        if (z .ge. 0) ileft=1
        if (z .lt. 0) ileft=0
c
cccc        call prinf('test for left turn, ileft=*',ileft,1)
c
        if (ileft .eq. 1) then
            call chunkreverse(k,nch5,chunks(1,1,nch+1),adjs(1,nch+1),
     1          ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
        endif


cc        call prin2('after reversal, chunks = *', chunks, 2*nch999*k)
cc        stop
        
c
c       try new rotation scheme - calculate the vectors joining
c       the ends of the bump, and the ends of the polygon cut
c
        if (ileft .eq. 1) then
            yminusz1=x2-x3
            yminusz2=y2-y3
        endif
c
        if (ileft .eq. 0) then
            yminusz1=x3-x2
            yminusz2=y3-y2
        endif
c
        w1=abs(2*b/a)
        w2=0
c
        wlen=sqrt(w1**2+w2**2)
        yzlen=sqrt(yminusz1**2+yminusz2**2)
c
        phi=atan2(yminusz2,yminusz1)
cccc        if (ileft .eq. 0) phi=-phi

c
c       now rotate and translate the corner into place
c
        cphi=cos(phi)
        sphi=sin(phi)
c
c        call prin2('x2 = *', x2, 1)
c        call prin2('x3 = *', x3, 1)
c        call prin2('y2 = *', y2, 1)
c        call prin2('y3 = *', y3, 1)
c        call prin2('yminusz1 = *', yminusz1, 1)
c        call prin2('yminusz2 = *', yminusz2, 1)
c        call prin2('w1 = *', w1, 1)
c        call prin2('w2 = *', w2, 1)
c        call prin2('wlen = *', wlen, 1)
c        call prin2('yzlen = *', yzlen, 1)
c        call prin2('phi = *', phi, 1)
c        call prin2('cphi = *', cphi, 1)
c        call prin2('sphi = *', sphi, 1)
c        call prin2('b = *', b, 1)

c        call prin2('x = *', x, 1)
c        call prin2('y = *', y, 1)



c        call prin2('befpre rotation, chunks = *', chunks, 2*k*2)
        
        do i = nch+1,nch+nch5

cc          call prinf('i = *', i, 1)
          do j = 1,k
c
cc            call prin2('before, chxy = *', chunks(1,j,i), 2)

cccc            chunks(2,j,i) = chunks(2,j,i) - b
c
            x7 = chunks(1,j,i)
            y7 = chunks(2,j,i) - b
cc            call prin2('x7 = *', x7, 1)
cc            call prin2('y7 = *', y7, 1)

cccc            call prin2('x7 = *', x7, 1)
            chunks(1,j,i)=cphi*x7-sphi*y7
            chunks(2,j,i)=sphi*x7+cphi*y7

cccc            call prin2('test pt = *', chunks(1,j,i), 2)
cccc            stop
c
            chunks(1,j,i) = chunks(1,j,i) + x
            chunks(2,j,i) = chunks(2,j,i) + y
c
cc            call prin2('after, chxy = *', chunks(1,j,i), 2)
cc            print *

            dx7=ders(1,j,i)
            dy7=ders(2,j,i)
            ders(1,j,i)=(cphi*dx7-sphi*dy7)
            ders(2,j,i)=(sphi*dx7+cphi*dy7)
c
            d2x7=ders2(1,j,i)
            d2y7=ders2(2,j,i)
            ders2(1,j,i)=cphi*d2x7-sphi*d2y7
            ders2(2,j,i)=sphi*d2x7+cphi*d2y7
c
          enddo
cccc          call prin2('after rotation, chunks = *', chunks, 2*k*2)
cccc          stop
        enddo


cccc        call prin2('after rotation, chunks = *', chunks, 2*k*nch999)
cccc        stop
        
c
        do i = nch+1,nch+nch5
          if (adjs(1,i) .gt. 0) adjs(1,i)=adjs(1,i)+nch
          if (adjs(2,i) .gt. 0) adjs(2,i)=adjs(2,i)+nch
        enddo
c
        do i = nch+1,nch+nch5
          if (adjs(1,i) .lt. 0) ileft=i
          if (adjs(2,i) .lt. 0) iright=i
        enddo
c
        adjs(1,ileft)=nch
        adjs(2,nch)=ileft
        nch=nch+nch5

cccc        call prin2('after rotation, chunks = *', chunks, 60)
cccc        stop
c

c
 5000 continue
 5100 continue



c
c       if the curve is closed update the adjacency info
c
        if (ifclosed .ne. 1) goto 5700
c
        do i=1,nch
          if (adjs(1,i) .lt. 0) ileft=i
          if (adjs(2,i) .lt. 0) iright=i
        enddo
c
        adjs(1,ileft)=iright
        adjs(2,iright)=ileft
c
 5700 continue

c
c       now just check that no neighboring chunks are off by
c       more than a factor of two in arclength
c
        maxiter = 100000
        do ijk = 1,maxiter
c
          nchold = nch
          ifdone = 1
          dlen = -1
          id_chunk = -1
c
            do i = 1,nchold
c
c             find the longest chunk that requires splitting
c
              i1=adjs(1,i)
              i2=adjs(2,i)

c
c             calculate chunk lengths
c
              call chunksize(k,chunks(1,1,i),ders(1,1,i),hs(i),
     1          xymid,rad,cmass,crad,rlself)
c
              if (i1 .gt. 0) then
                call chunksize(k,chunks(1,1,i1),ders(1,1,i1),hs(i1),
     1            xymid,rad,cmass,crad,rl1)
              endif
c
              if (i2 .gt. 0) then
                call chunksize(k,chunks(1,1,i2),ders(1,1,i2),hs(i2),
     1            xymid,rad,cmass,crad,rl2)
              endif
c
c             only check if self is larger than either of adjacent blocks
c             iterating a couple times will catch everything
c
              ifsplit=0
              sc=2.05d0
cccc              sc=2
c
              if (i1 .gt. 0) then
                if (rlself .gt. sc*rl1) ifsplit=1
              endif
c
              if (i2 .gt. 0) then
                if (rlself .gt. sc*rl2) ifsplit=1
              endif
c
              if (ifsplit .ne. 0) then
                ifdone = 0
                if (rlself .gt. dlen) then
                  id_chunk = i
                  dlen = rlself
                endif
              endif
c
          enddo
c
          if (ifdone .eq. 1) goto 9100
c
c         if not done, split id_chunk
c
          call chunksplit1(id_chunk, k, nch, chunks, adjs, ders,
     1      ders2, hs)
c
        enddo
c
        call prinf('bomb! maxiter too small in smoothpoly!*', pi, 0)
        stop
c
 9100 continue

c        iw=77
c        itype=1
c        nnn = k*nch
c        call zpyplot(iw, chunks, nnn, itype, 'before oversampling*')
c        call prinf('before oversampling, nch = *', nch, 1)
c
cccc        stop

!!!!        call prin2('before oversampling, chunks = *', chunks, 2*nnn)

c
c       and finally oversample the curve by a factor of nover
c
        if (nover .le. 1) return
        niter=nover-1
c
        do i = 1,niter       
          nchold = nch
          do jj = 1,nchold
            call chunksplit1(jj, k, nch, chunks, adjs, ders,
     1        ders2, hs)
          enddo
        enddo
c
c        iw=78
c        itype=1
c        nnn = k*nch
c        call zpyplot(iw, chunks, nnn, itype, 'after oversampling*')
c        call prinf('after oversampling, nch = *', nch, 1)
c
cccc        call prin2('after oversampling, chunks = *', chunks, 2*nnn)
c        stop

        deallocate( evec )
c
        return
        end
c
c
c
c
c
        subroutine fgauss(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)
c
c       wrapper for the routine crn_fconvgauss
c
        a=pars(1)
        b=pars(2)
        h=pars(3)
c
        call crn_fconvgauss(t,a,b,h,val,der,der2)
        x=t
        y=val
c
        dxdt=1
        dydt=der
c
        dxdt2=0
        dydt2=der2
c
        return
        end
c
c
c
c
c
        subroutine crn_fconvgauss(x, a, b, h, val, der, der2)
        implicit real *8 (a-h,o-z)
c
c       this routine computes the convolution
c
c         ( a*abs(x)+b ) \star 1/(sqrt(2*pi)*h^2) exp(-x^2/(2*h^2))
c
c       this effectively smoothes off the corner from the abs(x) function
c
c
        done=1
        two=2
        pi=4*atan(done)
c
c       . . . formulas are computed via maple
c
        x2=x/sqrt(two)/h
        call qerrfun(x2,verf)
        val=a*x*verf+b+sqrt(two/pi)*a*h*exp(-x*x/two/h/h)
c
        fnorm=1/sqrt(2*pi)*exp(-x*x/2)
        der=a*verf
        der2=a*sqrt(two/pi)/h*exp(-x*x/two/h/h)
c
        return
        end
c
c
c
c
c
