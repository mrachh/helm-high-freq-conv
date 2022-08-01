      subroutine  get_funcurv_geom_uni(a,b,nch,k,npts,srcinfo,
     1   srccoefs,ts1,qwts,norders,iptype,ixys,funcurv,ndd,dpars,ndz,
     2   zpars,ndi,ipars)
      implicit real *8 (a-h,o-z)
      integer nch,k,kg,nover,npts,nptsg,npts_over
      real *8 srcinfo(8,npts)
      real *8 srccoefs(6,npts)
      real *8 ts1(npts)
      real *8 qwts(npts)
      integer norders(nch),iptype(nch)
      integer ixys(nch+1)
      complex *16 zpars(ndz)
      real *8 dpars(ndd)
      integer ipars(ndi)
      external funcurv

      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)

      done = 1.0d0
      pi = atan(done)*4.0d0
      h = (b-a)/(nch+0.0d0)
      alpha = 1.0d0
      beta = 0.0d0

      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      itype = 2
      call legeexps(itype,k,ts,umat,vmat,wts)

      do ich=1,nch
        tstart = a+(ich-1.0d0)*h
        tend = a+(ich+0.0d0)*h
c
c  get source info
c

        do j=1,k
          ipt = (ich-1)*k + j
          tuse = (tstart+tend)/2 + (tend-tstart)*ts(j)/2
          ts1(ipt) = tuse
          call funcurv(tuse,ndd,dpars,ndz,zpars,ndi,ipars,
     1       srcinfo(1:8,ipt))
          
          srcinfo(3,ipt) = srcinfo(3,ipt)*h/2
          srcinfo(4,ipt) = srcinfo(4,ipt)*h/2
          srcinfo(5,ipt) = srcinfo(5,ipt)*h*h/2/2
          srcinfo(6,ipt) = srcinfo(6,ipt)*h*h/2/2
          qwts(ipt) = h/2*wts(j)
        enddo
        istart = (ich-1)*k+1

        call dgemm('n','t',6,k,k,alpha,srcinfo(1,istart),
     1    8,umat,k,beta,srccoefs(1,istart),6)
        norders(ich) = k
        iptype(ich) = 1
        ixys(ich) = (ich-1)*k+1
      enddo
      ixys(nch+1) = npts+1

      
      return
      end

      subroutine get_diamond(nch0,nch,k,npts,srcinfo,srccoefs,qwts,
     1  norders,iptype,ixys)
      implicit real *8 (a-h,o-z)
      real *8 srcinfo(8,npts),srccoefs(6,npts),qwts(npts)
      integer norders(nch),iptype(nch),ixys(nch+1)
      integer, allocatable :: ixys0(:)
      integer ipars(1)
      real *8 dpars(3)
      real *8, allocatable :: ts1(:)
      external absconv_geom

      allocate(ixys0(nch0+1))



      do i=1,nch+1
        ixys(i) = (i-1)*k + 1
      enddo
      ndd = 3
      ndi = 1
      dpars(1) = 0.1d0
      dpars(2) = -1.0d0
      dpars(3) = 0.0d0
      a = -1.0d0
      b = 1.0d0
      npts0 = nch0*k
      allocate(ts1(npts0))

      do iedge=1,4
        ii = (iedge-1)*nch0 + 1
        istart = ixys(ii) 
        call get_funcurv_geom_uni(a,b,nch0,k,npts0,srcinfo(1,istart),
     1   srccoefs(1,istart),ts1,qwts(istart),
     2   norders(ii),iptype(ii),ixys0,absconv_geom,ndd,dpars,ndz,
     2   zpars,ndi,iedge)
      enddo

      
      return
      end
c
c
c
c 



      subroutine circ_geom(t,ndd,dpars,ndz,zpars,ndi,ipars,srcinfo)
      implicit real *8 (a-h,o-z)
      real *8 dpars(ndd),srcinfo(8)
      integer ipars(ndi)
      complex *16 zpars(ndz)
      
      r = dpars(1)
      ct = cos(t)
      st = sin(t)
      srcinfo(1) = r*ct
      srcinfo(2) = r*st
      srcinfo(3) = -r*st
      srcinfo(4) = r*ct
      srcinfo(5) = -r*ct
      srcinfo(6) = -r*st
      srcinfo(7) = ct
      srcinfo(8) = st

      return
      end



      subroutine absconv_geom(t,ndd,dpars,ndz,zpars,ndi,ipars,srcinfo)
      implicit real *8 (a-h,o-z)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 srcinfo(8)
      
      h = dpars(1)
      a0 = dpars(2)
      b0 = dpars(3)
      iedge = ipars(1)


      call crn_fconvgauss(t,a0,b0,h,val,der,der2) 
      if(iedge.eq.1) then
        srcinfo(1) = -t
        srcinfo(2) = val +2.0d0
        srcinfo(3) = -1.0d0
        srcinfo(4) = der
        srcinfo(5) = 0.0d0
        srcinfo(6) = der2
      else if(iedge.eq.2) then
        srcinfo(1) = -val -2.0d0
        srcinfo(2) = -t
        srcinfo(3) = -der
        srcinfo(4) = -1.0d0
        srcinfo(5) = -der2
        srcinfo(6) = 0.0d0
      else if(iedge.eq.3) then
        srcinfo(1) = t
        srcinfo(2) = -val-2.0d0
        srcinfo(3) = 1.0d0
        srcinfo(4) = -der
        srcinfo(5) = 0.0d0
        srcinfo(6) = -der2
      else if(iedge.eq.4) then
        srcinfo(1) = val+2.0d0
        srcinfo(2) = t
        srcinfo(3) = der
        srcinfo(4) = 1.0d0
        srcinfo(5) = der2
        srcinfo(6) = 0.0d0
      endif

      ds = sqrt(srcinfo(3)**2 + srcinfo(4)**2)
      srcinfo(7) = srcinfo(4)/ds
      srcinfo(8) = -srcinfo(3)/ds


      return
      end

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
