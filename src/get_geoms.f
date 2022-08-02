      subroutine  get_funcurv_geom_uni(a,b,nch,k,npts,adjs,srcinfo,
     1   srccoefs,ts1,qwts,norders,iptype,ixys,funcurv,ndd,dpars,ndz,
     2   zpars,ndi,ipars)
      implicit real *8 (a-h,o-z)
      integer nch,k,kg,nover,npts,nptsg,npts_over
      integer adjs(2,nch)
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

        adjs(1,ich) = ich-1
        adjs(2,ich) = ich+1
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
          ds = sqrt(srcinfo(3,ipt)**2 + srcinfo(4,ipt)**2)
          qwts(ipt) = wts(j)*ds
        enddo
        istart = (ich-1)*k+1

        call dgemm('n','t',6,k,k,alpha,srcinfo(1,istart),
     1    8,umat,k,beta,srccoefs(1,istart),6)
        norders(ich) = k
        iptype(ich) = 1
        ixys(ich) = (ich-1)*k+1
      enddo
      ixys(nch+1) = npts+1
      adjs(1,1) = nch
      adjs(2,nch) = 1

      
      return
      end
c
c
c
c
c
c
      subroutine get_diamond(nch0,nch,k,npts,adjs,srcinfo,srccoefs,qwts,
     1  norders,iptype,ixys)
      implicit real *8 (a-h,o-z)
      integer adjs(2,nch)
      real *8 srcinfo(8,npts),srccoefs(6,npts),qwts(npts)
      integer norders(nch),iptype(nch),ixys(nch+1)
      integer, allocatable :: ixys0(:),adjs0(:,:)
      integer ipars(1)
      real *8 dpars(3)
      real *8, allocatable :: ts1(:)
      external absconv_geom

      allocate(ixys0(nch0+1),adjs0(2,nch0))



      do i=1,nch+1
        ixys(i) = (i-1)*k + 1
      enddo
      ndd = 3
      ndi = 1
      dpars(1) = 0.05d0
      dpars(2) = -1.0d0
      dpars(3) = 0.0d0
      a = -1.0d0
      b = 1.0d0
      npts0 = nch0*k
      allocate(ts1(npts0))

      do iedge=1,4
        ii = (iedge-1)*nch0 + 1
        istart = ixys(ii) 
        call get_funcurv_geom_uni(a,b,nch0,k,npts0,adjs0,
     1   srcinfo(1,istart),
     1   srccoefs(1,istart),ts1,qwts(istart),
     2   norders(ii),iptype(ii),ixys0,absconv_geom,ndd,dpars,ndz,
     2   zpars,ndi,iedge)
      enddo

      do i=1,nch
        adjs(1,i) = i-1
        adjs(2,i) = i+1
      enddo
      adjs(1,1) = nch
      adjs(2,nch) = 1

      
      return
      end
c
c
c
c 
c
c
c
c
      subroutine get_diamond_many(nch0,ncomp,rsc,shifts,
     1  nch,k,npts,adjs,srcinfo,srccoefs,qwts,
     1  norders,iptype,ixys)
      implicit real *8 (a-h,o-z)
      integer ncomp
      real *8 rsc(ncomp),shifts(2,ncomp)
      integer adjs(2,nch)
      real *8 srcinfo(8,npts),srccoefs(6,npts),qwts(npts)
      integer norders(nch),iptype(nch),ixys(nch+1)
      integer, allocatable :: ixys0(:),adjs0(:,:)
      integer ipars(1)
      real *8 dpars(3)
      real *8, allocatable :: ts1(:)

      nchcomp = 4*nch0
      npts0 = nchcomp*k
      do i=1,nch+1
        ixys(i) = (i-1)*k+1
      enddo


      allocate(ixys0(nchcomp+1),adjs0(2,nchcomp))
      do icomp=1,ncomp
        ichstart = (icomp-1)*nchcomp+1
        ichend = icomp*nchcomp
        istart = ixys(ichstart)
        iend = ixys(ichend+1)-1
        call get_diamond(nch0,nchcomp,k,npts0,adjs0,srcinfo(1,istart),
     1   srccoefs(1,istart),qwts(istart),norders(ichstart),
     2   iptype(ichstart),ixys0)
        do ipt=istart,iend
          srcinfo(1,ipt) = srcinfo(1,ipt)*rsc(icomp) + shifts(1,icomp)
          srcinfo(2,ipt) = srcinfo(2,ipt)*rsc(icomp) + shifts(2,icomp)
          srcinfo(3,ipt) = srcinfo(3,ipt)*rsc(icomp)
          srcinfo(4,ipt) = srcinfo(4,ipt)*rsc(icomp)
          srcinfo(5,ipt) = srcinfo(5,ipt)*rsc(icomp)
          srcinfo(6,ipt) = srcinfo(6,ipt)*rsc(icomp)
          srccoefs(1,ipt) = srccoefs(1,ipt)*rsc(icomp) + shifts(1,icomp)
          srccoefs(2,ipt) = srccoefs(2,ipt)*rsc(icomp) + shifts(2,icomp)
          srccoefs(3,ipt) = srccoefs(3,ipt)*rsc(icomp)
          srccoefs(4,ipt) = srccoefs(4,ipt)*rsc(icomp)
          srccoefs(5,ipt) = srccoefs(5,ipt)*rsc(icomp)
          srccoefs(6,ipt) = srccoefs(6,ipt)*rsc(icomp)
          qwts(ipt) = qwts(ipt)*rsc(icomp)
        enddo
        do ich=ichstart,ichend
          adjs(1,ich) = adjs0(1,ich-ichstart+1) + ichstart-1
          adjs(2,ich) = adjs0(2,ich-ichstart+1) + ichstart-1
        enddo
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
      
      a = dpars(1)
      b = dpars(2)
      ct = cos(t)
      st = sin(t)
      srcinfo(1) = a*ct
      srcinfo(2) = b*st
      srcinfo(3) = -a*st
      srcinfo(4) = b*ct
      srcinfo(5) = -a*ct
      srcinfo(6) = -b*st
      ds = sqrt(srcinfo(3)**2 + srcinfo(4)**2)
      srcinfo(7) = srcinfo(4)/ds
      srcinfo(8) = -srcinfo(3)/ds

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
c      subroutine load_cavity_zpars(a,b,m,zpars,ndz)
c      implicit real *8 (a-h,o-z)
c      complex *16 zpars(2*m)
c      real *8, allocatable :: s(:),r(:),th(:),rho(:),z(:)
c      complex *16, allocatable :: zwork(:)
c      complex *16 ima
c      data ima/(0.0d0,1.0d0)/
c
c      done = 1
c      pi = atan(done)*4
c
c      ndz = 2*m
c      lw = 10*ndz + 15
c      allocate(s(m),r(m),th(m),rho(m),z(m))
c      allocate(zwork(lw))
c      call zffti(ndz,zwork)
c      c = a
c      hh = c/sqrt(2.0d0)
c      aa = 1.0d0
c      bb = 0.0d0
c      do i=1,m
c        s(i) = (i-0.5d0)/(m+0.0d0)*pi
c        suse = (s(i)-pi/2)/a
c        call qerrfun(suse,rr)
c        r(i) = 1-a*rr
c        suse = (s(i)-pi/2)
c        call crn_fconvgauss(suse, aa, bb, hh, rr, tmp, tmp2)
c        th(i) = (b-a) + 2*(1.0d0-(b-a)/pi)*rr
c        rho(i) = r(i)*sin(th(i))
c        z(i) = r(i)*cos(th(i))
c      enddo
c      do i=1,m
c        zpars(i) = rho(i) + ima*z(i) 
c        zpars(i+m) = -rho(m-i+1) + ima*(m-i+1)
c      enddo
c      call zfftf(ndz,zpars,wsave)
c      do i=1,ndz
c        zpars(i) = zpars(i)/(m+0.0d0)
c      enddo
c
c
c      return
c      end
