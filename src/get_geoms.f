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
          srccoefs(1,ipt) = srccoefs(1,ipt)*rsc(icomp)
          srccoefs(2,ipt) = srccoefs(2,ipt)*rsc(icomp)
          srccoefs(3,ipt) = srccoefs(3,ipt)*rsc(icomp)
          srccoefs(4,ipt) = srccoefs(4,ipt)*rsc(icomp)
          srccoefs(5,ipt) = srccoefs(5,ipt)*rsc(icomp)
          srccoefs(6,ipt) = srccoefs(6,ipt)*rsc(icomp)
          qwts(ipt) = qwts(ipt)*rsc(icomp)
        enddo
        do ich=ichstart,ichend
          adjs(1,ich) = adjs0(1,ich-ichstart+1) + ichstart-1
          adjs(2,ich) = adjs0(2,ich-ichstart+1) + ichstart-1
          istart = ixys(ich)
          srccoefs(1:2,istart) = srccoefs(1:2,istart) +
     1       shifts(1:2,icomp)
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
        done=1.0d0
        two=2.0d0
        pi=4*atan(done)
c
c       . . . formulas are computed via maple
c
        x2=x/sqrt(two)/h
        call qerrfun(x2,verf)
        val=a*x*verf+b+sqrt(two/pi)*a*h*exp(-x*x/two/h/h)
c
        fnorm=1.0d0/sqrt(2.0d0*pi)*exp(-x*x/2)
        der=a*verf
        der2=a*sqrt(two/pi)/h*exp(-x*x/two/h/h)
c
        return
        end

c
c
c
c
      subroutine interp_dens(nch,norders,ixys,iptype,npts,
     1  solncoefs,ninterp,ts_interp,ich_interp,soln_interp)
c
c  given coefs of a density defined on a grid, and a set
c  of points identified through their local t coordinates
c  and chunk id, compute interpolated density
c
c

      implicit real *8 (a-h,o-z)
      integer nch,norders(nch),ixys(nch+1),iptype(nch),npts
      complex *16 solncoefs(npts),soln_interp(ninterp)
      integer ich_interp(ninterp)
      real *8 ts_interp(ninterp)

      real *8, allocatable :: pols(:)

      kmax = maxval(norders(1:nch))
      allocate(pols(kmax))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(pols,i,ich,istart,j)
      do i=1,ninterp
        ich = ich_interp(i)
        istart = ixys(ich)
        do j=1,norders(ich)
          pols(j) = 0
        enddo
        if(iptype(ich).eq.1) then 
          call legepols(ts_interp(i),norders(ich)-1,pols)
        endif
        soln_interp(i) = 0
        do j=1,norders(ich)
          soln_interp(i) = soln_interp(i) +
     1       solncoefs(istart+j-1)*pols(j)
        enddo
      enddo
C$OMP END PARALLEL DO

      
      return
      end



      subroutine get_circ_dens_error(dpars,nch1,k1,npts1,solncoefs1,
     1   nch2,k2,npts2,solncoefs2,k3,nch3,erra,errq)
c
c   Given two densities defined on two grids, and a third
c   reference grid to compute the error, evalute
c   the L2 error in the density
c
c   and errq is the optimality factor, i.e. project solncoefs1, 
c   onto the basis of solncoefs2 and then evaluate the error
c   between solncoefs1 and the projected version of solncoefs1
c
      implicit real *8 (a-h,o-z)
      integer nch1,k1,npts1,nch2,k2,npts2,k3,nch3
      complex *16 solncoefs1(npts1),solncoefs2(npts2)
      real *8 erra
      integer, allocatable :: adjs1(:,:),adjs2(:,:),adjs3(:,:),
     1  adjs22(:,:)
      real *8, allocatable :: srcinfo1(:,:),srcinfo2(:,:),srcinfo3(:,:),
     1   srcinfo22(:,:)
      real *8, allocatable :: srccoefs1(:,:),srccoefs2(:,:),
     1  srccoefs3(:,:),srccoefs22(:,:)
      real *8, allocatable :: ts1(:),ts2(:),ts3(:),ts22(:)
      real *8, allocatable :: qwts1(:),qwts2(:),qwts3(:),qwts22(:)
      integer, allocatable :: norders1(:),norders2(:),norders3(:),
     1   norders22(:)
      integer, allocatable :: iptype1(:),iptype2(:),iptype3(:),
     1   iptype22(:)
      integer, allocatable :: ixys1(:),ixys2(:),ixys3(:),ixys22(:)
      real *8, allocatable :: srcrad1(:),srcrad22(:)

      real *8, allocatable :: ts_interp1(:),ts_interp2(:)
      integer, allocatable :: ich_interp1(:),ich_interp2(:)
      real *8, allocatable :: dist1(:),dist2(:)

      real *8, allocatable :: ts_interp12(:)
      integer, allocatable :: ich_interp12(:)
      real *8, allocatable :: dist12(:)

      complex *16, allocatable :: soln1(:),soln2(:),soln12(:)
      complex *16, allocatable :: soln12coefs(:),soln123(:)
      real *8, allocatable :: t(:),w(:),umat(:,:),vmat(:,:)

      real *8 dpars(2),timeinfo(3)
      complex *16 zpars
      integer ipars

      external circ_geom

      done = 1.0d0
      pi = atan(done)*4

      npts3 = k3*nch3
      allocate(ts_interp1(npts3),ich_interp1(npts3))
      allocate(ts_interp2(npts3),ich_interp2(npts3))
      allocate(dist1(npts3),dist2(npts3))
      allocate(soln1(npts3),soln2(npts3))


      allocate(adjs1(2,nch1),srcinfo1(8,npts1),srccoefs1(6,npts1))
      allocate(ts1(npts1),qwts1(npts1),norders1(nch1),iptype1(nch1))
      allocate(ixys1(nch1+1),srcrad1(npts1))



      allocate(adjs2(2,nch2),srcinfo2(8,npts2),srccoefs2(6,npts2))
      allocate(ts2(npts2),qwts2(npts2),
     1    norders2(nch2),iptype2(nch2))
      allocate(ixys2(nch2+1))

      k22 = 30
      npts22 = nch2*k22

      allocate(adjs22(2,nch2),srcinfo22(8,npts22),srccoefs22(6,npts22))
      allocate(ts22(npts22),qwts22(npts22),
     1    norders22(nch2),iptype22(nch2))
      allocate(ixys22(nch2+1),srcrad22(npts22))

      allocate(adjs3(2,nch3),srcinfo3(8,npts3),srccoefs3(6,npts3))
      allocate(ts3(npts3),qwts3(npts3),norders3(nch3),iptype3(nch3))
      allocate(ixys3(nch3+1))
      
      ndd_curv = 2
      ndz_curv = 0
      ndi_curv = 0

      a = 0.0d0
      b = 2*pi
      call prin2('dpars=*',dpars,2)
      call get_funcurv_geom_uni(a,b,nch1,k1,npts1,adjs1,
     1  srcinfo1,srccoefs1,ts1,qwts1,norders1,iptype1,ixys1,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      
      call get_funcurv_geom_uni(a,b,nch2,k2,npts2,adjs2,
     1  srcinfo2,srccoefs2,ts2,qwts2,norders2,iptype2,ixys2,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)

      
      call get_funcurv_geom_uni(a,b,nch2,k22,npts22,adjs22,
     1  srcinfo22,srccoefs22,ts22,qwts22,norders22,iptype22,ixys22,
     1  circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)

      call get_funcurv_geom_uni(a,b,nch3,k3,npts3,adjs3,
     1  srcinfo3,srccoefs3,ts3,qwts3,norders3,iptype3,ixys3,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      
      do i=1,npts1
        srcrad1(i) = 0
      enddo

      do i=1,npts22
        srcrad22(i) = 0
      enddo


cc      call prin2('srccoefs1=*',srccoefs1,6*npts1)
      
      
      call findnearchunktarg_id_ts(nch1,norders1,ixys1,iptype1,
     1  npts1,srccoefs1,srcinfo1,srcrad1,8,npts3,srcinfo3,ich_interp1,
     2  ts_interp1,dist1,timeinfo,ier)
      
      call findnearchunktarg_id_ts(nch2,norders22,ixys22,iptype22,
     1  npts22,srccoefs22,srcinfo22,srcrad22,8,npts3,srcinfo3,
     2  ich_interp2,ts_interp2,dist2,timeinfo,ier)

c
c  
c
      allocate(ich_interp12(npts2),ts_interp12(npts2),dist12(npts2))
c      call findnearchunktarg_id_ts_brute(nch1,norders1,ixys1,iptype1,
c     1  npts1,srccoefs1,srcinfo1,8,npts2,srcinfo2,ich_interp12,
c     2  ts_interp12,dist12)
      
      call findnearchunktarg_id_ts(nch1,norders1,ixys1,iptype1,
     1  npts1,srccoefs1,srcinfo1,srcrad1,8,npts2,srcinfo2,ich_interp12,
     2  ts_interp12,dist12,timeinfo,ier)
      allocate(soln12(npts2),soln12coefs(npts2))
      
      call interp_dens(nch1,norders1,ixys1,iptype1,npts1,
     1  solncoefs1,npts2,ts_interp12,ich_interp12,soln12)
c
c
      allocate(t(k2),w(k2),umat(k2,k2),vmat(k2,k2))
      itype = 2
      call legeexps(itype,k2,t,umat,vmat,w)
c
c  convert soln12 to soln12 coefs
c
      alpha = 1.0d0
      beta = 0.0d0
      do i=1,nch2
        istart = ixys2(i)
        call dgemm('n','t',2,k2,k2,alpha,soln12(istart),2,
     1    umat,k2,beta,soln12coefs(istart),2)
      enddo

      call interp_dens(nch1,norders1,ixys1,iptype1,npts1,
     1  solncoefs1,npts3,ts_interp1,ich_interp1,soln1)
c
      call interp_dens(nch2,norders2,ixys2,iptype2,npts2,
     1  solncoefs2,npts3,ts_interp2,ich_interp2,soln2)
      
c
c
      allocate(soln123(npts3))
      call interp_dens(nch2,norders2,ixys2,iptype2,npts2,
     1  soln12coefs,npts3,ts_interp2,ich_interp2,soln123)
      
      erra = 0
      errq = 0
      ra = 0
      do i=1,npts3
        erra = erra + abs(soln1(i)-soln2(i))**2*qwts3(i)
        errq = errq + abs(soln1(i)-soln123(i))**2*qwts3(i)
        ra = ra + abs(soln1(i))**2*qwts3(i)
      enddo
      erra = sqrt(erra/ra)
      errq = sqrt(errq/ra)
      errq = errq/erra



      return
      end





      subroutine get_diamond_many_dens_error(ncomp,shifts,rsc,
     1   nch10,nch1,k1,npts1,solncoefs1,nch20,nch2,k2,npts2,
     1   solncoefs2,k3,nch30,nch3,erra,errq,ifwrite,iunit)
c
c   Given two densities defined on two grids, and a third
c   reference grid to compute the error, evalute
c   the L2 error in the density
c
c   and errq is the optimality factor, i.e. project solncoefs1, 
c   onto the basis of solncoefs2 and then evaluate the error
c   between solncoefs1 and the projected version of solncoefs1
c
      implicit real *8 (a-h,o-z)
      integer nch1,k1,npts1,nch2,k2,npts2,k3,nch3
      complex *16 solncoefs1(npts1),solncoefs2(npts2)
      real *8 erra
      integer, allocatable :: adjs1(:,:),adjs2(:,:),adjs3(:,:),
     1  adjs22(:,:)
      real *8, allocatable :: srcinfo1(:,:),srcinfo2(:,:),srcinfo3(:,:),
     1   srcinfo22(:,:)
      real *8, allocatable :: srccoefs1(:,:),srccoefs2(:,:),
     1  srccoefs3(:,:),srccoefs22(:,:)
      real *8, allocatable :: ts1(:),ts2(:),ts3(:),ts22(:)
      real *8, allocatable :: qwts1(:),qwts2(:),qwts3(:),qwts22(:)
      integer, allocatable :: norders1(:),norders2(:),norders3(:),
     1   norders22(:)
      integer, allocatable :: iptype1(:),iptype2(:),iptype3(:),
     1   iptype22(:)
      integer, allocatable :: ixys1(:),ixys2(:),ixys3(:),ixys22(:)
      real *8, allocatable :: srcrad1(:),srcrad22(:)

      real *8, allocatable :: ts_interp1(:),ts_interp2(:)
      integer, allocatable :: ich_interp1(:),ich_interp2(:)
      real *8, allocatable :: dist1(:),dist2(:)

      real *8, allocatable :: ts_interp12(:)
      integer, allocatable :: ich_interp12(:)
      real *8, allocatable :: dist12(:)

      complex *16, allocatable :: soln1(:),soln2(:),soln12(:)
      complex *16, allocatable :: soln12coefs(:),soln123(:)
      real *8, allocatable :: t(:),w(:),umat(:,:),vmat(:,:)

      real *8 dpars(2),timeinfo(3)
      complex *16 zpars
      integer ipars

      external circ_geom

      done = 1.0d0
      pi = atan(done)*4

      print *, "ncomp=",ncomp
      print *, "nch30=",nch30

      nch3 = nch30*ncomp*4
      npts3 = k3*nch3
      allocate(ts_interp1(npts3),ich_interp1(npts3))
      allocate(ts_interp2(npts3),ich_interp2(npts3))
      allocate(dist1(npts3),dist2(npts3))
      allocate(soln1(npts3),soln2(npts3))

      nch1 = nch10*ncomp*4
      npts1 = nch1*k1

      allocate(adjs1(2,nch1),srcinfo1(8,npts1),srccoefs1(6,npts1))
      allocate(ts1(npts1),qwts1(npts1),norders1(nch1),iptype1(nch1))
      allocate(ixys1(nch1+1),srcrad1(npts1))

      nch2 = nch20*ncomp*4
      npts2 = nch2*k2

      allocate(adjs2(2,nch2),srcinfo2(8,npts2),srccoefs2(6,npts2))
      allocate(ts2(npts2),qwts2(npts2),
     1    norders2(nch2),iptype2(nch2))
      allocate(ixys2(nch2+1))

      k22 = 30
      npts22 = nch2*k22

      allocate(adjs22(2,nch2),srcinfo22(8,npts22),srccoefs22(6,npts22))
      allocate(ts22(npts22),qwts22(npts22),
     1    norders22(nch2),iptype22(nch2))
      allocate(ixys22(nch2+1),srcrad22(npts22))

      allocate(adjs3(2,nch3),srcinfo3(8,npts3),srccoefs3(6,npts3))
      allocate(ts3(npts3),qwts3(npts3),norders3(nch3),iptype3(nch3))
      allocate(ixys3(nch3+1))
      
      ndd_curv = 2
      ndz_curv = 0
      ndi_curv = 0
      call get_diamond_many(nch10,ncomp,rsc,shifts,nch1,k1,npts1,
     1  adjs1,srcinfo1,srccoefs1,qwts1,norders1,iptype1,ixys1)
      
      call get_diamond_many(nch20,ncomp,rsc,shifts,nch2,k2,npts2,
     1  adjs2,srcinfo2,srccoefs2,qwts2,norders2,iptype2,ixys2)

      call get_diamond_many(nch20,ncomp,rsc,shifts,nch2,k22,npts22,
     1  adjs22,srcinfo22,srccoefs22,qwts22,norders22,iptype22,ixys22)

      call get_diamond_many(nch30,ncomp,rsc,shifts,nch3,k3,npts3,
     1  adjs3,srcinfo3,srccoefs3,qwts3,norders3,iptype3,ixys3)

      call prin2('srcinfo3=*',srcinfo3,24)
      call prin2('srcinfo2=*',srcinfo2,24)
      call prin2('srcinfo1=*',srcinfo1,24)
      call prin2('srcinfo22=*',srcinfo22,24)

      do i=1,npts1
        srcrad1(i) = 0
      enddo

      do i=1,npts22
        srcrad22(i) = 0
      enddo


cc      call prin2('srccoefs1=*',srccoefs1,6*npts1)
      
      call findnearchunktarg_id_ts(nch1,norders1,ixys1,iptype1,
     1  npts1,srccoefs1,srcinfo1,srcrad1,8,npts3,srcinfo3,ich_interp1,
     2  ts_interp1,dist1,timeinfo,ier)

      
      call findnearchunktarg_id_ts(nch2,norders22,ixys22,iptype22,
     1  npts22,srccoefs22,srcinfo22,srcrad22,8,npts3,srcinfo3,
     2  ich_interp2,ts_interp2,dist2,timeinfo,ier)

c
c  
c
      allocate(ich_interp12(npts2),ts_interp12(npts2),dist12(npts2))
c      call findnearchunktarg_id_ts_brute(nch1,norders1,ixys1,iptype1,
c     1  npts1,srccoefs1,srcinfo1,8,npts2,srcinfo2,ich_interp12,
c     2  ts_interp12,dist12)
      
      call findnearchunktarg_id_ts(nch1,norders1,ixys1,iptype1,
     1  npts1,srccoefs1,srcinfo1,srcrad1,8,npts2,srcinfo2,ich_interp12,
     2  ts_interp12,dist12,timeinfo,ier)
      allocate(soln12(npts2),soln12coefs(npts2))
      
      call interp_dens(nch1,norders1,ixys1,iptype1,npts1,
     1  solncoefs1,npts2,ts_interp12,ich_interp12,soln12)
c
c
      allocate(t(k2),w(k2),umat(k2,k2),vmat(k2,k2))
      itype = 2
      call legeexps(itype,k2,t,umat,vmat,w)
c
c  convert soln12 to soln12 coefs
c
      alpha = 1.0d0
      beta = 0.0d0
      do i=1,nch2
        istart = ixys2(i)
        call dgemm('n','t',2,k2,k2,alpha,soln12(istart),2,
     1    umat,k2,beta,soln12coefs(istart),2)
      enddo

      call interp_dens(nch1,norders1,ixys1,iptype1,npts1,
     1  solncoefs1,npts3,ts_interp1,ich_interp1,soln1)
c
      call interp_dens(nch2,norders2,ixys2,iptype2,npts2,
     1  solncoefs2,npts3,ts_interp2,ich_interp2,soln2)
      
c
c
      allocate(soln123(npts3))
      call interp_dens(nch2,norders2,ixys2,iptype2,npts2,
     1  soln12coefs,npts3,ts_interp2,ich_interp2,soln123)
      
      erra = 0
      errq = 0
      ra = 0
      do i=1,npts3
        erra = erra + abs(soln1(i)-soln2(i))**2*qwts3(i)
        errq = errq + abs(soln1(i)-soln123(i))**2*qwts3(i)
        ra = ra + abs(soln1(i))**2*qwts3(i)
      enddo
      erra = sqrt(erra/ra)
      errq = sqrt(errq/ra)
      errq = errq/erra

      if(ifwrite.eq.1) then 
        ra = 0
        do i=1,npts3 
          ra = ra + qwts3(i)
          write(iunit,*) ra,real(soln1(i)),imag(soln1(i)),
     1     real(soln2(i)),imag(soln2(i))
        enddo
      endif

      return
      end









      subroutine get_circ_dens_error2(dpars,nch1,k1,npts1,solncoefs1,
     1   nch2,k2,npts2,solncoefs2,k3,nch3,erra,errq)
c
c   Given two densities defined on two grids, and a third
c   reference grid to compute the error, evalute
c   the L2 error in the density
c
c   and errq is the optimality factor, i.e. project solncoefs1, 
c   onto the basis of solncoefs2 and then evaluate the error
c   between solncoefs1 and the projected version of solncoefs1
c
      implicit real *8 (a-h,o-z)
      integer nch1,k1,npts1,nch2,k2,npts2,k3,nch3
      complex *16 solncoefs1(npts1),solncoefs2(npts2)
      real *8 erra
      integer, allocatable :: adjs1(:,:),adjs2(:,:),adjs3(:,:),
     1  adjs22(:,:)
      real *8, allocatable :: srcinfo1(:,:),srcinfo2(:,:),srcinfo3(:,:),
     1   srcinfo22(:,:)
      real *8, allocatable :: srccoefs1(:,:),srccoefs2(:,:),
     1  srccoefs3(:,:),srccoefs22(:,:)
      real *8, allocatable :: ts1(:),ts2(:),ts3(:),ts22(:)
      real *8, allocatable :: qwts1(:),qwts2(:),qwts3(:),qwts22(:)
      integer, allocatable :: norders1(:),norders2(:),norders3(:),
     1   norders22(:)
      integer, allocatable :: iptype1(:),iptype2(:),iptype3(:),
     1   iptype22(:)
      integer, allocatable :: ixys1(:),ixys2(:),ixys3(:),ixys22(:)
      real *8, allocatable :: srcrad1(:),srcrad22(:)

      real *8, allocatable :: ts_interp1(:),ts_interp2(:)
      integer, allocatable :: ich_interp1(:),ich_interp2(:)
      real *8, allocatable :: dist1(:),dist2(:)

      real *8, allocatable :: ts_interp12(:)
      integer, allocatable :: ich_interp12(:)
      real *8, allocatable :: dist12(:)

      complex *16, allocatable :: soln1(:),soln2(:),soln12(:)
      complex *16, allocatable :: soln12coefs(:),soln123(:)
      real *8, allocatable :: t(:),w(:),umat(:,:),vmat(:,:)

      real *8 dpars(2),timeinfo(3)
      complex *16 zpars
      integer ipars

      external circ_geom

      done = 1.0d0
      pi = atan(done)*4

      npts3 = k3*nch3
      allocate(ts_interp1(npts3),ich_interp1(npts3))
      allocate(ts_interp2(npts3),ich_interp2(npts3))
      allocate(dist1(npts3),dist2(npts3))
      allocate(soln1(npts3),soln2(npts3))


      allocate(adjs1(2,nch1),srcinfo1(8,npts1),srccoefs1(6,npts1))
      allocate(ts1(npts1),qwts1(npts1),norders1(nch1),iptype1(nch1))
      allocate(ixys1(nch1+1),srcrad1(npts1))



      allocate(adjs2(2,nch2),srcinfo2(8,npts2),srccoefs2(6,npts2))
      allocate(ts2(npts2),qwts2(npts2),
     1    norders2(nch2),iptype2(nch2))
      allocate(ixys2(nch2+1))

      k22 = 30
      npts22 = nch2*k22

      allocate(adjs22(2,nch2),srcinfo22(8,npts22),srccoefs22(6,npts22))
      allocate(ts22(npts22),qwts22(npts22),
     1    norders22(nch2),iptype22(nch2))
      allocate(ixys22(nch2+1),srcrad22(npts22))

      allocate(adjs3(2,nch3),srcinfo3(8,npts3),srccoefs3(6,npts3))
      allocate(ts3(npts3),qwts3(npts3),norders3(nch3),iptype3(nch3))
      allocate(ixys3(nch3+1))
      
      ndd_curv = 2
      ndz_curv = 0
      ndi_curv = 0

      a = 0.0d0
      b = 2*pi
      call prin2('dpars=*',dpars,2)
      call get_funcurv_geom_uni(a,b,nch1,k1,npts1,adjs1,
     1  srcinfo1,srccoefs1,ts1,qwts1,norders1,iptype1,ixys1,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      
      call get_funcurv_geom_uni(a,b,nch2,k2,npts2,adjs2,
     1  srcinfo2,srccoefs2,ts2,qwts2,norders2,iptype2,ixys2,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)

      
      call get_funcurv_geom_uni(a,b,nch2,k22,npts22,adjs22,
     1  srcinfo22,srccoefs22,ts22,qwts22,norders22,iptype22,ixys22,
     1  circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)

      call get_funcurv_geom_uni(a,b,nch3,k3,npts3,adjs3,
     1  srcinfo3,srccoefs3,ts3,qwts3,norders3,iptype3,ixys3,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      
      do i=1,npts1
        srcrad1(i) = 0
      enddo

      do i=1,npts22
        srcrad22(i) = 0
      enddo


cc      call prin2('srccoefs1=*',srccoefs1,6*npts1)
      
      
      call findnearchunktarg_id_ts(nch1,norders1,ixys1,iptype1,
     1  npts1,srccoefs1,srcinfo1,srcrad1,8,npts3,srcinfo3,ich_interp1,
     2  ts_interp1,dist1,timeinfo,ier)
      
      call findnearchunktarg_id_ts(nch2,norders22,ixys22,iptype22,
     1  npts22,srccoefs22,srcinfo22,srcrad22,8,npts3,srcinfo3,
     2  ich_interp2,ts_interp2,dist2,timeinfo,ier)

c
c  
c
      allocate(ich_interp12(npts2),ts_interp12(npts2),dist12(npts2))
c      call findnearchunktarg_id_ts_brute(nch1,norders1,ixys1,iptype1,
c     1  npts1,srccoefs1,srcinfo1,8,npts2,srcinfo2,ich_interp12,
c     2  ts_interp12,dist12)
      
      call findnearchunktarg_id_ts(nch1,norders1,ixys1,iptype1,
     1  npts1,srccoefs1,srcinfo1,srcrad1,8,npts2,srcinfo2,ich_interp12,
     2  ts_interp12,dist12,timeinfo,ier)
      allocate(soln12(npts2),soln12coefs(npts2))
      
      call interp_dens(nch1,norders1,ixys1,iptype1,npts1,
     1  solncoefs1,npts2,ts_interp12,ich_interp12,soln12)
c
c
      allocate(t(k2),w(k2),umat(k2,k2),vmat(k2,k2))
      itype = 2
      call legeexps(itype,k2,t,umat,vmat,w)
c
c  convert soln12 to soln12 coefs
c
      alpha = 1.0d0
      beta = 0.0d0
      do i=1,nch2
        istart = ixys2(i)
        call dgemm('n','t',2,k2,k2,alpha,soln12(istart),2,
     1    umat,k2,beta,soln12coefs(istart),2)
      enddo

      call interp_dens(nch1,norders1,ixys1,iptype1,npts1,
     1  solncoefs1,npts3,ts_interp1,ich_interp1,soln1)
c
      call interp_dens(nch2,norders2,ixys2,iptype2,npts2,
     1  solncoefs2,npts3,ts_interp2,ich_interp2,soln2)
      
c
c
      allocate(soln123(npts3))
      call interp_dens(nch2,norders2,ixys2,iptype2,npts2,
     1  soln12coefs,npts3,ts_interp2,ich_interp2,soln123)
      
      erra = 0
      errq = 0
      ra = 0
      do i=1,npts3
        erra = erra + abs(soln1(i)-soln2(i))**2*qwts3(i)
        errq = errq + abs(soln1(i)-soln123(i))**2*qwts3(i)
        ra = ra + abs(soln1(i))**2*qwts3(i)
      enddo
      erra = sqrt(erra)
      errq = sqrt(errq)
      errq = errq/erra



      return
      end







      subroutine load_cavity_zpars(a,b,m,zpars,ndz)
      implicit real *8 (a-h,o-z)
      complex *16 zpars(2*m)
      real *8, allocatable :: s(:),r(:),th(:),rho(:),z(:)
      complex *16, allocatable :: zwork(:)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/

      done = 1
      pi = atan(done)*4

      ndz = 2*m
      lw = 10*ndz + 15
      print *, "ndz=",ndz
      allocate(s(m),r(m),th(m),rho(m),z(m))
      allocate(zwork(lw))
      call zffti(ndz,zwork)
      c = a
      hh = c/sqrt(2.0d0)
      aa = 1.0d0
      bb = 0.0d0
      do i=1,m
        s(i) = (i-0.5d0)/(m+0.0d0)*pi
        suse = (s(i)-pi/2)/a
        call qerrfun(suse,rr)
        r(i) = 1-a*rr
        suse = (s(i)-pi/2)
        call crn_fconvgauss(suse, aa, bb, hh, rr, tmp, tmp2)
        th(i) = (b-a) + 2*(1.0d0-(b-a)/pi)*rr
        rho(i) = r(i)*sin(th(i))
        z(i) = r(i)*cos(th(i))*1.2d0

      enddo

      do i=1,m
        zpars(i) = rho(i) + ima*z(i) 
        zpars(i+m) = -rho(m-i+1) + ima*z(m-i+1)
      enddo

      call zfftf(ndz,zpars,zwork)

      call prin2('zpars=*',zpars,8)
      do i=1,ndz
        zpars(i) = zpars(i)/(m+0.0d0)
      enddo


      return
      end




      subroutine funcurv_zfft(t,ndd,dpars,ndz,zpars,ndi,ipars,
     1       srcinfo)
      implicit real *8 (a-h,o-z)
      real *8 dpars(ndd),srcinfo(8),t,pi,done
      complex *16 zpars(ndz),z,zp,zpp,ima
      integer ndd,ndz,ndi,iparss(ndi)
      data ima/(0.0d0,1.0d0)/

      done = 1.0d0
      pi = atan(done)*4
      m = ndz/2
      z = 0
      zp = 0
      zpp = 0
      do i=0,m
        z = z + zpars(i+1)*exp(ima*i*t)
        zp = zp + zpars(i+1)*exp(ima*i*t)*ima*(i+0.0d0)
        zpp = zpp - zpars(i+1)*exp(ima*i*t)*(i+0.0d0)**2
      enddo
      do i=-m,-1
        z = z + zpars(i+1+ndz)*exp(ima*i*t)
        zp = zp + zpars(i+1+ndz)*exp(ima*i*t)*ima*(i+0.0d0)
        zpp = zpp - zpars(i+1+ndz)*exp(ima*i*t)*(i+0.0d0)**2
      enddo

      srcinfo(1) = real(z)
      srcinfo(2) = imag(z)
      srcinfo(3) = real(zp)
      srcinfo(4) = imag(zp)
      srcinfo(5) = real(zpp)
      srcinfo(6) = imag(zpp)

      ds = sqrt(srcinfo(3)**2 + srcinfo(4)**2)
      srcinfo(7) = srcinfo(4)/ds
      srcinfo(8) = -srcinfo(3)/ds


      return
      end

