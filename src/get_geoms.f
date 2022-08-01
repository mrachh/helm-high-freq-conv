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

      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      do ich=1,nch
        tstart = a+(ich-1.0d0)*h
        tend = a+(ich+0.0d0)*h
c
c  get source info
c

        do j=1,k
          ipt = (ich-1)*k + j
          tuse = tstart + (tend-tstart)*(ts(j)+1)/2
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

