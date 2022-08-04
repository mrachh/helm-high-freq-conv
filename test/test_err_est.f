      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcinfog(:,:),qwtsg(:),srccoefsg(:,:)
      real *8, allocatable :: srcover(:,:),wover(:),srccoefsover(:,:)
      real *8, allocatable :: ts1(:),ts1g(:),tsover(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tsg(:),umatg(:,:),vmatg(:,:),wtsg(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8 umo,vmo
      complex *16, allocatable :: sigma(:),sigmag(:)
      complex *16, allocatable :: sigmacoefs(:),sigmacoefsg(:)

      integer, allocatable :: adjs(:,:),adjsg(:,:),adjso(:,:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),ixysg(:)
      integer, allocatable :: iptypeg(:),iptypeo(:)
      integer, allocatable :: nordersg(:)
      integer, allocatable :: novers(:),ixyso(:)
      complex *16 zk,zpars(3),ima,z1,z2
      complex *16 fjvals(0:100),fhvals(0:100),fjders(0:100)
      real *8 dpars(2)
      complex *16 fhders(0:100),zfac
      real *8 xy_in(2),xy_out(2)
      data ima/(0.0d0,1.0d0)/
      external circ_geom

      

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

      zk = 100.0d0 + 0.0d0*ima
      zk = 1.0d0
      zpars(1) = zk
      zpars(2) = -ima*zk
      zpars(3) = 1.0d0

      imode = 7

      alpha = 1.0d0
      beta = 0.0d0
c
c  k is the boundary discretization order
c
c  kg is the galerkin polynomial order
c
c  nover is the oversampled source order
c

      k = 22
      kg = 16
      nover = 30
      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      allocate(tsg(kg),umatg(kg,kg),vmatg(kg,kg),wtsg(kg))
      call legeexps(itype,kg,tsg,umatg,vmatg,wtsg)


      nch = 10
      nchg = 8
      ncho = 12
      npts = nch*k
      nptsg = nchg*kg
      npts_over = ncho*nover
      allocate(srcinfo(8,npts),qwts(npts))
      allocate(srccoefs(6,npts),ts1(npts))
      allocate(norders(nch),iptype(nch),ixys(nch+1))
      allocate(adjs(2,nch),adjsg(2,nchg),adjso(2,ncho))


      allocate(srcinfog(8,nptsg),srccoefsg(6,nptsg))
      allocate(qwtsg(nptsg))
      allocate(nordersg(nchg),ixysg(nch+1),ts1g(nptsg))
      allocate(iptypeg(nchg))

      allocate(srccoefsover(6,npts_over))
      allocate(srcover(8,npts_over),wover(npts_over))
      allocate(novers(ncho),ixyso(ncho+1))
      allocate(tsover(npts_over),iptypeo(ncho))



      a = 0.0d0
      b = 2*pi
      dpars(1) = 1.0d0
      dpars(2) = 1.0d0
      ndd_curv = 2
      ndz_curv = 0
      ndi_curv = 0
      call get_funcurv_geom_uni(a,b,nch,k,npts,adjs,srcinfo,
     1  srccoefs,ts1,qwts,norders,iptype,ixys,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      print *, "nch=",nch
      print *, "kg=",kg
      print *, "nptsg=",nptsg
      call get_funcurv_geom_uni(a,b,nchg,kg,nptsg,adjsg,srcinfog,
     1  srccoefsg,ts1g,qwtsg,nordersg,iptypeg,ixysg,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      call prin2('a=*',a,1)
      call prin2('b=*',b,1)

      print *, "nch=",nch
      print *, "nover=",nover
      print *, "a=",a
      print *, "b=",b
      print *, "npts_over=",npts_over
      call get_funcurv_geom_uni(a,b,ncho,nover,npts_over,adjso,srcover,
     1  srccoefsover,tsover,wover,novers,iptypeo,ixyso,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)



      h = 2*pi/(nch+0.0d0)
      print *, "nptsg=",nptsg
      allocate(sigma(npts),sigmacoefs(npts))
      allocate(sigmag(nptsg),sigmacoefsg(nptsg))

c
c  get density info
c
      do ich=1,nch

        do j=1,k
          ipt = (ich-1)*k + j
          tuse = ts1(ipt)
          sigma(ipt) = exp(ima*(imode+0.0d0)*tuse)
          write(33,*) tuse,real(sigma(ipt)),imag(sigma(ipt))
        enddo
        istart = (ich-1)*k+1

        call dgemm('n','t',2,k,k,alpha,sigma(istart),
     1     2,umat,k,beta,sigmacoefs(istart),2)
      enddo
      do ich=1,nchg
        do j=1,kg
          ipt = (ich-1)*kg + j
          tuse = ts1g(ipt) 
          sigmag(ipt) = exp(ima*(imode+0.0d0)*tuse)
          write(34,*) tuse,real(sigmag(ipt)),imag(sigmag(ipt))
        enddo
        istart = (ich-1)*kg + 1
        call dgemm('n','t',2,kg,kg,alpha,sigmag(istart),
     1     2,umatg,kg,beta,sigmacoefsg(istart),2)
      enddo

      call get_circ_dens_error(dpars,nch,k,npts,sigmacoefs,
     1  nchg,kg,nptsg,sigmacoefsg,nover,ncho,erra)
      print *, "erra=",erra

       

      return
      end
