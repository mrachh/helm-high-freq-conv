      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcinfog(:,:),qwtsg(:),srccoefsg(:,:)
      real *8, allocatable :: srcover(:,:),wover(:),srccoefsover(:,:)
      real *8, allocatable :: ts1(:),ts1g(:),tsover(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tsg(:),umatg(:,:),vmatg(:,:),wtsg(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8 umo,vmo
      complex *16, allocatable :: sigmag(:),sigmacoefsg(:),solncoefsg(:)
      complex *16, allocatable :: solng(:)
      integer, allocatable :: row_ptrg(:),col_indg(:),iquadg(:)
      complex *16, allocatable :: wnearg(:),wnearcoefsg(:)
      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_pts(:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),ixysg(:)
      integer, allocatable :: nordersg(:)
      integer, allocatable :: novers(:),ixyso(:),adjs(:,:)
      complex *16 zk,zpars(3),ima,z1,z2,ztmp
      complex *16 pottarg,pottargex
      real *8 xyin(2),xyout(2)
      real *8 dpars(2)
      real *8, allocatable :: errs(:)
      complex *16 zfac
      real *8 xy_in(2),xy_out(2)
      data ima/(0.0d0,1.0d0)/
      external circ_geom

      

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

      zk = 100.0d0 + 0.0d0*ima
      zk = 1.0d0
      ndz = 3
      zpars(1) = zk
      zpars(2) = -ima*zk
      zpars(3) = 1.0d0

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


      nch = 20
      npts = nch*k
      npts_over = nch*nover
      nptsg = nch*kg
      allocate(srcinfo(8,npts),qwts(npts))
      allocate(srccoefs(6,npts),ts1(npts))
      allocate(norders(nch),iptype(nch),ixys(nch+1),adjs(2,nch))


      allocate(srcinfog(8,nptsg),srccoefsg(6,nptsg))
      allocate(qwtsg(nptsg))
      allocate(nordersg(nch),ixysg(nch+1),ts1g(nptsg))

      allocate(srccoefsover(6,npts_over))
      allocate(srcover(8,npts_over),wover(npts_over))
      allocate(novers(nch),ixyso(nch+1))
      allocate(tsover(npts_over))



      a = 0.0d0
      b = 2*pi
      dpars(1) = 1.0d0
      dpars(2) = 1.2d0
      ndd_curv = 2
      ndz_curv = 0
      ndi_curv = 0
      call get_funcurv_geom_uni(a,b,nch,k,npts,adjs,srcinfo,
     1  srccoefs,ts1,qwts,norders,iptype,ixys,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      print *, "nch=",nch
      print *, "kg=",kg
      print *, "nptsg=",nptsg
      call get_funcurv_geom_uni(a,b,nch,kg,nptsg,adjs,srcinfog,
     1  srccoefsg,ts1g,qwtsg,nordersg,iptype,ixysg,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      call prin2('a=*',a,1)
      call prin2('b=*',b,1)

      print *, "nch=",nch
      print *, "nover=",nover
      print *, "a=",a
      print *, "b=",b
      print *, "npts_over=",npts_over
      call get_funcurv_geom_uni(a,b,nch,nover,npts_over,adjs,srcover,
     1  srccoefsover,tsover,wover,novers,iptype,ixyso,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)



      h = 2*pi/(nch+0.0d0)
      print *, "nptsg=",nptsg
      allocate(sigmag(nptsg),sigmacoefsg(nptsg),solncoefsg(nptsg))
      allocate(solng(nptsg))

      xyin(1) = 0.1d0
      xyin(2) = -0.33d0

      xyout(1) = 3.5d0
      xyout(2) = 3.3d0

c
c  get density info
c
      do ich=1,nch
        do j=1,kg
          ipt = (ich-1)*kg + j
          call h2d_slp(srcinfog(1,ipt),2,xyin,ndd,dpars,1,zk,ndi,
     1       ipars,sigmag(ipt))
          solncoefsg(ipt) = 0.0d0
        enddo
        istart = (ich-1)*kg + 1
        call dgemm('n','t',2,kg,kg,alpha,sigmag(istart),
     1     2,umatg,kg,beta,sigmacoefsg(istart),2)
      enddo
      ra = sum(wover)
      ra2 = sum(qwts)
      ra3 = sum(qwtsg)


c
c  test near quad at second patch from first patch
c
c

      nnzg = 3*kg*nch
      nquadg = nnzg*kg

      allocate(row_ptrg(nptsg+1),col_indg(nnzg),wnearg(nquadg))
      allocate(wnearcoefsg(nquadg))

      print *, "starting trid quad"
      print *, "nch=",nch
      print *, "ndz=",ndz
      ndd = 0
      ndi = 0
      ndz = 3
      call prin2('zk=*',zk,2)
      print *, "npts=",npts
      print *, "nptsg=",nptsg
      call cpu_time(t1)
C$       t1 = omp_get_wtime()     
      call get_helm_dir_trid_quad_corr(zk,nch,k,kg,npts,nptsg,adjs,
     1   srcinfo,srcinfog,ndz,zpars,nnzg,row_ptrg,col_indg,nquadg,
     2   wnearg,wnearcoefsg)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      call prin2('total quad gen time=*',t2-t1,1)
      call prinf('nnzg=*',nnzg,1)


      allocate(iquadg(nnzg+1))
      call get_iquad_rsc2d(nch,ixysg,nptsg,nnzg,row_ptrg,col_indg,
     1   iquadg)
      iquadtype = 1
      eps = 0.51d-14

      niter = 0
      numit = max(ceiling(10*abs(zk)),200)
      allocate(errs(numit+1))
      ifinout = 1
      call helm_comb_dir_galerkin_solver2d(nch,kg,ixysg,nptsg,
     1  srcinfog,eps,zpars,nnzg,row_ptrg,col_indg,iquadg,
     2  nquadg,wnearcoefsg,novers(1),npts_over,ixyso,srcover,
     3  wover,numit,ifinout,sigmacoefsg,eps,niter,errs,rres,
     4  solncoefsg)

      call h2d_slp(xyout,2,xyin,ndd,dpars,1,zk,ndi,ipars,pottargex)
      do i=1,nch
        istart = (i-1)*kg+1
        call dgemm('n','t',2,kg,kg,alpha,solncoefsg(istart),
     1     2,vmatg,kg,beta,solng(istart),2)
      enddo

      pottarg = 0
      do i=1,nptsg
        call h2d_comb(srcinfog(1,i),2,xyout,ndd,dpars,ndz,zpars,
     1     ndi,ipars,ztmp)
        pottarg = pottarg + ztmp*qwtsg(i)*solng(i)
      enddo
      call prin2_long('pottarg=*',pottarg,2)
      call prin2_long('pottargex=*',pottargex,2)
      print *, "erra=",abs(pottarg-pottargex)/abs(pottargex)
      
       

      return
      end
