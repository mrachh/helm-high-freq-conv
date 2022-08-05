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
      real *8 shifts(2,100),rsc(100)
      complex *16 zk,zpars(3),ima,z1,z2,ztmp
      complex *16 pottarg,pottargex
      real *8 xyin(2),xyout(2)
      real *8 dpars(1)
      real *8, allocatable :: errs(:)
      complex *16 zfac
      real *8 xy_in(2),xy_out(2)
      data ima/(0.0d0,1.0d0)/
      external circ_geom

      

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

      zk = 100.0d0 + 0.0d0*ima
      zk = 10.0d0*sqrt(2.0d0)
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

      rfac = sqrt(2.0d0)
      ncomp = 4
      rsc(1) = 0.8d0*rfac*pi/2
      rsc(2) = 0.8d0*rfac*pi/2
      rsc(3) = 0.8d0*rfac*pi/2
      rsc(4) = 0.8d0*rfac*pi/2

      shifts(1,1) = rfac*pi
      shifts(2,1) = rfac*pi

      shifts(1,2) = rfac*pi
      shifts(2,2) = -rfac*pi

      shifts(1,3) = -rfac*pi
      shifts(2,3) = rfac*pi

      shifts(1,4) = -rfac*pi
      shifts(2,4) = -rfac*pi

      nch0 = ceiling(0.4*4*abs(zk)/sqrt(2.0d0))
      nch = 4*nch0*ncomp
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

      print *, "nch0=",nch0
      print *, "ncomp=",ncomp
      print *, "nch=",nch
      call get_diamond_many(nch0,ncomp,rsc,shifts,nch,
     1 k,npts,adjs,srcinfo,srccoefs,
     1 qwts,norders,iptype,ixys)
      call get_diamond_many(nch0,ncomp,rsc,shifts,nch,
     1 kg,nptsg,adjs,srcinfog,srccoefsg,
     1 qwtsg,nordersg,iptype,ixysg)
      call get_diamond_many(nch0,ncomp,rsc,shifts,nch,
     1 nover,npts_over,adjs,srcover,
     1 srccoefsover,wover,novers,iptype,ixyso)
      




      print *, "nptsg=",nptsg
      allocate(sigmag(nptsg),sigmacoefsg(nptsg),solncoefsg(nptsg))
      allocate(solng(nptsg))

      xyin(1) = 0.1d0*rsc(1) + shifts(1,1)
      xyin(2) = -0.33d0*rsc(1) + shifts(2,1)

      xyout(1) = 3.5d1
      xyout(2) = 3.3d1

c
c  get density info
c
      thet = pi/4
      do ich=1,nch
        do j=1,kg
          ipt = (ich-1)*kg + j
          call h2d_slp(srcinfog(1,ipt),2,xyin,ndd,dpars,1,zk,ndi,
     1       ipars,sigmag(ipt))
c          sigmag(ipt) = exp(ima*zk*(srcinfog(1,ipt)*cos(thet)+ 
c     1         srcinfog(2,ipt)*sin(thet)))
          solncoefsg(ipt) = 0.0d0
        enddo
        istart = (ich-1)*kg + 1
        call dgemm('n','t',2,kg,kg,alpha,sigmag(istart),
     1     2,umatg,kg,beta,sigmacoefsg(istart),2)
      enddo
      ra = sum(wover)
      ra2 = sum(qwts)
      print *, "ra=",ra
      print *, "ra2=",ra2
      print *, "diff=",abs(ra2-ra)

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
     1   srcinfo,srcinfog,ndz,zpars,nnzg,row_ptrg,col_indg,
     2   nquadg,wnearg,wnearcoefsg)
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
