      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),srccoefs(:,:),qwts(:),ts1(:)
      integer, allocatable :: norders(:),iptype(:),ixys(:)
      integer, allocatable :: adjs(:,:)
      real *8 xypt(2),shifts(2)
      complex *16, allocatable :: zpars_curv(:)

      complex *16 zpars(3)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:),wnearcoefs(:)
      complex *16, allocatable :: sigma(:),sigmacoefs(:)
      complex *16, allocatable :: soln(:),solncoefs(:)

      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)

      real *8 xyin(2),xyout(2)
      real *8, allocatable :: errs(:)
      complex *16 zk,potex,pottarg,ztmp

      external funcurv_zfft

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4


      nch = 30
c      k = 10
      k = 16
      npts = k*nch

      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      zk = 1.0d0
      ndz = 3
      zpars(1) = zk
      zpars(2) = -ima*zk
      zpars(3) = 1.0d0

      mm = 200
      ndz_curv = 2*mm

      ndd_curv = 0
      ndi_curv = 0

      allocate(zpars_curv(ndz_curv))

      aa = 0.2d0
      bb = pi/12

      call load_cavity_zpars(aa,bb,mm,zpars_curv,ndz_curv)

      allocate(adjs(2,nch),srcinfo(8,npts),srccoefs(6,npts))
      allocate(qwts(npts),norders(nch),iptype(nch),ixys(nch+1))
      allocate(ts1(npts))

      a = 0
      b = 2.0d0*pi
      call get_funcurv_geom_uni(a,b,nch,k,npts,adjs,srcinfo,srccoefs,
     1  ts1,qwts,norders,iptype,ixys,funcurv_zfft,ndd_curv,dpars,
     2  ndz_curv,zpars_curv,ndi_curv,ipars)
      
      ra = 0
      ifwrite = 0
      do i=1,npts
        if(ifwrite.eq.1) then
         write(37,*) srcinfo(1,i),srcinfo(2,i),srcinfo(7,i),srcinfo(8,i)
        endif
        ra = ra + qwts(i)
      enddo

      call prin2('length of curve=*',ra,1)
c
c  get density info
c
      alpha = 1.0d0
      beta = 0.0d0

      xyin(1) = 0.0d0
      xyin(2) = -2.2d0

      xyout(1) = 6.1d0
      xyout(2) = 3.2d0


      allocate(sigma(npts),sigmacoefs(npts))
      allocate(soln(npts),solncoefs(npts))


      do ich=1,nch
        do j=1,k
          ipt = (ich-1)*k + j
          call h2d_slp(srcinfo(1,ipt),2,xyin,ndd,dpars,1,zk,ndi,
     1       ipars,sigma(ipt))
          solncoefs(ipt) = 0
          soln(ipt) = 0
        enddo
        istart = (ich-1)*k + 1
        call dgemm('n','t',2,k,k,alpha,sigma(istart),2,umat,
     1    k,beta,sigmacoefs(istart),2)
      enddo
      

      nnz = 3*k*nch
      nquad = nnz*k
      allocate(row_ptr(npts+1),col_ind(nnz),wnear(nquad))
      allocate(wnearcoefs(nquad))

      call get_helm_dir_trid_quad_corr(zk,nch,k,k,npts,npts,ads,
     1 srcinfo,srcinfo,ndz,zpars,nnz,row_ptr,col_ind,nquad,wnear,
     2 wnearcoefs)

      allocate(iquad(nnz+1))
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,iquad)

      eps = 0.51d-12
      niter = 0
      numit = 300
      ifinout = 1
      allocate(errs(numit+1))
      call helm_comb_dir_galerkin_solver2d(nch,k,ixys,npts,
     1  srcinfo,eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,
     2  wnearcoefs,norders,npts,ixys,srcinfo,qwts,numit,ifinout,
     3  sigmacoefs,eps,niter,errs,rres,solncoefs)
      
      potex = 0
      call h2d_slp(xyout,2,xyin,ndd,dpars,1,zk,ndi,ipars,potex)
      
      do ich=1,nch
        istart = (ich-1)*k+1
        call dgemm('n','t',2,k,k,alpha,solncoefs(istart),2,vmat,
     1   k,beta,soln(istart),2)
      enddo

      pottarg = 0
      do i=1,npts
        call h2d_comb(srcinfo(1,i),2,xyout,ndd,dpars,ndz,zpars,
     1    ndi,ipars,ztmp)
        pottarg = pottarg + ztmp*qwts(i)*soln(i)
      enddo
      call prin2_long('pottarg=*',pottarg,2)
      call prin2_long('potex=*',potex,2)

      erra = abs(potex-pottarg)/abs(potex)
      call prin2('relative error=*',erra,1)

      


      stop
      end
