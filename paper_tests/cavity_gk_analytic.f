      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),srccoefs(:,:),qwts(:),ts1(:)
      integer, allocatable :: norders(:),iptype(:),ixys(:)
      integer, allocatable :: adjs(:,:)
      real *8 xypt(2),shifts(2)
      complex *16, allocatable :: zpars_curv(:)

      complex *16 zpars(3)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:),wnearcoefs(:)
      complex *16, allocatable :: sigma(:),sigmacoefs(:),sigmatest(:)
      complex *16, allocatable :: soln(:),solncoefs(:)

      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: srcrad(:)
      real *8, allocatable :: pmat(:,:)
      real *8 srctmp(8)

      real *8 xyin(2),xyout(2),xyuse(2),timeinfo(3)

      integer, allocatable :: ich_id_targ(:)
      real *8, allocatable :: targs(:,:),ts_targ(:),potplot(:)

      complex *16, allocatable :: pottarg_plot(:)
      real *8, allocatable :: errs(:)
      real *8, allocatable :: pols(:)
      complex *16 zk,potex,pottarg,ztmp,pottarg2

      external funcurv_zfft

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4


      nch = 50
      k = 16
      ngk = 4
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
      
c
c  get density info
c
      alpha = 1.0d0
      beta = 0.0d0

      xyin(1) = 0.01d0
      xyin(2) = -1.1d0

      xyout(1) = 6.1d0
      xyout(2) = 3.2d0

      npts_gk = nch*ngk 
      allocate(sigma(npts), sigmacoefs(npts_gk), sigmatest(npts))
      allocate(soln(npts), solncoefs(npts_gk), pols(ngk))

      erra = 0
      ra = 0
      do ich=1,nch
        do j=1,k
          ipt = (ich-1)*k + j
          call h2d_slp(srcinfo(1,ipt),2,xyin,ndd,dpars,1,zk,ndi,
     1       ipars,sigma(ipt))
          soln(ipt) = 0
        enddo

        do j=1,ngk
          ipt = (ich-1)*ngk + j
          sigmacoefs(ipt) = 0
          solncoefs(ipt) = 0
        enddo

        
        istart = (ich-1)*k + 1
        jstart = (ich-1)*ngk + 1
        call get_galerkin_proj(2,ngk,k,ts,sigma(istart),qwts(istart), 
     1      sigmacoefs(jstart))
        do j=1,k
          call legepols(ts(j),ngk-1,pols)
          jpt = (ich-1)*k + j
          sigmatest(jpt) = 0
          do l=1,ngk
            lpt = (ich-1)*ngk + l
            sigmatest(jpt) = sigmatest(jpt) + pols(l)*sigmacoefs(lpt)
          enddo
          erra = erra + abs(sigmatest(jpt) - sigma(jpt))**2*qwts(jpt)
          ra = ra + abs(sigma(jpt))**2*qwts(jpt)
        enddo
      enddo

      call prin2('sigmacoefs=*',sigmacoefs,2*ngk)
      erra = sqrt(erra/ra)
      call prin2('error in projector=*',erra,1)
c
c
c  now test layer potential evaluator and compare with existing test
c


      nb = ngk*nch
      nnz = 3*ngk*nch
      nquad = nnz*ngk
      allocate(row_ptr(nb+1),col_ind(nnz))
      allocate(wnear(nquad))

      call prinf('ngk=*',ngk,1)
      call prinf('nb=*',nb,1)

      ra = 0
      do i=1,npts
        ra = ra + qwts(i)
      enddo
      call prin2('length of curve=*',ra,1)

      call get_helm_dir_trid_quad_corr(zk,nch,k,ngk,npts,nb,adjs,
     1 srcinfo,srccoefs,ndz,zpars,nnz,row_ptr,col_ind,nquad,wnear)
      
      call prinf('nnz=*',nnz,1)
      call prinf('nb=*',nb,1)
      call prinf('row_ptr=*',row_ptr,nb+1)
      call prinf('col_ind=*',col_ind,nnz)

      call prin2('wnear=*',wnear,24)

      allocate(iquad(nnz+1))

      do i=1,nnz+1
        iquad(i) = (i-1)*ngk+1
      enddo

      eps = 0.51d-12
      niter = 0
      numit = 300
      ifinout = 1
      allocate(errs(numit+1))
      call helm_comb_dir_galerkin_solver2d(nch,k,ngk,ixys,npts,nb,
     1  srcinfo,adjs,qwts,eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,
     2  wnear,norders,npts,ixys,srcinfo,qwts,numit,ifinout,
     3  sigmacoefs,eps,niter,errs,rres,solncoefs)
       
      call prin2('solncoefs=*',solncoefs,24)

      potex = 0
      call h2d_slp(xyout,2,xyin,ndd,dpars,1,zk,ndi,ipars,potex)
      allocate(pmat(ngk,k))
      do i=1,k
        call legepols(ts(i),ngk-1,pmat(1,i))
      enddo

      do ich=1,nch
        istart = (ich-1)*ngk+1
        jstart = (ich-1)*k+1
        call dgemm('n','n',2,k,ngk,alpha,solncoefs(istart),2,pmat,
     1   ngk,beta,soln(jstart),2)
      enddo


      pottarg = 0
      pottarg2 = 0
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
