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
      complex *16, allocatable :: sigmacoefspw(:),sigmapw(:)
      complex *16, allocatable :: solncoefspw(:),solnpw(:)

      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: srcrad(:)
      real *8 srctmp(8)

      real *8 xyin(2),xyuse(2),timeinfo(3)
      real *8, allocatable :: xyout(:,:)

      integer, allocatable :: ich_id_targ(:)
      real *8, allocatable :: targs(:,:),ts_targ(:)

      complex *16, allocatable :: pottarg(:),pottargex(:),potplot(:)
      complex *16, allocatable :: uinc(:)

      real *8, allocatable :: errs(:)
      complex *16 zk,ztmp,pottarg2,ima

      data ima/(0.0d0,1.0d0)/

      external funcurv_zfft

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4



      open(unit=33,file='input_data.dat')
      open(unit=34,file='target_data.dat')

      read(33,*) tmp
      zk = tmp
      read(33,*) nch
      read(33,*) ntarg
      read(33,*) thet
      read(33,*) is_analytic_test

      close(33)

      print *, ntarg

      allocate(targs(2,ntarg),ich_id_targ(ntarg),ts_targ(ntarg))
      allocate(pottarg(ntarg),pottargex(ntarg),potplot(ntarg))
      allocate(uinc(ntarg))

      do i=1,ntarg
        read(34,*) targs(1,i),targs(2,i)
      enddo
      close(34)


      k = 16
      npts = k*nch

      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

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


      allocate(sigma(npts),sigmacoefs(npts))
      allocate(soln(npts),solncoefs(npts))

      allocate(sigmacoefspw(npts),solncoefspw(npts))
      allocate(solnpw(npts),sigmapw(npts))


      print *, "thet=",thet
      ct = cos(thet)
      st = sin(thet)
      do ich=1,nch
        do j=1,k
          ipt = (ich-1)*k + j
          call h2d_slp(srcinfo(1,ipt),2,xyin,ndd,dpars,1,zk,ndi,
     1       ipars,sigma(ipt))
          sigmapw(ipt) = -exp(ima*zk*(srcinfo(1,ipt)*ct +
     1        srcinfo(2,ipt)*st))
          solncoefs(ipt) = 0
          soln(ipt) = 0
          solnpw(ipt) = 0
          solncoefspw(ipt) = 0
        enddo
        istart = (ich-1)*k + 1
        call dgemm('n','t',2,k,k,alpha,sigma(istart),2,umat,
     1    k,beta,sigmacoefs(istart),2)
        call dgemm('n','t',2,k,k,alpha,sigmapw(istart),2,umat,
     1    k,beta,sigmacoefspw(istart),2)
      enddo

      call prin2('sigmacoefs=*',sigmacoefs,2*k)


      nnz = 3*k*nch
      nquad = nnz*k
      allocate(row_ptr(npts+1),col_ind(nnz),wnear(nquad))
      allocate(wnearcoefs(nquad))

      call get_helm_dir_trid_quad_corr(zk,nch,k,k,npts,npts,adjs,
     1 srcinfo,srcinfo,ndz,zpars,nnz,row_ptr,col_ind,nquad,wnear,
     2 wnearcoefs)

      allocate(iquad(nnz+1))
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,iquad)

      eps = 0.51d-12
      niter = 0
      numit = 300
      ifinout = 1
      allocate(errs(numit+1))

      if(is_analytic_test.eq.1) then
        call helm_comb_dir_galerkin_solver2d(nch,k,ixys,npts,
     1    srcinfo,eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,
     2    wnearcoefs,norders,npts,ixys,srcinfo,qwts,numit,ifinout,
     3    sigmacoefs,eps,niter,errs,rres,solncoefs)
          do ich=1,nch
            istart = (ich-1)*k+1
            call dgemm('n','t',2,k,k,alpha,solncoefs(istart),2,vmat,
     1        k,beta,soln(istart),2)
          enddo

      endif
      call helm_comb_dir_galerkin_solver2d(nch,k,ixys,npts,
     1    srcinfo,eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,
     2    wnearcoefs,norders,npts,ixys,srcinfo,qwts,numit,ifinout,
     3    sigmacoefspw,eps,niter,errs,rres,solncoefspw)
       

      do ich=1,nch
        istart = (ich-1)*k+1
        call dgemm('n','t',2,k,k,alpha,solncoefspw(istart),2,vmat,
     1   k,beta,solnpw(istart),2)
      enddo

      call prin2('soln=*',soln,2*k)
      do i=1,ntarg
        ich_id_targ(i) = -1
        ts_targ(i) = 0
      enddo

      eps_plot = 1.0d-10
      if(is_analytic_test.eq.1) then
        do i=1,ntarg
          call h2d_slp(targs(1,i),2,xyin,ndd,dpars,1,zk,ndi,
     1       ipars,pottargex(i))
        enddo


        call lpcomp_helm_comb_dir_2d(nch,norders,ixys,iptype,npts,
     1    srccoefs,srcinfo,2,ntarg,targs,ich_id_targ,ts_targ,eps_plot,
     2     zpars,soln,pottarg)
        
        erra = 0
        ra = 0
        do i=1,ntarg
           ra = ra + abs(pottargex(i))**2
           erra = erra + abs(pottargex(i)-pottarg(i))**2
        enddo
        erra = sqrt(erra/ra)
        print *, "erra=",erra
      endif

      do i=1,ntarg
        uinc(i) = -exp(ima*zk*(targs(1,i)*ct + targs(2,i)*st))
      enddo
      call lpcomp_helm_comb_dir_2d(nch,norders,ixys,iptype,npts,
     1    srccoefs,srcinfo,2,ntarg,targs,ich_id_targ,ts_targ,eps_plot,
     2    zpars,solnpw,potplot)
      open(unit=35,file='potential.dat',status='replace')
      do i=1,ntarg
        write(35,*) real(potplot(i)),imag(potplot(i)),real(uinc(i)),
     1      imag(uinc(i))
      enddo
      close(35)

      if(is_analytic_test.eq.1) then
        open(unit=36,file='potential_analytic.dat',status='replace')
        do i=1,ntarg
           write(36,*) real(pottarg(i)),imag(pottarg(i)),
     1         real(pottargex(i)),imag(pottargex(i))
        enddo
      endif

        

      stop
      end
