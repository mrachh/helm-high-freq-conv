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

      complex *16, allocatable :: solnwp(:),solncoefswp(:)
      complex *16, allocatable :: sigmawp(:),sigmacoefswp(:)

      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: srcrad(:)
      real *8 srctmp(8)

      real *8 xyin(2),xyout(2),xyuse(2),timeinfo(3)

      integer, allocatable :: ich_id_targ(:)
      real *8, allocatable :: targs(:,:),ts_targ(:),potplot(:)

      complex *16, allocatable :: pottarg_plot(:)
      real *8, allocatable :: errs(:)
      complex *16 zk,potex,pottarg,ztmp

      external funcurv_zfft

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4


      nch = 200
c      k = 10
      k = 16
      npts = k*nch

      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      zk = 40.0d0
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
      ifwrite = 1
      do i=1,npts
        if(ifwrite.eq.1) then
         write(37,*) srcinfo(1,i),srcinfo(2,i),srcinfo(7,i),srcinfo(8,i)
        endif
        ra = ra + qwts(i)
      enddo


      xyuse(1) = -0.799277640571460d0
      xyuse(2) = -0.040786944817011d0

      allocate(srcrad(npts))

      do i=1,npts
        srcrad(i) = 1.0d0
      enddo

      call findnearchunktarg_id_ts(nch,norders,ixys,iptype,npts,
     1   srccoefs,srcinfo,srcrad,2,1,xyuse,ich_id,tst,dist_targ,
     2   timeinfo,ier) 
      call prinf('ich_id=*',ich_id,1)
      call prin2_long('ts_targ=*',tst,1)
      call prin2('dist_targ=*',dist_targ,1)
      call prinf('ier=*',ier,1)

      ttuse = (ich_id-1)*2*pi/nch + 2*pi/nch*(tst+1.0d0)/2
      call funcurv_zfft(ttuse,ndd_curv,dpars,ndz_curv,zpars_curv,
     1  ndi_curv,ipars,srctmp)
      call prin2_long('srctmp=*',srctmp,2)
      call prin2_long('ttuse=*',ttuse,1)


      ttuse = 0.3913220708069514D+01




      call prin2('length of curve=*',ra,1)
c
c  get density info
c
      alpha = 1.0d0
      beta = 0.0d0

      xyin(1) = 0.01d0
      xyin(2) = -1.1d0

      xyout(1) = 6.1d0
      xyout(2) = 3.2d0


      allocate(sigma(npts),sigmacoefs(npts))
      allocate(soln(npts),solncoefs(npts))

      allocate(sigmawp(npts),sigmacoefswp(npts))
      allocate(solnwp(npts),solncoefswp(npts))

      s0 = 1.0d0/sqrt(real(zk))
      thet = 0.45d0*pi

      do ich=1,nch
        do j=1,k
          ipt = (ich-1)*k + j
          call h2d_slp(srcinfo(1,ipt),2,xyin,ndd,dpars,1,zk,ndi,
     1       ipars,sigma(ipt))
          sigmawp(ipt) = exp(-(ts1(ipt)-ttuse)**2/s0**2)*
     1       exp(ima*zk*cos(thet)*ts1(ipt))
          solncoefs(ipt) = 0
          soln(ipt) = 0
          solncoefswp(ipt) = 0
          solnwp(ipt) = 0
        enddo
        istart = (ich-1)*k + 1
        call dgemm('n','t',2,k,k,alpha,sigma(istart),2,umat,
     1    k,beta,sigmacoefs(istart),2)
        call dgemm('n','t',2,k,k,alpha,sigmawp(istart),2,umat,
     1    k,beta,sigmacoefswp(istart),2)
      enddo

      if(ifwrite.eq.1) then
        ra = 0
        do i=1,npts
          ra = ra + qwts(i)
          write(38,*) ra,real(sigmawp(i)),imag(sigmawp(i))

        enddo

      endif
      

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
      call helm_comb_dir_galerkin_solver2d(nch,k,ixys,npts,
     1  srcinfo,eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,
     2  wnearcoefs,norders,npts,ixys,srcinfo,qwts,numit,ifinout,
     3  sigmacoefs,eps,niter,errs,rres,solncoefs)
      
      call helm_comb_dir_galerkin_solver2d(nch,k,ixys,npts,
     1  srcinfo,eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,
     2  wnearcoefs,norders,npts,ixys,srcinfo,qwts,numit,ifinout,
     3  sigmacoefswp,eps,niter,errs,rres,solncoefswp)
      
      potex = 0
      call h2d_slp(xyout,2,xyin,ndd,dpars,1,zk,ndi,ipars,potex)
      
      do ich=1,nch
        istart = (ich-1)*k+1
        call dgemm('n','t',2,k,k,alpha,solncoefs(istart),2,vmat,
     1   k,beta,soln(istart),2)
        call dgemm('n','t',2,k,k,alpha,solncoefswp(istart),2,vmat,
     1   k,beta,solnwp(istart),2)
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

      nlatx = 4001
      nlaty = 3001
      ntarg = nlatx*nlaty

      xmin = -4.0d0
      xmax = 4.0d0

      ymin = -3.0d0
      ymax = 3.0d0
      
      allocate(targs(2,ntarg))
      do ix=1,nlatx
        do iy=1,nlaty
          ipt = (ix-1)*nlaty + iy
          targs(1,ipt) = xmin + (xmax-xmin)*(ix-1.0d0)/(nlatx-1.0d0)
          targs(2,ipt) = ymax + (ymin-ymax)*(iy-1.0d0)/(nlaty-1.0d0)
        enddo
      enddo

      allocate(pottarg_plot(ntarg),ich_id_targ(ntarg))
      allocate(ts_targ(ntarg),potplot(ntarg))

      do i=1,ntarg
        ich_id_targ(i) = -1
        ts_targ(i) = 0
      enddo
      eps = 1.0d-7
      call lpcomp_helm_comb_dir_2d(nch,norders,ixys,iptype,npts,
     1  srccoefs,srcinfo,2,ntarg,targs,ich_id_targ,ts_targ,eps,
     2  zpars,solnwp,pottarg_plot)
      
      do i=1,ntarg
        potplot(i) = abs(pottarg_plot(i))
      enddo

      call pyimage4(35,nlaty,nlatx,potplot,1,srcinfo(1,1:npts),
     1  srcinfo(2,1:npts),npts,xmin,xmax,ymin,ymax,5)
      


      stop
      end
