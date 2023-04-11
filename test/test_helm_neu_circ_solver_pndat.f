      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcinfog(:,:),qwtsg(:),srccoefsg(:,:)
      real *8, allocatable :: srcinfo2(:,:),qwts2(:),srccoefs2(:,:)
      real *8, allocatable :: srcover(:,:),wover(:),srccoefsover(:,:)
      real *8, allocatable :: ts1(:),ts1g(:),tsover(:),ts12(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tsg(:),umatg(:,:),vmatg(:,:),wtsg(:)
      real *8, allocatable :: ts2(:),umat2(:,:),vmat2(:,:),wts2(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8 umo,vmo
      complex *16, allocatable :: sigmag(:),sigmacoefsg(:),solncoefsg(:)
      complex *16, allocatable :: solnikcoefsg(:)
      complex *16, allocatable :: solng(:),solnikg(:)
      complex *16, allocatable :: sigma(:),sigma2(:)
      complex *16, allocatable :: sigmacoefs(:),sigmacoefs2(:)
      complex *16, allocatable :: soln(:),soln2(:)
      complex *16, allocatable :: solnik(:),solnik2(:)
      complex *16, allocatable :: solncoefs(:),solncoefs2(:)
      complex *16, allocatable :: solnikcoefs(:),solnikcoefs2(:)
      integer, allocatable :: row_ptrg(:),col_indg(:),iquadg(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      integer, allocatable :: row_ptr2(:),col_ind2(:),iquad2(:)
      complex *16, allocatable :: wnearg(:,:),wnearcoefsg(:,:)
      complex *16, allocatable :: wnear(:,:),wnear2(:,:)
      complex *16, allocatable :: wnearcoefs(:,:),wnearcoefs2(:,:)
      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_pts(:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),ixysg(:)
      integer, allocatable :: nordersg(:)
      integer, allocatable :: norders2(:),ixys2(:),iptype2(:),adjs2(:,:)
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

      ik = 7
      ippw = 10
      nppw = 20

      zk = 10.0d0*2**(ik-1) + 0.0d0
cc      zk = 1.0d0
      ndz = 3
      zpars(1) = zk
      zpars(2) = zk

      alpha = 1.0d0
      beta = 0.0d0
c
c  k is the boundary discretization order
c  k2 is different discretization order for verification
c
c  kg is the galerkin polynomial order for projection
c  
c
c  nover is the oversampled source order
c

      k = 22
      k2 = 20
      kg = 1
      nover = 30
      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      allocate(ts2(k2),umat2(k2,k2),vmat2(k2,k2),wts2(k2))
      call legeexps(itype,k2,ts2,umat2,vmat2,wts2)

      allocate(tsg(kg),umatg(kg,kg),vmatg(kg,kg),wtsg(kg))
      call legeexps(itype,kg,tsg,umatg,vmatg,wtsg)

 
      rexp = (ippw-1.0d0)/(nppw-1.0d0)/3.0d0
      if(ik.gt.5) rexp = (ippw-1.0d0)/(nppw-1.0d0)/5.0d0
      dppw = 2*abs(zk)**rexp
      nch = ceiling(abs(zk)*dppw)
      npts = nch*k
      npts2 = nch*k2
      npts_over = nch*nover
      nptsg = nch*kg
      allocate(srcinfo(8,npts),qwts(npts))
      allocate(srccoefs(6,npts),ts1(npts))
      allocate(norders(nch),iptype(nch),ixys(nch+1),adjs(2,nch))

      allocate(srcinfo2(8,npts2),qwts2(npts2))
      allocate(srccoefs2(6,npts2),ts12(npts2))
      allocate(norders2(nch),iptype2(nch),ixys2(nch+1),adjs2(2,nch))


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
      dpars(2) = 1.0d0
      ndd_curv = 2
      ndz_curv = 0
      ndi_curv = 0
      call get_funcurv_geom_uni(a,b,nch,k,npts,adjs,srcinfo,
     1  srccoefs,ts1,qwts,norders,iptype,ixys,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)

      call get_funcurv_geom_uni(a,b,nch,k2,npts2,adjs2,srcinfo2,
     1  srccoefs2,ts12,qwts2,norders2,iptype2,ixys2,circ_geom,
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
      allocate(solng(nptsg),solnikcoefsg(nptsg),solnikg(nptsg))

      allocate(sigma(npts),sigmacoefs(npts))
      allocate(soln(npts),solncoefs(npts))
      allocate(solnik(npts),solnikcoefs(npts))

      allocate(sigma2(npts2),sigmacoefs2(npts2))
      allocate(soln2(npts2),solncoefs2(npts2))
      allocate(solnik2(npts2),solnikcoefs2(npts2))

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
          zkuse = ceiling(real(zk)) - 1.0d0
          sigmag(ipt) = exp(ima*zkuse*ts1g(ipt))
          solncoefsg(ipt) = 0.0d0
          solnikcoefsg(ipt) = 0.0d0
        enddo
        istart = (ich-1)*kg + 1
        call dgemm('n','t',2,kg,kg,alpha,sigmag(istart),
     1     2,umatg,kg,beta,sigmacoefsg(istart),2)
      enddo
      call prin2('sigmacoefsg=*',sigmacoefsg,4)

      do ich=1,nch
        do j=1,k
          ipt = (ich-1)*k + j
          zkuse = ceiling(real(zk)) - 1.0d0
          sigma(ipt) = exp(ima*zkuse*ts1(ipt))
          solncoefs(ipt) = 0.0d0
          solnikcoefs(ipt) = 0.0d0
        enddo
        istart = (ich-1)*k + 1
        call dgemm('n','t',2,k,k,alpha,sigma(istart),
     1     2,umat,k,beta,sigmacoefs(istart),2)
        sigmacoefs(istart) = sigmacoefs(istart)-sigmacoefsg(ich)
      enddo
      call prin2('sigmacoefs=*',sigmacoefs,4*k)


      do ich=1,nch
        do j=1,k2
          ipt = (ich-1)*k2 + j
          zkuse = ceiling(real(zk)) - 1.0d0
          sigma2(ipt) = exp(ima*zkuse*ts12(ipt))
          solncoefs2(ipt) = 0.0d0
          solnikcoefs2(ipt) = 0.0d0
        enddo
        istart = (ich-1)*k2 + 1
        call dgemm('n','t',2,k2,k2,alpha,sigma2(istart),
     1     2,umat2,k2,beta,sigmacoefs2(istart),2)
        sigmacoefs2(istart) = sigmacoefs2(istart)-sigmacoefsg(ich)
      enddo
      call prin2('sigmacoefs2=*',sigmacoefs2,4*k2)





      ra = sum(wover)
      ra2 = sum(qwts)
      ra3 = sum(qwtsg)


c
c  test near quad at second patch from first patch
c
c

      nnz = 3*k*nch
      nquad = nnz*k

      allocate(row_ptr(npts+1),col_ind(nnz),wnear(nquad,4))
      allocate(wnearcoefs(nquad,4))

      nnz2 = 3*k2*nch
      nquad2 = nnz2*k2

      allocate(row_ptr2(npts2+1),col_ind2(nnz2),wnear2(nquad2,4))
      allocate(wnearcoefs2(nquad2,4))

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
      call get_helm_neu_trid_quad_corr(zk,nch,k,k,npts,npts,adjs,
     1   srcinfo,srcinfo,ndz,zpars,nnz,row_ptr,col_ind,nquad,
     2   wnear,wnearcoefs)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()    
      call prin2('total quad gen time=*',t2-t1,1)
      call cpu_time(t1)
C$       t1 = omp_get_wtime()     
      call get_helm_neu_trid_quad_corr(zk,nch,k,k2,npts,npts2,adjs,
     1   srcinfo,srcinfo2,ndz,zpars,nnz2,row_ptr2,col_ind2,nquad2,
     2   wnear2,wnearcoefs2)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      call prin2('total quad gen time=*',t2-t1,1)
      call prinf('nnz2=*',nnz2,1)
      allocate(iquad(nnz+1))
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,
     1   iquad)

      allocate(iquad2(nnz2+1))
      call get_iquad_rsc2d(nch,ixys2,npts2,nnz2,row_ptr2,col_ind2,
     1   iquad2)


      iquadtype = 1
      eps = 0.51d-8

      niter = 0
      numit = max(ceiling(10*abs(zk)),200)
      allocate(errs(numit+1))
      ifinout = 1
      call helm_comb_neu_galerkin_solver2d(nch,k,ixys,npts,
     1  srcinfo,eps,zpars,nnz,row_ptr,col_ind,iquad,
     2  nquad,wnearcoefs,novers(1),npts_over,ixyso,srcover,
     3  wover,numit,ifinout,sigmacoefs,eps,niter,errs,rres,
     4  solncoefs,solnikcoefs)

      iquadtype = 1
      eps = 0.51d-8

      niter = 0
      ifinout = 1
      call helm_comb_neu_galerkin_solver2d(nch,k2,ixys2,npts2,
     1  srcinfo2,eps,zpars,nnz2,row_ptr2,col_ind2,iquad2,
     2  nquad2,wnearcoefs2,novers(1),npts_over,ixyso,srcover,
     3  wover,numit,ifinout,sigmacoefs2,eps,niter,errs,rres,
     4  solncoefs2,solnikcoefs2)
      
      call prin2('solncoefs=*',solncoefs,2*k)
      call prin2('solncoefs2=*',solncoefs2,2*k2)

      call prin2('solnikcoefs=*',solnikcoefs,2*k)
      call prin2('solnikcoefs2=*',solnikcoefs2,2*k2)

      solncoefsg = 0

      nchover = 2*nch
      call get_circ_dens_error2(dpars,nch,k,npts,solncoefs,
     1   nch,kg,nptsg,solncoefsg,nover,nchover,err_dens,err_q)
      call get_circ_dens_error2(dpars,nch,k2,npts2,solncoefs2,
     1   nch,kg,nptsg,solncoefsg,nover,nchover,err_dens2,err_q2)
      print *, "err_dens=",err_dens
      print *, "err_dens2=",err_dens2
      err_est = abs(err_dens-err_dens2)
      print *, "err_q=",err_q
      print *, "err_q2=",err_q2
      drat = (nch+0.0d0)/real(zk)
      open(unit=33,file='circ_neu_data/circ_res_pn.txt',access='append')
       write(33,'(2x,e11.5,1x,i5,3(2x,e11.5),2(2x,i4),2(2x,e11.5))') 
     1     real(zk),nch,drat,dppw,err_est,0,0,err_dens,err_q
      

       

      return
      end
