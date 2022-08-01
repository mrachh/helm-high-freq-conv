      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcinfog(:,:),qwtsg(:),srccoefsg(:,:)
      real *8, allocatable :: srcover(:,:),wover(:),srccoefsover(:,:)
      real *8, allocatable :: ts1(:),ts1g(:),tsover(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tsg(:),umatg(:,:),vmatg(:,:),wtsg(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8 umo,vmo
      complex *16, allocatable :: sigma(:),potex(:),pot(:)
      complex *16, allocatable :: sigmag(:),potexg(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      integer, allocatable :: row_ptrg(:),col_indg(:),iquadg(:)
      complex *16, allocatable :: wnear(:),wnearcoefs(:)
      complex *16, allocatable :: wnearg(:),wnearcoefsg(:)
      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_pts(:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),ixysg(:)
      integer, allocatable :: nordersg(:)
      integer, allocatable :: novers(:),ixyso(:)
      complex *16 zk,zpars(3),ima,z1,z2
      complex *16, allocatable :: sigmacoefsg(:),potcoefsg(:)
      complex *16, allocatable :: sigmacoefs(:),potcoefs(:)
      complex *16, allocatable :: potexcoefsg(:),potexcoefs(:)
      complex *16 fjvals(0:100),fhvals(0:100),fjders(0:100)
      real *8 dpars(1)
      complex *16 fhders(0:100),zfac
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

      imode = 7

      nterms = imode + 5
      rscale = 1.0d0
      ifder = 1
      call jbessel2d(nterms,zk,rscale,fjvals,ifder,fjders)
      call h2dall(nterms,zk,rscale,fhvals,ifder,fhders)


      zfac = pi/2*ima*zk*fjders(imode)*fhvals(imode)*zpars(3)
      zfac = zfac + pi/2*ima*fjvals(imode)*fhvals(imode)*zpars(2)

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


      nch = 5
      npts = nch*k
      npts_over = nch*nover
      nptsg = nch*kg
      allocate(srcinfo(8,npts),qwts(npts))
      allocate(srccoefs(6,npts),ts1(npts))
      allocate(norders(nch),iptype(nch),ixys(nch+1))


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
      ndd_curv = 1
      ndz_curv = 0
      ndi_curv = 0
      call get_funcurv_geom_uni(a,b,nch,k,npts,srcinfo,
     1  srccoefs,ts1,qwts,norders,iptype,ixys,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      print *, "nch=",nch
      print *, "kg=",kg
      print *, "nptsg=",nptsg
      call get_funcurv_geom_uni(a,b,nch,kg,nptsg,srcinfog,
     1  srccoefsg,ts1g,qwtsg,nordersg,iptype,ixysg,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)
      call prin2('a=*',a,1)
      call prin2('b=*',b,1)

      print *, "nch=",nch
      print *, "nover=",nover
      print *, "a=",a
      print *, "b=",b
      print *, "npts_over=",npts_over
      call get_funcurv_geom_uni(a,b,nch,nover,npts_over,srcover,
     1  srccoefsover,tsover,wover,novers,iptype,ixyso,circ_geom,
     2  ndd_curv,dpars,ndz_curv,zpars,ndi_curv,ipars)



      h = 2*pi/(nch+0.0d0)
      print *, "nptsg=",nptsg
      allocate(sigma(npts),sigmacoefsg(nptsg),pot(npts))
      allocate(sigmag(nptsg))
      allocate(potcoefsg(nptsg))
      allocate(sigmacoefs(npts),potcoefs(npts))
      allocate(potexcoefs(npts))

      allocate(potex(npts),potexg(nptsg),potexcoefsg(nptsg))

c
c  get density info
c
      do ich=1,nch

        do j=1,k
          ipt = (ich-1)*k + j
          tuse = ts1(ipt)
          sigma(ipt) = exp(ima*(imode+0.0d0)*tuse)
          potex(ipt) = zfac*sigma(ipt)
        enddo
        istart = (ich-1)*k+1

        call dgemm('n','t',2,k,k,alpha,sigma(istart),
     1     2,umat,k,beta,sigmacoefs(istart),2)
        call dgemm('n','t',2,k,k,alpha,potex(istart),
     1     2,umat,k,beta,potexcoefs(istart),2)
        do j=1,kg
          ipt = (ich-1)*kg + j
          tuse = ts1g(ipt) 
          sigmag(ipt) = exp(ima*(imode+0.0d0)*tuse)
          potexg(ipt) = zfac*sigmag(ipt)
        enddo
        istart = (ich-1)*kg + 1
        call dgemm('n','t',2,kg,kg,alpha,sigmag(istart),
     1     2,umatg,kg,beta,sigmacoefsg(istart),2)
        call dgemm('n','t',2,kg,kg,alpha,potexg(istart),
     1     2,umatg,kg,beta,potexcoefsg(istart),2)

      enddo
      ra = sum(wover)
      ra2 = sum(qwts)
      call prin2('Error in perimeter of circle=*',abs(ra-2*pi),1)
      call prin2('Error in perimeter of circle=*',abs(ra2-2*pi),1)
      call prin2('sigma=*',sigma,32)
      call prin2('sigmacoefs=*',sigmacoefs,32)


c
c  test near quad at second patch from first patch
c
c
      nnz = 3*k*nch

      nnzg = 3*kg*nch
      nquadg = nnz*kg
      allocate(row_ptr(npts+1),col_ind(nnz),wnear(nquad))
      allocate(wnearcoefs(nquad))

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
      call get_helm_dir_trid_quad_corr(zk,nch,k,kg,npts,nptsg,srcinfo,
     1   srcinfog,ndz,zpars,nnzg,row_ptrg,col_indg,nquadg,wnearg,
     2   wnearcoefsg)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      call prin2('total quad gen time=*',t2-t1,1)
      call prinf('nnzg=*',nnzg,1)
      call prinf('row_ptrg=*',row_ptrg,nptsg+1)


      allocate(iquadg(nnzg+1))
      call get_iquad_rsc2d(nch,ixysg,nptsg,nnzg,row_ptrg,col_indg,
     1   iquadg)
      iquadtype = 1
      eps = 0.51d-14

      potcoefsg = 0
      print *, "novers=",novers(1)
      call lpcomp_galerkin_helm2d(nch,kg,ixysg,nptsg,
     1  srcinfog,eps,ndz,zpars,nnzg,row_ptrg,col_indg,iquadg,
     2  nquadg,wnearcoefsg,sigmacoefsg,novers(1),npts_over,ixyso,
     3  srcover,wover,potcoefsg)
      erra = 0
      ra = 0
      do i=1,nptsg
        potcoefsg(i) = potcoefsg(i) + sigmacoefsg(i)/2*zpars(3)
        erra = erra + abs(potcoefsg(i)-potexcoefsg(i))**2
        if(i.le.5) print *, real(potcoefsg(i)),real(potexcoefsg(i)),
     1      real(potcoefsg(i))/real(potexcoefsg(i))
        ra = ra + abs(potexcoefsg(i))**2
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in pot galerkin=*',erra,1)

       

      return
      end
