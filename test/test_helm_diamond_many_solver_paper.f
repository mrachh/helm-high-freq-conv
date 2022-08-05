      implicit real *8 (a-h,o-z)
      complex *16, allocatable :: solnrefg(:),solnref2g(:),solncoefsg(:)
      
      real *8, allocatable :: err_est(:,:),err_dens(:,:),err_q(:,:)
      real *8, allocatable :: err_est_ref(:),err_est_ref2(:)
      real *8, allocatable :: err_dens_ref(:),err_q_ref(:)
      integer, allocatable :: niter_ref(:),niter_analytic_ref(:)
      integer, allocatable :: niter_ref2(:),niter_analytic_ref2(:)
      integer, allocatable :: niter(:,:),niter_analytic(:,:)
      real *8 dpars(2)
      real *8 shifts(2,100),rsc(100)
      complex *16 zk

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

      dpars(1) = 1.0d0
      dpars(2) = 1.0d0

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


      
c      k = 22
c      nover = 30
c      kref = k
c      kref2 = k-4
      k = 4
      nover = 10
      kref = k
      kref2 = k-2
      nkd = 2
      nppw = 3

      kerr = 32
      allocate(err_est_ref(nkd),err_est_ref2(nkd),err_dens_ref(nkd))
      allocate(niter_ref(nkd),niter_ref2(nkd))
      allocate(niter_analytic_ref(nkd),niter_analytic_ref2(nkd))
      allocate(err_q_ref(nkd))
c
c  reference solutions tested, 
c    and then write
c    nppw error estimator
c

      allocate(err_est(nppw,nkd),err_dens(nppw,nkd),err_q(nppw,nkd))
      allocate(niter(nppw,nkd),niter_analytic(nppw,nkd))
      open(unit=33,file='diamond_data/diamond_res.txt',access='append')
      do ik=1,nkd
        print *, "ik=",ik
        zk = (10.0d0*2**(ik-1))*sqrt(2.0d0) + 0.0d0
        zk = zk/4
        nchref0 = ceiling(0.4d0*4*abs(zk)/sqrt(2.0d0))
        ncherr0 = 2*nchref0
        nchref20 = nchref0-2

        nchref = ncomp*nchref0*4
        ncherr = ncomp*ncherr*4
        nchref2 = ncomp*nchref20*4 

        nptsrefg = nchref*kref
        nptsref2g = nchref2*kref2
        print *, "zk=",zk
        print *, "nchref=",nchref
        print *, "nchref2=",nchref2
        print *, "kref=",kref
        print *, "kref2=",kref2
        allocate(solnrefg(nptsrefg),solnref2g(nptsref2g))
        call get_diamond_sols(ncomp,shifts,rsc,zk,k,kref,nover,
     1     nchref0,nchref,nptsrefg,
     1     solnrefg,err_est_ref(ik),niter_analytic_ref(ik),
     1     niter_ref(ik))

        call get_diamond_sols(ncomp,shifts,rsc,zk,k,kref2,nover,
     1     nchref02,nchref2,nptsref2g,
     1     solnref2g,err_est_ref2(ik),niter_analytic_ref2(ik),
     1     niter_ref2(ik))
        
        print *, "nchref=",nchref

        call get_diamond_many_dens_error(ncomp,shifts,rsc,nchref0,
     1    nchref,kref,nptsrefg,
     1    solnrefg,nchref20,nchref2,kref2,nptsref2g,solnref2g,
     1    kerr,ncherr0,ncherr,
     1    err_dens_ref(ik),err_q_ref(ik))

        do ippw=1,nppw
          print *, ippw
          rexp = (ippw-1.0d0)/(nppw-1.0d0)*0.5d0
          dppw = 2*abs(zk)**rexp
          nch10 = ceiling(0.4d0*dppw*abs(zk)/sqrt(2.0d0))
          nch1 = ncomp*nch10*4
          k1 = 1
          npts1 = nch1*k1
          allocate(solncoefsg(npts1))

           call get_diamond_sols(ncomp,shifts,rsc,zk,k,k1,nover,
     1       nch10,nch1,npts1,
     1       solncoefsg,err_est(ippw,ik),niter_analytic(ippw,ik),
     1       niter(ippw,ik))
           print *, "here1"

        
           print *, "nchref=",nchref
           print *, "kref=",kref
           print *, "nch1=",nch1
           print *, "kerr=",kerr
           print *, "ncherr=",ncherr
           call get_diamond_many_dens_error(ncomp,shifts,rsc,nchref0,
     1       nchref,kref,nptsrefg,
     1       solnrefg,nch10,nch1,k1,npts1,solncoefsg,kerr,ncherr0,
     1       ncherr,err_dens(ippw,ik),err_q(ippw,ik))
           print *, "here2"
           drat = (nch1+0.0d0)/abs(zk)
       write(33,'(2x,e11.5,1x,i5,3(2x,e11.5),2(2x,i4),2(2x,e11.5))') 
     1     real(zk),
     1     nch1,drat,dppw,err_est(ippw,ik),
     1     niter_analytic(ippw,ik),niter(ippw,ik),err_dens(ippw,ik),
     2     err_q(ippw,ik)
          deallocate(solncoefsg)
        enddo
        deallocate(solnrefg,solnref2g)
      enddo
      call prin2('err_est=*',err_est_ref,nkd)
      call prin2('err_est2=*',err_est_ref2,nkd)
      call prin2('err_dens_ref=*',err_dens_ref,nkd)
      call prinf('niter_ref=*',niter_ref,nkd)
      call prinf('niter_ref2=*',niter_ref2,nkd)
      call prinf('niter_analytic_ref=*',niter_analytic_ref,nkd)
      call prinf('niter_analytic_ref2=*',niter_analytic_ref2,nkd)
      call prin2('err_q_ref=*',err_q_ref,nkd)
      close(33)
      


      return
      end


      subroutine get_diamond_sols(ncomp,shifts,rsc,zk,k,kg,
     1   nover,nch0,nch,nptsg,solncoefsg,
     1   err_est,niter1,niter2)
      
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcinfog(:,:),qwtsg(:),srccoefsg(:,:)
      real *8, allocatable :: srcover(:,:),wover(:),srccoefsover(:,:)
      real *8, allocatable :: ts1(:),ts1g(:),tsover(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tsg(:),umatg(:,:),vmatg(:,:),wtsg(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8 umo,vmo
      complex *16, allocatable :: sigmag(:),sigmacoefsg(:)
      complex *16, allocatable :: sigmapwg(:),sigmacoefspwg(:)
      complex *16 solncoefsg(nptsg)
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
      real *8 shifts(2,ncomp),rsc(ncomp)
      real *8, allocatable :: errs(:)
      complex *16 zfac
      real *8 xy_in(2),xy_out(2)
      data ima/(0.0d0,1.0d0)/
      external circ_geom

      

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

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

      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      allocate(tsg(kg),umatg(kg,kg),vmatg(kg,kg),wtsg(kg))
      call legeexps(itype,kg,tsg,umatg,vmatg,wtsg)
      nch = ncomp*nch0*4
      print *, "nch=",nch
      print *, "k=",k
      print *, "kg=",kg
      print *, "nover=",nover
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
      ndd_curv = 2
      ndz_curv = 0
      ndi_curv = 0
      
      
      call get_diamond_many(nch0,ncomp,rsc,shifts,nch,k,
     1  npts,adjs,srcinfo,srccoefs,qwts,norders,iptype,ixys)
      print *, "nch=",nch
      print *, "kg=",kg
      print *, "nptsg=",nptsg
      call get_diamond_many(nch0,ncomp,rsc,shifts,nch,kg,
     1  nptsg,adjs,srcinfog,srccoefsg,qwtsg,nordersg,iptype,ixysg)
      call prin2('a=*',a,1)
      call prin2('b=*',b,1)

      print *, "nch=",nch
      print *, "nover=",nover
      print *, "a=",a
      print *, "b=",b
      print *, "npts_over=",npts_over
      call get_diamond_many(nch0,ncomp,rsc,shifts,nch,nover,
     1  npts_over,adjs,srcover,
     1  srccoefsover,wover,novers,iptype,ixyso)



      h = 2*pi/(nch+0.0d0)
      print *, "nptsg=",nptsg
      allocate(sigmag(nptsg),sigmacoefsg(nptsg))
      allocate(sigmacoefspwg(nptsg),sigmapwg(nptsg))
      allocate(solng(nptsg))

      xyin(1) = 0.1d0*rsc(1) + shifts(1,1)
      xyin(2) = -0.33d0*rsc(1) + shifts(2,1)

      xyout(1) = 3.5d1
      xyout(2) = 3.3d1

      thet = pi/4

c
c  get density info
c
      do ich=1,nch
        do j=1,kg
          ipt = (ich-1)*kg + j
          call h2d_slp(srcinfog(1,ipt),2,xyin,ndd,dpars,1,zk,ndi,
     1       ipars,sigmag(ipt))
          zkuse = ceiling(real(zk)) - 1.0d0
          xx = srcinfog(1,ipt)
          yy = srcinfog(2,ipt)
          sigmapwg(ipt) = exp(ima*zk*(xx*cos(thet) + yy*sin(thet)))
          solncoefsg(ipt) = 0.0d0
        enddo
        istart = (ich-1)*kg + 1
        call dgemm('n','t',2,kg,kg,alpha,sigmag(istart),
     1     2,umatg,kg,beta,sigmacoefsg(istart),2)
        call dgemm('n','t',2,kg,kg,alpha,sigmapwg(istart),
     1     2,umatg,kg,beta,sigmacoefspwg(istart),2)
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
      eps = 0.51d-9
      eps_gmres = 0.51d-9

      niter = 0
      numit = max(ceiling(10*abs(zk)),200)
      allocate(errs(numit+1))
      ifinout = 1
      print *, "starting gmres"
      call helm_comb_dir_galerkin_solver2d(nch,kg,ixysg,nptsg,
     1  srcinfog,eps,zpars,nnzg,row_ptrg,col_indg,iquadg,
     2  nquadg,wnearcoefsg,novers(1),npts_over,ixyso,srcover,
     3  wover,numit,ifinout,sigmacoefsg,eps_gmres,niter,errs,rres,
     4  solncoefsg)
      niter1 = niter

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
      erra = abs(pottarg-pottargex)/abs(pottargex)
      err_est = erra 
      
      call helm_comb_dir_galerkin_solver2d(nch,kg,ixysg,nptsg,
     1  srcinfog,eps,zpars,nnzg,row_ptrg,col_indg,iquadg,
     2  nquadg,wnearcoefsg,novers(1),npts_over,ixyso,srcover,
     3  wover,numit,ifinout,sigmacoefspwg,eps_gmres,niter,errs,rres,
     4  solncoefsg)
      niter2 = niter
      print *, "erra=",erra

      
       

      return
      end
