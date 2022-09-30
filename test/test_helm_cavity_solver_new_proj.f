      implicit real *8 (a-h,o-z)

      nk = 0
      nppw = 20

      nks = 1
      nke = 1
      ippws = 1
      ippwe = 16
      ifwrite = 1
      ifplot = 1
      do ik=nks,nke
        ntlat = 901 
        do ippw=ippws,ippwe
          call cavity(ik,ippw,ifwrite,ifplot,ntlat)
        enddo
      enddo



      stop
      end

      subroutine cavity(ik,ippw,ifwrite,ifplot,ntlat)


      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcinfog(:,:),qwtsg(:),srccoefsg(:,:)
      real *8, allocatable :: srcover(:,:),wover(:),srccoefsover(:,:)
      real *8, allocatable :: ts1(:),ts1g(:),tsover(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tsg(:),umatg(:,:),vmatg(:,:),wtsg(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8 umo,vmo,pols(100)
      complex *16, allocatable :: sigmag(:),sigmacoefsg(:),solncoefsg(:)
      complex *16, allocatable :: solng(:)
      complex *16, allocatable :: sigma(:),sigmacoefs(:),solncoefs(:)
      complex *16, allocatable :: soln(:),soln_use(:)
      complex *16, allocatable :: sigmacoefs_full(:),solncoefs_full(:)
      complex *16, allocatable :: sigmacoefsg_full(:),solncoefsg_full(:)
      integer, allocatable :: row_ptrg(:),col_indg(:),iquadg(:)
      complex *16, allocatable :: wnearg(:),wnearcoefsg(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:),wnearcoefs(:)
      complex *16, allocatable :: vnval(:),vval(:),pnvval(:)
      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_pts(:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),ixysg(:)
      integer, allocatable :: nordersg(:)
      integer, allocatable :: novers(:),ixyso(:),adjs(:,:)
      real *8 shifts(2,100),rsc(100)
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: isin(:),isin0(:)
      complex *16, allocatable :: pottarg_plot(:),pottargex_plot(:)
      real *8, allocatable :: potplot(:)
      real *8, allocatable :: xscat(:),yscat(:)
      integer, allocatable :: nptcomp(:)
      complex *16 zk,zpars(3),ima,z1,z2,ztmp
      complex *16 pottarg,pottargex
      real *8 xyin(2),xyout(2)
      real *8 dpars(1)
      real *8, allocatable :: errs(:)
      complex *16 zfac
      real *8 xy_in(2),xy_out(2)

      complex *16, allocatable :: zpars_curv(:)

      data ima/(0.0d0,1.0d0)/
      external circ_geom,funcurv_zfft
      character *300 fname_sol,fplot1,fplot2,fname4,fname_res,dir_name1
      character *300 dir_name2

      

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4
      dir_name1 = 'cavity_data/'
      dir_name2='/mnt/home/mrachh/ceph/helm-high-freq-conv/' //
     1  'cavity-data/'
      write(fname_res,'(a,a)') trim(dir_name1),
     1  'cavity_res_sep30.txt'

      open(unit=133,file=trim(fname_res),access='append')

      zk = 10.0d0*ik + ima*0.0d0
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

      k = 16
      kg = 14
      nover = 24
      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      allocate(tsg(kg),umatg(kg,kg),vmatg(kg,kg),wtsg(kg))
      call legeexps(itype,kg,tsg,umatg,vmatg,wtsg)

      allocate(tover(nover),wtsover(nover))
      itype = 1
      call legeexps(itype,nover,tover,umo,vmo,wtsover)

      nch = ceiling(2*pi*abs(zk)*13.14d0*(2+0.5*ippw))

      write(fname_sol,'(a,a,i1,a,i5.5,a)') trim(dir_name2),
     1   'sol_ik',ik,'_nch',nch,'.bin'
      print *, trim(fname_sol)
      print *, trim(fname_res)
      open(unit=79,file=trim(fname_sol),form='unformatted')

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

      mm = 200
      ndz_curv = 2*mm

      ndd_curv = 0
      ndi_curv = 0

      allocate(zpars_curv(ndz_curv))

      aa = 0.2d0
      bb = pi/12

      call load_cavity_zpars(aa,bb,mm,zpars_curv,ndz_curv)

      a = 0
      b = 2.0d0*pi
      print *, "nch=",nch
      call get_funcurv_geom_uni(a,b,nch,k,npts,adjs,srcinfo,srccoefs,
     1  ts1,qwts,norders,iptype,ixys,funcurv_zfft,ndd_curv,dpars,
     2  ndz_curv,zpars_curv,ndi_curv,ipars)
      
      call get_funcurv_geom_uni(a,b,nch,kg,nptsg,adjs,srcinfog,
     1  srccoefsg,ts1g,qwtsg,nordersg,iptype,ixysg,
     2  funcurv_zfft,ndd_curv,dpars,ndz_curv,zpars_curv,ndi_curv,ipars)
      
      call get_funcurv_geom_uni(a,b,nch,nover,npts_over,adjs,srcover,
     1  srccoefsover,tsover,wover,novers,iptype,ixyso,
     2  funcurv_zfft,ndd_curv,dpars,ndz_curv,zpars_curv,ndi_curv,ipars)
      
      ra = sum(qwts)
      print *, "ra=",ra
      nppw = floor(nch*2*pi/ra/abs(zk))
      dppw = nch*2*pi/ra/abs(zk)
      print *, "nppw=",nppw
      print *, "dppw=",dppw
      



      print *, "nptsg=",nptsg
      allocate(sigmag(nptsg),sigmacoefsg(nch),solncoefsg(nch))
      allocate(solng(nptsg))
      allocate(sigma(npts),sigmacoefs(nch),solncoefs(nch))
      allocate(sigmacoefs_full(npts),solncoefs_full(npts))
      allocate(sigmacoefsg_full(nptsg),solncoefsg_full(nptsg))
      allocate(soln(npts),soln_use(npts))

c
c  get density info
      
      s0 = 1.0d0/sqrt(real(zk))
      thet = 0.45d0*pi
      ttuse = 0.3913220708069514D+01

      
      do ich=1,nch
        sigmacoefsg(ich) = 0
        ra = 0
        do j=1,kg
          ipt = (ich-1)*kg + j
          sigmag(ipt) = exp(-(ts1g(ipt)-ttuse)**2/s0**2)*
     1       exp(ima*zk*cos(thet)*ts1g(ipt))
          sigmacoefsg(ich) = sigmacoefsg(ich) + sigmag(ipt)*qwtsg(ipt)
          ra = ra + qwtsg(ipt)
        enddo
        sigmacoefsg(ich) = sigmacoefsg(ich)/ra
        istart = (ich-1)*kg+1
        call dgemm('n','t',2,kg,kg,alpha,sigmag(istart),2,umatg,
     1     kg,beta,sigmacoefsg_full(istart),2)
      enddo

      do ich=1,nch
        sigmacoefs(ich) = 0
        ra  = 0
        do j=1,k
          ipt = (ich-1)*k + j
          sigmag(ipt) = exp(-(ts1(ipt)-ttuse)**2/s0**2)*
     1       exp(ima*zk*cos(thet)*ts1(ipt))
          sigmacoefs(ich) = sigmacoefs(ich) + sigma(ipt)*qwts(ipt)
          ra = ra + qwts(ipt)
        enddo
        sigmacoefs(ich) = sigmacoefs(ich)/ra
        istart = (ich-1)*k+1
        call dgemm('n','t',2,k,k,alpha,sigma(istart),2,umat,
     1     k,beta,sigmacoefs_full(istart),2)
      enddo
      call prin2('sigmacoefsg=*',sigmacoefsg,24)
      call prin2('sigmacoefs=*',sigmacoefs,24)
c

      nnz = 3*k*nch
      nquad = nnz*k

      allocate(row_ptr(npts+1),col_ind(nnz),wnear(nquad))
      allocate(wnearcoefs(nquad))

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
      print *, "k=",k
      print *, "nnz=",nnz
      print *, "nquad=",nquad
      call cpu_time(t1)
C$       t1 = omp_get_wtime()     
      call get_helm_dir_trid_quad_corr(zk,nch,k,k,npts,npts,adjs,
     1   srcinfo,srcinfo,ndz,zpars,nnz,row_ptr,col_ind,
     2   nquad,wnear,wnearcoefs)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      call prin2('total quad gen time=*',t2-t1,1)
      call prinf('nnzg=*',nnz,1)

      allocate(iquad(nnz+1))
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,iquad)
      call cpu_time(t1)
C$       t1 = omp_get_wtime()     
      call get_helm_dir_trid_quad_corr(zk,nch,kg,kg,nptsg,nptsg,adjs,
     1   srcinfog,srcinfog,ndz,zpars,nnzg,row_ptrg,col_indg,
     2   nquadg,wnearg,wnearcoefsg)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      call prin2('total quad gen time=*',t2-t1,1)
      call prinf('nnzg=*',nnzg,1)

      allocate(iquadg(nnzg+1))
      call get_iquad_rsc2d(nch,ixysg,nptsg,nnzg,row_ptrg,col_indg,
     1   iquadg)
      iquadtype = 1
      eps = 0.51d-5

      niter = 0
      numit = max(ceiling(10*abs(zk)),200)
      allocate(errs(numit+1))
      ifinout = 1
c      call helm_comb_dir_galerkin_solver2d_new_proj(nch,kg,ixysg,nptsg,
c     1  srccoefsg,srcinfog,qwtsg,eps,zpars,nnzg,row_ptrg,col_indg,
c     2  iquadg,
c     2  nquadg,wnearg,novers(1),npts_over,ixyso,srcover,
c     3  wover,numit,ifinout,sigmacoefsg,eps,niter,errs,rres,
c     4  solncoefsg)

      call helm_comb_dir_galerkin_solver2d_new_proj(nch,k,ixys,npts,
     1  srccoefs,srcinfo,qwts,eps,zpars,nnz,row_ptr,col_ind,
     2  iquad,
     2  nquad,wnear,novers(1),npts_over,ixyso,srcover,
     3  wover,numit,ifinout,sigmacoefs,eps,niter,errs,rres,
     4  solncoefs)

c      call helm_comb_dir_galerkin_solver2d(nch,kg,ixysg,nptsg,
c     1  srcinfog,eps,zpars,nnzg,row_ptrg,col_indg,
c     2  iquadg,
c     2  nquadg,wnearcoefsg,novers(1),npts_over,ixyso,srcover,
c     3  wover,numit,ifinout,sigmacoefsg_full,eps,niter,errs,rres,
c     4  solncoefsg_full)

      call helm_comb_dir_galerkin_solver2d(nch,k,ixys,npts,
     1  srcinfo,eps,zpars,nnz,row_ptr,col_ind,
     2  iquad,
     2  nquad,wnearcoefs,novers(1),npts_over,ixyso,srcover,
     3  wover,numit,ifinout,sigmacoefs_full,eps,niter,errs,rres,
     4  solncoefs_full)

c      call prin2('solncoefsg=*',solncoefsg,24)
      call prin2('solncoefs=*',solncoefs,24)

c      call prin2('solncoefsg_full=*',solncoefsg_full,2*kg+2)
      call prin2('solncoefs_full=*',solncoefs_full,2*k+2)

      allocate(vval(npts_over),vnval(npts_over),pnvval(npts_over))
      do i=1,nch
        do j=1,nover
          ipt = (i-1)*nover + j
          vnval(ipt) = solncoefs(i)
          call legepols(tover(j),k-1,pols)
          vval(ipt) = 0
          pnvval(ipt) = 0
          do l=1,k
            lpt = (i-1)*k + l
            vval(ipt) = vval(ipt) + solncoefs_full(lpt)*pols(l)
          enddo
        enddo
      enddo

      do i=1,nch
        ra = 0
        ztmp = 0
        do j=1,nover
          ipt = (i-1)*nover + j
          ztmp = ztmp + vval(ipt)*wover(ipt)
          ra = ra + wover(ipt)
        enddo
        ztmp =ztmp/ra
        do j=1,nover
          ipt = (i-1)*nover+j
          pnvval(ipt) = ztmp
        enddo
      enddo
      call prin2('pnnval=*',pnvval,2)
      call prin2('vval=*',vval,2*nover)
      erra = 0
      errq = 0
      rnum = 0
      ra = 0
      do i=1,npts_over
        erra = erra + abs(vval(i)-vnval(i))**2*wover(i)
        rnum = rnum + abs(vval(i)-pnvval(i))**2*wover(i)
      enddo
      erra = sqrt(erra)
      rnum = sqrt(rnum)
      errq = rnum/erra

      print *, "erra=",erra
      print *, "errq=",errq
      write(79) zk
      write(79) nch
      write(79) k
      write(79) kg
      write(79) nover
      write(79) sigmacoefs
      write(79) sigmacoefs_full
      write(79) solncoefs
      write(79) solncoefs_full
      write(79) erra
      write(79) errq


      
       drat = (nch+0.0d0)/abs(zk)
       if(ifwrite.eq.1) then
       write(133,'(2x,e11.5,1x,i5,3(2x,e11.5),2(2x,i4),2(2x,e11.5))') 
     1     real(zk),
     1     nch,drat,dppw,0.0d0,0,0,erra,errq
      endif

      if(ifplot.eq.1) then
        ntarg = ntlat*ntlat

        allocate(targs(2,ntarg))

        xmin = -4
        xmax = 4
        ymin = -4
        ymax = 4
        allocate(isin(ntarg),isin0(ntarg))
        do ix=1,ntlat
          do iy=1,ntlat
            ipt = (ix-1)*ntlat + iy
            targs(1,ipt) = xmin + (xmax-xmin)*(ix-1.0d0)/(ntlat-1.0d0)
            targs(2,ipt) = ymax + (ymin-ymax)*(iy-1.0d0)/(ntlat-1.0d0)
            isin(ipt) = 0
          enddo
        enddo


        write(79) ntlat
        write(79) ntarg
        write(79) targs

        allocate(potplot(ntarg))
        allocate(pottarg_plot(ntarg),pottargex_plot(ntarg))
        allocate(ich_id(ntarg),ts_pts(ntarg))

        do i=1,nch
          do j=1,k
            ipt = (i-1)*k + j
            soln(ipt) = solncoefs(i)
          enddo
        enddo


      
        eps = 1.0d-7
        call lpcomp_helm_comb_dir_2d(nch,norders,ixys,iptype,npts,
     1    srccoefs,srcinfo,2,ntarg,targs,ich_id,ts_pts,eps,zpars,
     2    soln,pottarg_plot)
        write(79) pottarg_plot
      
        do i=1,ntarg
          potplot(i) = abs(pottarg_plot(i))
        enddo

      
      soln_use = 0
      do i=1,nch
        istart = (i-1)*k+1
        call dgemm('n','t',2,k,k,alpha,solncoefs_full(istart),
     1     2,vmat,k,beta,soln_use(istart),2)
      enddo
      
        eps = 1.0d-7
        call lpcomp_helm_comb_dir_2d(nch,norders,ixys,iptype,npts,
     1    srccoefs,srcinfo,2,ntarg,targs,ich_id,ts_pts,eps,zpars,
     2    soln_use,pottarg_plot)
        write(79) soln_use
        write(79) pottarg_plot
      
        do i=1,ntarg
          potplot(i) = abs(pottarg_plot(i))
        enddo

      endif
      close(133)
      close(79)
      


      return
      end



      subroutine get_nch0_nch(ik,ippw,nch0,nch)
      implicit real *8 (a-h,o-z)
      if(ik.eq.0) then
        if(ippw.eq.21) nch = 512
        if(ippw.eq.22) nch = 544
      endif

      if(ik.eq.2) then
        if(ippw.eq.21) nch = 608
        if(ippw.eq.22) nch = 640
        if(ippw.eq.23) nch = 656
      endif

      if(ik.eq.3) then
        if(ippw.eq.21) nch = 1232
        if(ippw.eq.22) nch = 1264
        if(ippw.eq.23) nch = 1184
        if(ippw.eq.24) nch = 1200
        if(ippw.eq.25) nch = 1216
        if(ippw.eq.26) nch = 1280
        if(ippw.eq.27) nch = 1296
        if(ippw.eq.28) nch = 1312
        if(ippw.eq.29) nch = 1328
        if(ippw.eq.30) nch = 2832*2
      endif

      if(ik.eq.4) then
        if(ippw.eq.21) nch = 2496-32
        if(ippw.eq.22) nch = 2496-16
        if(ippw.eq.23) nch = 2496+16
        if(ippw.eq.24) nch = 2496+32
        if(ippw.eq.25) nch = 2496+48
        if(ippw.eq.26) nch = 2496+64
        if(ippw.eq.27) nch = 2496+16*5
        if(ippw.eq.28) nch = 2496+16*6
        if(ippw.eq.29) nch = 2496+16*7
        if(ippw.eq.30) nch = 2496+16*8
        if(ippw.eq.31) nch = 7136
        if(ippw.eq.32) nch = 8000
        if(ippw.eq.33) nch = 6256
      endif



      nch0 = nch/16
      
      return
      end