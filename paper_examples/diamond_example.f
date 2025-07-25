      implicit real *8 (a-h,o-z)
      integer ippws(34)

      do i=1,20
        ippws(i) = i
      enddo
      
      ippws(31) = 96
      ippws(32) = 196
      ippws(33) = 296
      ippws(34) = 396

      nk = 0

      nks = 4
      nke = 4
      ifwrite = 1
      ifplot = 1
      ngk = 1
      do ik=nks,nke
        ntlat = 901 
        do ip=16,16
          ippw = ippws(ip)
          call diamond(ik,ippw,ngk,ifwrite,ifplot,ntlat)
        enddo
      enddo



      stop
      end

      subroutine diamond(ik,ippw,ngk,ifwrite,ifplot,ntlat)


      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      integer, allocatable :: norders(:),ixys(:),iptype(:)
      real *8, allocatable :: srcover(:,:),wover(:),srccoefsover(:,:)
      integer, allocatable :: novers(:),ixyso(:),adjs(:,:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8, allocatable :: ts1(:),tsover(:)
      real *8 umo,vmo,pols(100)
      complex *16, allocatable :: sigma(:), sigmacoefs(:,:)
      complex *16, allocatable :: soln(:), solncoefs(:,:)
      complex *16, allocatable :: soln_full(:)
      complex *16, allocatable :: sigmacoefs_full(:,:)
      complex *16, allocatable :: solncoefs_full(:,:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      integer, allocatable :: row_ptr_full(:),col_ind_full(:)
      integer, allocatable :: iquad_full(:)
      complex *16, allocatable :: wnear(:),wnear_full(:)
      complex *16, allocatable :: vnval(:),vval(:),pnvval(:)
      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_pts(:)
      real *8 shifts(2,100),rsc(100)
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: isin(:),isin0(:)
      complex *16, allocatable :: pottarg_plot(:),pottargex_plot(:)
      real *8, allocatable :: potplot(:)
      real *8, allocatable :: xscat(:),yscat(:)
      integer, allocatable :: nptcomp(:)
      complex *16 zk,zpars(3),ima,z1,z2,ztmp(2)
      complex *16 pottarg,pottargex
      real *8 xyin(2),xyout(2)
      real *8 dpars(1)
      real *8, allocatable :: errs(:),errs_full(:)
      complex *16 zfac
      real *8 xy_in(2),xy_out(2)
      integer ncomp

      complex *16, allocatable :: zpars_curv(:)

      data ima/(0.0d0,1.0d0)/
      external circ_geom,funcurv_zfft
      character *300 fname_sol,fplot1,fplot2,fname4,fname_res,dir_name1
      character *300 dir_name2

      

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

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


      dir_name1 = 'results/'
      dir_name2='/Users/mrachh/flatiron_data_ccmlin043/'//
     1     'ceph/helm-high-freq-conv/' //
     1     'diamond-data-jul2025-pw/'
      write(fname_res,'(a,a)') trim(dir_name1),
     1  'diamond_res_jul2025_pw.txt'

      open(unit=133,file=trim(fname_res),access='append')

      zk = (10.0d0*2**(ik-1))*sqrt(2.0d0) + 0.0d0

      ndz = 3
      zpars(1) = zk
      zpars(2) = -ima*zk
      zpars(3) = 1.0d0

      alpha = 1.0d0
      beta = 0.0d0

c
c  k is the boundary discretization order
c
c  ngk is the galerkin polynomial order
c  ngk_full is the galerkin polynomial order for reference solution
c
c  nover is the oversampled source order
c

      k = 16
      nover = 24
      ngk_full = 12
      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      allocate(tover(nover),wtsover(nover))
      itype = 1
      call legeexps(itype,nover,tover,umo,vmo,wtsover)

      if(ippw.le.20) drat = 25 + (ippw-1)*5.0d0
      print *, "drat=",drat
      
      nch0 = ceiling(drat*abs(zk)/4/ncomp/ngk)
      nch = ncomp*nch0*4
      print *, nch
      print *, real(zk)

      write(fname_sol,'(a,a,i2.2,a,i5.5,a,i1,a)') trim(dir_name2),
     1   'sol_ik',ik,'_nch',nch,'_ngk',ngk,'.bin'
      print *, trim(fname_sol)
      print *, trim(fname_res)
      open(unit=79,file=trim(fname_sol),form='unformatted')

      npts = nch*k
      npts_over = nch*nover
      allocate(srcinfo(8,npts),qwts(npts))
      allocate(srccoefs(6,npts),ts1(npts))
      allocate(norders(nch),iptype(nch),ixys(nch+1),adjs(2,nch))


      allocate(srccoefsover(6,npts_over))
      allocate(srcover(8,npts_over),wover(npts_over))
      allocate(novers(nch),ixyso(nch+1))
      allocate(tsover(npts_over))

      call get_diamond_many(nch0, ncomp, rsc, shifts, nch, k,
     1  npts, adjs, srcinfo, srccoefs, qwts, norders, iptype, 
     2  ixys)
      
      call get_diamond_many(nch0, ncomp, rsc, shifts, nch, nover,
     1  npts_over, adjs, srcover, srccoefsover, wover, novers, iptype, 
     2  ixyso)
      
      ra = sum(qwts)
      print *, "ra=",ra
      nppw = floor(ngk*nch*2*pi/ra/abs(zk))
      dppw = ngk*nch*2*pi/ra/abs(zk)
      print *, "nppw=",nppw
      print *, "dppw=",dppw
      


      allocate(sigma(npts),sigmacoefs(ngk,nch),solncoefs(ngk,nch))
      allocate(sigmacoefs_full(ngk_full,nch))
      allocate(solncoefs_full(ngk_full,nch))
      allocate(soln(npts),soln_full(npts))

c
c  get density info
      
      thet = 5*2*pi/360.0d0

      do ich=1,nch
        do j=1,ngk
          sigmacoefs(j,ich) = 0
          solncoefs(j,ich) = 0
        enddo
        do j=1,ngk_full
          sigmacoefs_full(j,ich) = 0
          solncoefs_full(j,ich) = 0
        enddo
        do j=1,k
          ipt = (ich-1)*k + j
          sigma(ipt) = exp(ima*zk*(srcinfo(1,ipt)*cos(thet)+ 
     1         srcinfo(2,ipt)*sin(thet)))
        enddo
        istart = (ich-1)*k+1
        call get_galerkin_proj(2,ngk,k,ts,sigma(istart),qwts(istart), 
     1      sigmacoefs(1,ich))
        call get_galerkin_proj(2,ngk_full,k,ts,sigma(istart),
     1      qwts(istart),sigmacoefs_full(1,ich))
      enddo

      call prin2('sigmacoefs=*',sigmacoefs(1,1),4*ngk)
      call prin2('sigmacoefs_full=*',sigmacoefs_full(1,1),4*ngk_full)
c

      nb = ngk*nch
      nnz = 3*ngk*nch
      nquad = nnz*ngk

      nb_full = ngk_full*nch
      nnz_full = 3*ngk_full*nch
      nquad_full = nnz_full*ngk_full

      allocate(row_ptr(nb+1),col_ind(nnz),wnear(nquad))
      allocate(row_ptr_full(nb_full+1),col_ind_full(nnz_full))
      allocate(wnear_full(nquad_full))


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
      call get_helm_dir_trid_quad_corr(zk,nch,k,ngk,npts,nb,adjs,
     1   srcinfo,srccoefs,ndz,zpars,nnz,row_ptr,col_ind,
     2   nquad,wnear)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      call prin2('total quad gen time=*',t2-t1,1)
      call prinf('nnzg=*',nnz,1)

      call cpu_time(t1)
C$       t1 = omp_get_wtime()     
      call get_helm_dir_trid_quad_corr(zk,nch,k,ngk_full,npts,
     1   nb_full,adjs,srcinfo,srccoefs,ndz,zpars,nnz_full,
     2   row_ptr_full,col_ind_full,nquad_full,wnear_full)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      call prin2('total quad gen time=*',t2-t1,1)

      allocate(iquad(nnz+1))
      allocate(iquad_full(nnz_full+1))


      do i=1,nnz+1
        iquad(i) = (i-1)*ngk+1
      enddo

      do i=1,nnz_full+1
        iquad_full(i) = (i-1)*ngk_full + 1
      enddo


      eps = 0.51d-7
      niter = 0
      numit = max(ceiling(10*abs(zk)),200)

      niter_full = 0
      allocate(errs(numit+1))
      allocate(errs_full(numit+1))
      ifinout = 1

      call helm_comb_dir_galerkin_solver2d(nch,k,ngk,ixys,npts,nb,
     1  srcinfo,adjs,qwts,eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,
     2  wnear,novers(1),npts_over,ixyso,srcover,wover,numit,ifinout,
     3  sigmacoefs,eps,niter,errs,rres,solncoefs)

      call helm_comb_dir_galerkin_solver2d(nch,k,ngk_full,ixys,npts,
     1  nb_full,srcinfo,adjs,qwts,eps,zpars,nnz_full,row_ptr_full,
     2  col_ind_full,iquad_full,nquad_full,wnear_full,novers(1),
     3  npts_over,ixyso,srcover,wover,numit,ifinout,sigmacoefs_full,
     4  eps,niter_full,errs,rres,solncoefs_full)

      call prin2('solncoefs=*',solncoefs,24)
      call prin2('solncoefs_full=*',solncoefs_full,2*k+2)

      allocate(vval(npts_over),vnval(npts_over),pnvval(npts_over))
      do i=1,nch
        do j=1,nover
          ipt = (i-1)*nover + j
          call legepols(tover(j),k-1,pols)
          vnval(ipt) = 0
          vval(ipt) = 0
          pnvval(ipt) = 0
          do l=1,ngk
            vnval(ipt) = vnval(ipt) + pols(l)*solncoefs(l,i)
            pnvval(ipt) = pnvval(ipt) + pols(l)*solncoefs_full(l,i)
          enddo
          do l=1,ngk_full
            vval(ipt) = vval(ipt) + pols(l)*solncoefs_full(l,i)
          enddo
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
      write(79) ncomp
      write(79) rsc(1:ncomp)
      write(79) shifts(1:2,1:ncomp)
      write(79) zk
      write(79) nch0
      write(79) nch
      write(79) thet
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
       write(133,'(2x,e11.5,1x,i5,4(2x,e11.5),2x,i1)') 
     1     real(zk),nch,drat,dppw,erra,errq,ngk
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

        do icomp=1,ncomp
           istart = 4*k*nch0*(icomp-1)+1
           nchuse = 4*nch0
           npts0 = nchuse*k
           call chunk_interior(nchuse,norders,ixys,iptype,npts0,
     1        srcinfo(1,istart),srccoefs(1,istart),targs,2,ntarg,isin0)
           do i=1,ntarg
             isin(i) = isin(i) + isin0(i)
           enddo
        enddo

        write(79) ntlat
        write(79) ntarg
        write(79) targs
        write(79) isin

        allocate(potplot(ntarg))
        allocate(pottarg_plot(ntarg),pottargex_plot(ntarg))
        allocate(ich_id(ntarg),ts_pts(ntarg))

        do i=1,nch
          do j=1,k
            ipt = (i-1)*k + j
            call legepols(ts(j),k-1,pols)
            soln(ipt) = 0
            soln_full(ipt) = 0
            do l=1,ngk
              soln(ipt) = soln(ipt) + solncoefs(l,i)*pols(l)
            enddo
            do l=1,k
              soln_full(ipt) = soln_full(ipt) + 
     1           solncoefs_full(l,i)*pols(l)
            enddo
          enddo
        enddo

        do i=1,ntarg
          ich_id(i) = -1
          ts_pts(i) = 0
        enddo


      
        eps = 1.0d-7
        call lpcomp_helm_comb_dir_2d(nch,norders,ixys,iptype,npts,
     1    srccoefs,srcinfo,2,ntarg,targs,ich_id,ts_pts,eps,zpars,
     2    soln,pottarg_plot)
        write(79) pottarg_plot
      
        do i=1,ntarg
          potplot(i) = abs(pottarg_plot(i))
        enddo

      
      
        eps = 1.0d-7
        call lpcomp_helm_comb_dir_2d(nch,norders,ixys,iptype,npts,
     1    srccoefs,srcinfo,2,ntarg,targs,ich_id,ts_pts,eps,zpars,
     2    soln_full,pottarg_plot)
        write(79) soln_full
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
