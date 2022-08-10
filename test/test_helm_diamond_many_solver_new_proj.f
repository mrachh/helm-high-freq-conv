      implicit real *8 (a-h,o-z)

      nk = 2
      nppw = 20
      ifwrite = 1
      ifplot = 0
      do ik=1,nk
        do ippw=1,nppw
          call diamond_many(ik,ippw,ifwrite,iplot)
        enddo
      enddo



      stop
      end

      subroutine diamond_many(ik,ippw,ifwrite,ifplot)


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
      complex *16, allocatable :: soln(:)
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
      data ima/(0.0d0,1.0d0)/
      external circ_geom

      

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

      open(unit=133,file='diamond_data/diamond_res_newproj2.txt',
     1    access='append')
      ifwrite = 0
      ifplot = 1

      zk = 100.0d0 + 0.0d0*ima
      zk = (10.0d0*2**(ik-1))*sqrt(2.0d0) + 0.0d0
      if(ik.eq.0) zk = 55.0d0
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
      nover = 40
      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      allocate(tsg(kg),umatg(kg,kg),vmatg(kg,kg),wtsg(kg))
      call legeexps(itype,kg,tsg,umatg,vmatg,wtsg)

      allocate(tover(nover),wtsover(nover))
      itype = 1
      call legeexps(itype,nover,tover,umo,vmo,wtsover)

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

c      nch0 = ceiling(0.4*4*abs(zk)/sqrt(2.0d0))*4
c      nch = 4*nch0*ncomp

      if(ippw.le.10) drat = 10 + ippw*1.5d0
      if(ippw.gt.10) drat = 25 + (ippw-10)*2.5d0

      nch0 = ceiling(drat*abs(zk)/4/ncomp)
      nch = ncomp*nch0*4
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
      allocate(sigmag(nptsg),sigmacoefsg(nch),solncoefsg(nch))
      allocate(solng(nptsg))
      allocate(sigma(npts),sigmacoefs(nch),solncoefs(nch))
      allocate(sigmacoefs_full(npts),solncoefs_full(npts))
      allocate(sigmacoefsg_full(nptsg),solncoefsg_full(nptsg))
      allocate(soln(npts))

c
c  get density info
c
      thet = 0.087d0
      do ich=1,nch
        sigmacoefsg(ich) = 0
        ra = 0
        do j=1,kg
          ipt = (ich-1)*kg + j
          sigmag(ipt) = exp(ima*zk*(srcinfog(1,ipt)*cos(thet)+ 
     1         srcinfog(2,ipt)*sin(thet)))
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
          sigma(ipt) = exp(ima*zk*(srcinfo(1,ipt)*cos(thet)+ 
     1         srcinfo(2,ipt)*sin(thet)))
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
      eps = 0.51d-7

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
      
c      nchuse0 = nch0
c      nchuse = ncomp*4*nchuse0
c      ifwrite = 0
c      iunit = 39
c      kuse = 2*k
c
c      call get_diamond_many_dens_error(ncomp,shifts,rsc,
c     1  nch0,nch,k,npts,solncoefs_full,nch0,nch,1,nch,
c     2  solncoefs,kuse,nchuse0,nchuse,erra1,errq1,ifwrite,iunit)
c      print *, "erra=",erra1
c      print *, "errq=",errq1
       drat = (nch+0.0d0)/abs(zk)
       if(ifwrite.eq.1) then
       write(133,'(2x,e11.5,1x,i5,3(2x,e11.5),2(2x,i4),2(2x,e11.5))') 
     1     real(zk),
     1     nch,drat,dppw,0.0d0,0,0,erra,errq
      endif

      if(ifplot.eq.1) then
        ntlat = 601
        ntarg = ntlat*ntlat

        allocate(targs(2,ntarg))

        xmin = -10
        xmax = 10
        ymin = -10
        ymax = 10
        allocate(isin(ntarg),isin0(ntarg))
        do ix=1,ntlat
          do iy=1,ntlat
            ipt = (ix-1)*ntlat + iy
            targs(1,ipt) = xmin + (xmax-xmin)*(ix-1.0d0)/(ntlat-1.0d0)
            targs(2,ipt) = ymax + (ymin-ymax)*(iy-2.0d0)/(ntlat-1.0d0)
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

        allocate(potplot(ntarg))
        do i=1,ntarg
          potplot(i) = isin(i)
        enddo
        allocate(nptcomp(ncomp))
        do i=1,ncomp
          nptcomp(i) = 4*nch0*k
        enddo
        call pyimage4(33,ntlat,ntlat,potplot,ncomp,srcinfo(1,1:npts),
     1    srcinfo(2,1:npts),nptcomp,xmin,xmax,ymin,ymax,5)

        allocate(pottarg_plot(ntarg),pottargex_plot(ntarg))
        allocate(ich_id(ntarg),ts_pts(ntarg))
        do i=1,ntarg
          call h2d_slp(targs(1,i),2,xyin,ndd,dpars,1,zk,ndi,ipars,
     1       pottargex_plot(i))
          ich_id(i) = -1
          ts_pts(i) = 0
        enddo

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
      
        do i=1,ntarg
          potplot(i) = abs(pottarg_plot(i))
        enddo
        call pyimage4(34,ntlat,ntlat,potplot,ncomp,srcinfo(1,1:npts),
     1    srcinfo(2,1:npts),nptcomp,xmin,xmax,ymin,ymax,5)

      
        do i=1,ntarg
          potplot(i) = real(pottarg_plot(i))
        enddo
        call pyimage4(35,ntlat,ntlat,potplot,ncomp,srcinfo(1,1:npts),
     1    srcinfo(2,1:npts),nptcomp,xmin,xmax,ymin,ymax,5)

      
        errmax = -16.0d0
        do i=1,ntarg
          potplot(i) =
     1      log(abs(pottarg_plot(i)-pottargex_plot(i)))/log(10.0d0)
          if(isin(i).eq.-4) then
            if(potplot(i).ge.errmax) errmax = potplot(i)
          endif
        enddo
        print *, "max log10 ext error=",errmax
        call pyimage4(36,ntlat,ntlat,potplot,ncomp,srcinfo(1,1:npts),
     1    srcinfo(2,1:npts),nptcomp,xmin,xmax,ymin,ymax,5)

      endif
      close(133)
      


      return
      end
