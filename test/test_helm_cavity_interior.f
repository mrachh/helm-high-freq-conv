      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcinfog(:,:),qwtsg(:),srccoefsg(:,:)
      real *8, allocatable :: srcover(:,:),wover(:),srccoefsover(:,:)
      real *8, allocatable :: ts1(:),ts1g(:),tsover(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tsg(:),umatg(:,:),vmatg(:,:),wtsg(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8 umo,vmo,pols(100)
      complex *16, allocatable :: sigmag(:),sigmacoefsg(:,:)
      complex *16, allocatable :: solncoefsg(:,:)
      complex *16, allocatable :: solng(:)
      complex *16, allocatable :: sigma(:),sigmacoefs(:,:)
      complex *16, allocatable :: solncoefs(:,:)
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
      complex *16 zk,zpars(3),ima,z1,z2,ztmp(2)
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

      nch = 2024

      write(fname_sol,'(a,a,i2.2,a,i5.5,a)') trim(dir_name2),
     1   'cavity_interior.bin'
      print *, trim(fname_sol)
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

      ntlat = 901

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

      call chunk_interior(nch,norders,ixys,iptype,npts,srcinfo,
     1   srccoefs,targs,2,ntarg,isin)


      write(79) ntlat
      write(79) ntarg
      write(79) targs
      write(79) isin
      write(79) qwts

      close(79)
      


      
      stop
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
