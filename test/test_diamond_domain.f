      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),srccoefs(:,:),qwts(:)
      integer, allocatable :: norders(:),iptype(:),ixys(:)
      integer, allocatable :: adjs(:,:)
      real *8 xypt(2),shifts(2)

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4
      nch0 = 100
c      nch0 = 16
      nch = 4*nch0
c      k = 10
      k = 22
      npts = k*nch

      rfac = sqrt(2.0d0)
      rsc = 0.8d0*rfac*pi/2
      shifts(1) = rfac*pi
      shifts(2) = rfac*pi
      allocate(adjs(2,nch),srcinfo(8,npts),srccoefs(6,npts))
      allocate(qwts(npts),norders(nch),iptype(nch),ixys(nch+1))
      
      call get_diamond(nch0,nch,k,npts,adjs,srcinfo,srccoefs,qwts,
     1  norders,iptype,ixys)
      call prin2('srccoefs=*',srccoefs,6*k)
      do i=1,npts
        srcinfo(1:6,i) = srcinfo(1:6,i)*rsc
        qwts(i) = qwts(i)*rsc
        srccoefs(1:6,i) = srccoefs(1:6,i)*rsc
        srcinfo(1:2,i) = srcinfo(1:2,i) + shifts(1:2)
cc        srccoefs(1:2,i) = srccoefs(1:2,i) + shifts(1:2)
        write(37,*) srcinfo(1,i),srcinfo(2,i),srcinfo(7,i),srcinfo(8,i)
      enddo
      do i=1,nch
        istart = ixys(i)
        srccoefs(1:2,istart) = srccoefs(1:2,istart) + shifts(1:2)

      enddo
      call prin2('srcinfo=*',srcinfo,8)
      call prin2('srccoefs=*',srccoefs,6*k)
      
      xypt(1) = 0.6219884159152896D+01
      xypt(2) = 0.6220188067690530D+01

      write(38,*) xypt(1),xypt(2)
      call prinf('nch=*',nch,1)
      call prinf('norders=*',norders,nch)
      call prinf('ixys=*',ixys,nch+1)
      call prinf('iptype=*',iptype,nch)
      call prinf('npts=*',npts,1)
      call prinf('adjs=*',adjs,2*nch)
      ndt = 2
      nt = 1
      call findnearchunktarg_id_ts_brute(nch,norders,ixys,iptype,
     1  npts,srccoefs,srcinfo,ndt,nt,xypt,ich_id,ts,dist)
      call prin2('dist=*',dist,1)
      call prinf('ich_id=*',ich_id,1)
      call prin2('ts=*',ts,1)
      


      stop
      end
