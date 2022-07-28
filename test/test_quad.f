      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: tover(:),wtsover(:)
      real *8 umo,vmo
      complex *16, allocatable :: sigma(:),potex(:),pot(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:),wnearcoefs(:)
      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_pts(:)
      integer, allocatable :: norders(:),ixys(:),iptype(:)
      integer, allocatable :: novers(:),ixyso(:)
      complex *16 zk,zpars(3),ima,z1,z2
      complex *16, allocatable :: sigmacoefs(:),potcoefs(:)
      complex *16 fjvals(0:100),fhvals(0:100),fjders(0:100)
      complex *16 fhders(0:100),zfac
      real *8 xy_in(2),xy_out(2)
      data ima/(0.0d0,1.0d0)/

      

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

      zk = 1.1d0 + 0.0d0*ima
      zpars(1) = zk
      zpars(2) = -ima*zk
      zpars(3) = 1.0d0

      imode = 7

      nterms = imode + 5
      rscale = 1.0d0
      ifder = 1
      call jbessel2d(nterms,zk,rscale,fjvals,ifder,fjders)
      call h2dall(nterms,zk,rscale,fhvals,ifder,fhders)

      call prin2('fjvals=*',fjvals,2*(imode+1))
      call prin2('fjders=*',fjders,2*(imode+1))
      call prin2('fhvals=*',fhvals,2*(imode+1))
      call prin2('fhders=*',fhders,2*(imode+1))

      zfac = pi/2*ima*zk*fjders(imode)*fhvals(imode)*zpars(3)
      zfac = zfac + pi/2*ima*fjvals(imode)*fhvals(imode)*zpars(2)

      alpha = 1.0d0
      beta = 0.0d0

      k = 16
      nover = 24
      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      allocate(tover(nover),wtsover(nover))
      itype = 1
      call legeexps(itype,nover,tover,umo,vmo,wtsover)


      nch = 14
      allocate(srcinfo(8,k*nch),qwts(k*nch))
      allocate(srccoefs(6,k*nch))
      allocate(srcover(8,nover*nch),wover(nover*nch))
      allocate(norders(nch),iptype(nch),ixys(nch+1))
      allocate(novers(nch),ixyso(nch+1))
      h = 2*pi/(nch+0.0d0)
      npts = nch*k
      npts_over = nch*nover
      allocate(sigma(npts),sigmacoefs(npts),pot(npts),potcoefs(npts))
      allocate(potex(npts))
      do ich=1,nch
        tstart = (ich-1.0d0)*h
        tend = (ich+0.0d0)*h
        do j=1,k
          ipt = (ich-1)*k + j
          tuse = tstart + (tend-tstart)*(ts(j)+1)/2
          srcinfo(1,ipt) = cos(tuse)
          srcinfo(2,ipt) = sin(tuse)
          srcinfo(3,ipt) = -sin(tuse)*h/2
          srcinfo(4,ipt) = cos(tuse)*h/2
          srcinfo(5,ipt) = -cos(tuse)*h*h/2/2
          srcinfo(6,ipt) = -sin(tuse)*h*h/2/2
          srcinfo(7,ipt) = cos(tuse)
          srcinfo(8,ipt) = sin(tuse)
          sigma(ipt) = exp(ima*(imode+0.0d0)*tuse)
          potex(ipt) = zfac*sigma(ipt)
          qwts(ipt) = h/2*wts(j)
        enddo
        istart = (ich-1)*k+1

        call dgemm('n','t',6,k,k,alpha,srcinfo(1,istart),
     1    8,umat,k,beta,srccoefs(1,istart),6)
        call prin2('tover*',tover,nover)
        call prin2('wts=*',wtsover,nover)
        do j=1,nover
          ipt = (ich-1)*nover + j
          tuse = tstart + (tend-tstart)*(tover(j)+1)/2
          srcover(1,ipt) = cos(tuse)
          srcover(2,ipt) = sin(tuse)
          srcover(3,ipt) = -sin(tuse)*h/2
          srcover(4,ipt) = cos(tuse)*h/2
          srcover(5,ipt) = -cos(tuse)*h*h/2/2
          srcover(6,ipt) = -sin(tuse)*h*h/2/2
          srcover(7,ipt) = cos(tuse)
          srcover(8,ipt) = sin(tuse)
          wover(ipt) = h/2*wtsover(j)
        enddo
        norders(ich) = k
        iptype(ich) = 1
        ixys(ich) = (ich-1)*k+1
        novers(ich) = nover
        ixyso(ich) = (ich-1)*nover + 1
      enddo
      ixys(nch+1) = npts+1
      ixyso(nch+1) = npts_over+1
      ra = sum(wover)
      ra2 = sum(qwts)
      call prin2('Error in perimeter of circle=*',abs(ra-2*pi),1)
      call prin2('Error in perimeter of circle=*',abs(ra2-2*pi),1)


      call prin2('srccoefs=*',srccoefs,96)

      isrc = 5
      itarg = 3


      rdotns = -0.5d0
      rdotnt = 0.5d0


      ndd = 0
      ndi = 0
      ndz = 3
      
      call h2d_comb_stab(srcinfo(1,isrc),8,srcinfo(1,itarg),
     1         rdotns,rdotnt,ndd,dpars,ndz,zpars,ndi,ipars,z1)

      call h2d_comb(srcinfo(1,isrc),8,srcinfo(1,itarg),
     1         ndd,dpars,ndz,zpars,ndi,ipars,z2)
      
      call prin2('stab kernel eval=*',z1,2)
      call prin2('kernel eval=*',z2,2)
      call prin2('error in stabilized kernel evaluator=*',abs(z1-z2),1)


c
c
c  test self quad at the first patch
c
c

c
c  test near quad at second patch from first patch
c
c
      nnz = 3*k*nch
      nquad = nnz*k

      allocate(row_ptr(npts+1),col_ind(nnz),wnear(nquad))
      allocate(wnearcoefs(nquad))

      allocate(ich_id(npts),ts_pts(npts))
      call get_chunk_id_ts(nch,norders,ixys,iptype,npts,ich_id,ts_pts)

      print *, "starting trid quad"

      call get_helm_dir_trid_quad_corr(zk,nch,k,npts,srcinfo,
     1   ixys,ich_id,ts_pts,ndz,zpars,nnz,row_ptr,col_ind,nquad,
     2   wnear,wnearcoefs)
      allocate(iquad(nnz+1))
      call prinf('nch=*',nch,1)
      call prinf('npts=*',npts,1)
      call prinf('nnz=*',nnz,1)
c      call prinf('ixys=*',ixys,nch+1)
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,iquad)
      iquadtype = 1
      eps = 0.51d-11
c      call prinf('iptype=*',iptype,nch)
c      call prin2('srcvals=*',srcinfo,8*npts)
c      call prin2('srccoefs=*',srccoefs,6*npts)
c      call prinf('ich_id=*',ich_id,npts)
c      call prin2('ts_pts=*',ts_pts,npts)
c      call prin2('eps=*',eps,1)
c      call prin2('zpars=*',zpars,6)
c      call prinf('iquadtype=*',iquadtype,1)
c      call prinf('nnz=*',nnz,1)
c      call prinf('row_ptr=*',row_ptr,npts+1)
c      call prinf('col_ind=*',col_ind,nnz)
c      call prinf('iquad=*',iquad,nnz+1)
c      call prinf('nquad=*',nquad,1)
c      call getnearquad_helm_comb_dir_2d(nch,norders,ixys,iptype,
c     1  npts,srccoefs,srcinfo,8,npts,srcinfo,ich_id,ts_pts,
c     2  eps,zpars,iquadtype,nnz,row_ptr,col_ind,iquad,nquad,
c     3  wnear)

      call prin2('zpars=*',zpars,6)
      call lpcomp_helm_comb_dir_addsub_2d(nch,norders,ixys,
     1  iptype,npts,srccoefs,srcinfo,8,npts,srcinfo,eps,
     2  zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,
     3  novers,npts_over,ixyso,srcover,wover,pot)
      
      erra = 0
      ra = 0
      do i=1,npts
        pot(i) = pot(i) + sigma(i)/2
        erra = erra + abs(pot(i)-potex(i))**2
        ra = ra + abs(potex(i))**2
        if(i.le.5) print *, i,real(pot(i)),real(potex(i)),
     1     real(pot(i))/real(potex(i))
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in pot=*',erra,1)

c      call lpcomp_galerkin_helm2d(nch,k,ixys,npts,
c     1  srcinfo,8,npts,srcinfo,eps,zpars,nnz,row_ptr,col_ind,iquad,
c     2  nquad,wnear,sigma,nover,npts_over,ixyso,srcover,wover,pot)

       

      return
      end




      subroutine get_helm_dir_trid_quad_corr(zk,nch,k,npts,srcinfo,ixys,
     1   ich_id,ts_pts,ndz,zpars,nnz,row_ptr,col_ind,nquad,wnear,
     2   wnearcoefs)
      implicit real *8 (a-h,o-z)
      complex *16 zk
      integer nch,k,npts
      real *8 srcinfo(8,npts)
      integer ixys(nch+1),ich_id(npts)
      real *8 ts_pts(npts)
      integer nnz,nquad,row_ptr(npts+1),col_ind(nnz)
      complex *16 wnear(nquad),wnearcoefs(nquad)
      real *8, allocatable :: tadj(:),wadj(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: srcoverslf(:,:),woverslf(:)
      real *8, allocatable :: xintermat(:,:),umat(:,:),vmat(:,:)
      real *8, allocatable :: ts(:),wts(:),pmat(:,:)
      real *8, allocatable :: tslf0(:),wslf0(:)
      real *8, allocatable :: tslf(:,:),wslf(:,:),xslfmat(:,:,:)
      real *8, allocatable :: pslfmat(:,:,:)
      complex *16, allocatable :: zints(:,:),zints2(:,:)
      complex *16, allocatable :: zpmat(:,:), zpslfmat(:,:,:)
      complex *16, allocatable :: fkern(:,:), fkernslf(:)
      complex *16 zpars(3),ima,zalpha,zbeta
      data ima/(0.0d0,1.0d0)/
c
c  get legendre nodes and weights
c

      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))

      itype = 2
      call legeexps(itype,k,ts,umat,vmat,wts)
      

      do i=1,npts+1
        row_ptr(i) = (i-1)*3+1
      enddo
c
c  set up row_ptr and col_ind for this case
c
      do ich=1,nch
        il = ich-1
        ir = ich+1
        if(il.le.0) il = nch
        if(ir.gt.nch) ir = 1
        do j=1,k
          ipt = (ich-1)*k + j
          col_ind(row_ptr(ipt)) = il
          col_ind(row_ptr(ipt)+1) = ich
          col_ind(row_ptr(ipt)+2) = ir
        enddo
      enddo
      print *, "done computing row_ptr and col_ind"
cc      call prinf('row_ptr=*',row_ptr,npts+1)
cc      call prinf('col_ind=*',col_ind,nnz)

c
c  Compute diagonal interpolation matrices 
c
      nslf = 24
      allocate(tslf0(nslf),wslf0(nslf))
      call load_selfquad_ipv0_iord3(tslf0,wslf0,nslf0)
      print *, "nslf=",nslf
      print *, "nslf0=",nslf0

      print *, "done loading self quad"


      allocate(tslf(2*nslf0,k),wslf(2*nslf0,k))
      allocate(xslfmat(k,2*nslf0,k))
      allocate(pslfmat(k,2*nslf0,k))
      allocate(zpslfmat(k,2*nslf0,k))
      alpha = 1.0d0
      beta = 0.0d0
      zalpha = 1.0d0
      zbeta = 0.0d0

      tl = -1.0d0
      tr = 1.0d0
      do inode=1,k
        tm = ts(inode)
        do l=1,nslf0
          tslf(l,inode) = (tslf0(l)+1)*(tm-tl)/2 + tl
          wslf(l,inode) = wslf0(l)*(tm-tl)/2
          tslf(l+nslf0,inode) = (tslf0(l)+1)*(tr-tm)/2 + tm
          wslf(l+nslf0,inode) = wslf0(l)*(tr-tm)/2
          call legepols(tslf(l,inode),k-1,pslfmat(1,l,inode))
          call legepols(tslf(l+nslf0,inode),k-1,
     1        pslfmat(1,l+nslf0,inode))
        enddo
        call dgemm('t','n',k,2*nslf0,k,alpha,umat,k,
     1     pslfmat(1,1,inode),k,beta,xslfmat(1,1,inode),k)
      enddo
c 
c  note vector assign
c  

      zpslfmat = pslfmat
      print *, "done computing self mat"

c
c  Compute off diagonal interpolation matrices 
c
      m = 20
      iref = min(3,ceiling(2*log(k+0.0d0)/log(2.0d0)))
      iref = 4
      print *, "iref=",iref
      nmid = 5
      nchquadadj = 2*(iref+1) + nmid
      
      mquad = m*nchquadadj

      allocate(tadj(mquad),wadj(mquad))
      call get_lege_adj_quad(m,nmid,iref,mquad,tadj,wadj)
      print *, "done getting adj quad"
      
      allocate(xintermat(k,mquad),pmat(k,mquad))
      allocate(zpmat(k,mquad))
      do j=1,mquad
        call legepols(tadj(j),k-1,pmat(1,j))
        zpmat(1:k,j) = pmat(1:k,j)
      enddo
      call dgemm('t','n',k,mquad,k,alpha,umat,k,pmat,k,beta,xintermat,
     1   k)
      print *, "done getting xintermat"
      
c
c
c
c  zints(i,j) are integrals against polynomials P_{i}(t), if
c     j \in[1,k] then target on self, if j\in [k+1,2*k], then
c     target on left panel, and if j \in[2*k+1,3*k] then target
c     on right panel
c
      allocate(zints(k,3*k),zints2(k,3*k))
      ndd = 0
      ndi = 0

      allocate(srcover(8,mquad),wover(mquad))
      allocate(fkern(mquad,2*k))
      allocate(fkernslf(2*nslf0))
      allocate(srcoverslf(8,2*nslf0),woverslf(2*nslf0))
      rdotns = -0.5d0
      rdotnt = 0.5d0
      do ich=1,nch
        istart = (ich-1)*k+1
        il = ich-1
        ir = ich+1
        if(il.le.0) il = nch
        if(ir.gt.nch) ir = 1

c  start self quadrature now
        do inode=1,k
          itarg = (ich-1)*k + inode
          call dgemm('n','n',6,2*nslf0,k,alpha,srcinfo(1,istart),8,
     1       xslfmat(1,1,inode),k,beta,srcoverslf,8)
          do j=1,2*nslf0
            ds = sqrt(srcoverslf(3,j)**2 + srcoverslf(4,j)**2)
            srcoverslf(7,j) = srcoverslf(4,j)/ds
            srcoverslf(8,j) = -srcoverslf(3,j)/ds
            woverslf(j) = ds*wslf(j,inode)
            call h2d_comb_stab(srcoverslf(1,j),8,srcinfo(1,itarg),
     1         rdotns,rdotnt,ndd,dpars,ndz,zpars,ndi,ipars,fkernslf(j)) 
            fkernslf(j) = fkernslf(j)*woverslf(j)
          enddo
          call zgemv('n',k,2*nslf0,zalpha,zpslfmat(1,1,inode),k,
     1       fkernslf,1,zbeta,zints(1,inode),1)
        enddo
c
c  start off diagonal quadrature now
c
        call dgemm('n','n',6,mquad,k,alpha,srcinfo(1,istart),8,
     1    xintermat,k,beta,srcover,8)
        do j=1,mquad
          ds = sqrt(srcover(3,j)**2 + srcover(4,j)**2)
          srcover(7,j) = srcover(4,j)/ds
          srcover(8,j) = -srcover(3,j)/ds
          wover(j) = ds*wadj(j)
        enddo

        do j=1,k
          itargl = (il-1)*k + j 
          itargr = (ir-1)*k + j
          do l=1,mquad
            call h2d_comb_stab(srcover(1,l),8,srcinfo(1,itargl),rdotns,
     1         rdotnt,ndd,dpars,ndz,zpars,ndi,ipars,fkern(l,j))
            fkern(l,j) = fkern(l,j)*wover(l)
            call h2d_comb_stab(srcover(1,l),8,srcinfo(1,itargr),rdotns,
     1         rdotnt,ndd,dpars,ndz,zpars,ndi,ipars,fkern(l,j+k)) 
            fkern(l,j+k) = fkern(l,j+k)*wover(l)
          enddo
        enddo
        call zgemm('n','n',k,2*k,mquad,zalpha,zpmat,k,fkern,mquad,
     1     zbeta,zints(1,k+1),k)
        call zrmatmatt(3*k,k,zints,k,umat,zints2)
c
c
c   insert quadrature at correct location
c
        do j=1,k
          itarg = (ich-1)*k + j
          icind = (row_ptr(itarg)+1-1)*k + 1
          do l=1,k
            wnearcoefs(icind + l-1) = zints(l,j)
            wnear(icind + l-1) = zints2(l,j)
          enddo
          itarg = (il-1)*k+j
          icind = (row_ptr(itarg)-1)*k+1
          do l=1,k
            wnearcoefs(icind + l-1) = zints(l,j+k)
            wnear(icind + l-1) = zints2(l,j+k)
          enddo
          itarg = (ir-1)*k+j
          icind = (row_ptr(itarg)+2-1)*k+1
          do l=1,k
            wnearcoefs(icind + l-1) = zints(l,j+2*k)
            wnear(icind + l-1) = zints2(l,j+2*k)
          enddo
        enddo
      enddo


      return
      end

      

      subroutine get_lege_adj_quad(m,nmid,iref,mquad,t,w)
      implicit real *8 (a-h,o-z)
      real *8 t(mquad),w(mquad),umat,vmat
      real *8, allocatable :: ts(:),ws(:)
      real *8, allocatable :: tchse(:)

      itype =1
      allocate(ts(m),ws(m))
      call legeexps(itype,m,ts,umat,vmat,ws)

      nch = 2*(iref+1) + nmid
      allocate(tchse(nch+1))

      tchse(1) = -1.0d0
      hpan = 2.0d0/(nmid+2)
      rpan = hpan*(0.5d0)**(iref)
      tchse(2) = tchse(1) + rpan
      do i=2,iref+1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan*2
      enddo
      
      do i=1,iref+2,iref+1+nmid
        tchse(i+1) = tchse(i) + hpan
      enddo
      rpan = hpan/2
      do i=iref+2,nch-1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan/2
      enddo
      tchse(nch+1) = 1.0d0

      do ich=1,nch
        ipstart = (ich-1)*m
        h = tchse(ich+1) - tchse(ich)
        do i=1,m
          t(ipstart+i) = tchse(ich) + h/2*(ts(i)+1)
          w(ipstart+i) = ws(i)*h/2
        enddo
      enddo
      
      
      return
      end
