c   Note: this file has many changes from helm_galerkin_routs.f
c   which had some version of hybrid galerkin/collocation scheme
c
c   trid_quad_corr only returns wnearcoefs, and only uses
c   srcinfo as provided
c
c   similar changes for the guru routine
c
c   helm_comb_dir_galerkin_solver2d has also been appropriately
c   changed. Do not try to use these routines with the test
c   codes in test, these will only work with the test codes
c   in paper_tests, and the examples in examples
c
c
c 
      subroutine get_helm_dir_trid_quad_corr(zk,nch,k,ngk,npts,nb,
     1   adjs,srcinfo,srccoefs,ndz,zpars,nnz,row_ptr,col_ind,
     2   nquad,wnearcoefs)
c
c  returns the galerkin quadrature corrections of order ngk 
c  where the geometry is discretized with kth order panels
c
c  the quadrature returned acts in coefficients space
c
c  Input arguments:
c    - zk: complex *16
c        Helmholtz wavenumber
c    - nch: integer
c        number of chunks
c    - k: integer
c        discretization order for geometry information
c    - ngk: integer
c        discretization order for basis functions
c    - npts: integer
c        total number of points on the boundary = nch*k
c    - nb: integer
c        total number of basis functions
c    - adjs: integer(2,nch)
c        adjacency information
c    - srcinfo: real *8(8,npts)
c        source info, srcinfo(1:2,:), xy
c                     srcinfo(3:4,:), dxydt
c                     srcinfo(5:6,:), d2xydt2
c                     srcinfo(7:8,:), normals info
c    - srccoefs: real *8(8,npts)
c        source coefficients
c    - ndz: integer
c        number of complex parameters (must be 3)
c    - zpars: complex *16(ndz)
c        complex kernel parameters, zpars(1) = zk, 
c        zpars(2) must be single layer strength
c        zpars(3) must be double layer strength
c    - nnz: integer
c        number of non-zero entries in row-sparse
c        compressed format as viewed from number of chunks
c        to number of chunks*number of basis functions
c  Output arguments:
c    - row_ptr: integer(nch+1)
c        fmm3dbie row pointer array
c    - col_ind: integer(nnz)
c    - nquad: integer
c    - wnearcoefs: complex *16(nquad)
c        quadrature correction array
c      
c 
c

      implicit real *8 (a-h,o-z)
      complex *16 zk
      integer nch,k,npts,ngk,nb
      integer adjs(2,nch)
      real *8 srcinfo(8,npts),srccoefs(6,npts)
      integer nnz,nquad,row_ptr(nb+1),col_ind(nnz)
      complex *16 wnearcoefs(nquad)
      real *8 dpars(1)
      complex *16 zpars(ndz)
      integer ipars
      procedure (), pointer :: fker,fkerstab
      external h2d_comb,h2d_comb_stab

      ndd = 0
      ndi = 0
      fker => h2d_comb
      fkerstab => h2d_comb_stab


      call get_helm_guru_trid_quad_corr(zk,nch,k,ngk,npts,nb,
     1   adjs,srcinfo,srccoefs,fker,fkerstab,ndd,dpars,
     2   ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,nquad,
     3   wnearcoefs)
      return
      end

c
c
c
c
c
      subroutine get_helm_guru_trid_quad_corr(zk,nch,k,ngk,npts,
     1   nb,adjs,srcinfo,srccoefs,fker,fkerstab,ndd,dpars,ndz,zpars,
     2   ndi,ipars,nnz,row_ptr,col_ind,nquad,wnear)
      implicit real *8 (a-h,o-z)
      complex *16 zk
      integer nch,k,npts
      integer adjs(2,nch)
      real *8 srcinfo(8,npts),srccoefs(6,npts)
      integer nnz,nquad,row_ptr(nb+1),col_ind(nnz)
      complex *16 wnear(nquad)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

      real *8, allocatable :: srcinfog(:,:),qwtsg(:),tsg(:)

      real *8, allocatable :: tadj(:),wadj(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: srcoverslf(:,:),woverslf(:)
      real *8, allocatable :: xintermat(:,:),umat(:,:),vmat(:,:)
      real *8, allocatable :: ts(:),wts(:),pmat(:,:)
      real *8, allocatable :: tslf0(:),wslf0(:)
      real *8, allocatable :: tslf(:,:),wslf(:,:),xslfmat(:,:,:)
      real *8, allocatable :: pslfmat(:,:,:)
      complex *16, allocatable :: zints(:,:), ztmp(:), zints2(:,:)
      complex *16, allocatable :: zpmat(:,:), zpslfmat(:,:,:)
      complex *16, allocatable :: zpmat2(:,:), zpslfmat2(:,:,:)
      complex *16, allocatable :: fkern(:,:), fkernslf(:)
      real *8, allocatable :: xint(:,:,:),rdlpcoefs(:,:),rspcoefs(:,:)
      real *8, allocatable :: rdotns_all(:),rdotnt_all(:)
      real *8, allocatable :: pmattmp(:,:),pols(:)

      complex *16 ima,zalpha,zbeta
      external fker,fkerstab
      data ima/(0.0d0,1.0d0)/

c
c  we use a slightly complicated approach here
c  first, we determine nodes on the surface which would
c  integrate log*smooth + smooth functions. Compute
c  integrals of legendre polynomials upto degree ngk
c  at those points, then integrate out the first ngk basis 
c  functions against the ngk basis functions for the outer
c  integral
c

c
c  initialize blas parameters
c
c
      alpha = 1.0d0
      beta = 0.0d0
      zalpha = 1.0d0
      zbeta = 0.0d0

c
c  get log*smooth + smooth nodes + weights 
c

      nslf = 24
      allocate(tslf0(nslf),wslf0(nslf))
      call load_selfquad_ipv0_iord3(tslf0,wslf0,nslf0)

      print *, "nslf=",nslf
      print *, "nslf0=",nslf0

      print *, "done loading self quad"

      kg = nslf
      allocate(tsg(kg),pmattmp(k,kg))
      do i=1,kg
        tsg(i) = tslf0(i)
        call legepols(tsg(i),k-1,pmattmp(1,i))
      enddo

      call prinf('kg=*',kg,1)
      call prin2('tsg=*',tsg,nslf)


c
c  set targets to be log*smooth + smooth nodes
c
c
      nptsg = nch*kg
      allocate(srcinfog(8,nptsg),qwtsg(nptsg))

      

      do ich=1,nch
        istart = (ich-1)*k + 1
        jstart = (ich-1)*kg + 1
        call dgemm('n','n',6,kg,k,alpha,srccoefs(1,istart),6,
     1       pmattmp,k,beta,srcinfog(1,jstart),8)
        do j=1,kg
          jpt = (ich-1)*kg + j
          ds = sqrt(srcinfog(3,jpt)**2 + srcinfog(4,jpt)**2)
          srcinfog(7,jpt) = srcinfog(4,jpt)/ds
          srcinfog(8,jpt) = -srcinfog(3,jpt)/ds
          qwtsg(jpt) = ds*wslf0(j)
        enddo
      enddo

c
c
c  initialize row_ptr, col_ind, and row_ptrg and col_indg
c
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))

      itype = 2
      call legeexps(itype,k,ts,umat,vmat,wts)
      
      
C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nb+1
        row_ptr(i) = (i-1)*3+1
      enddo
C$OMP END PARALLEL DO      

c
c  set up row_ptr and col_ind for this case
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ich,il,ir,j,ipt)
      do ich=1,nch
        il = adjs(1,ich)
        ir = adjs(2,ich)
        do j=1,ngk
          ipt = (ich-1)*ngk + j
          col_ind(row_ptr(ipt)) = il
          col_ind(row_ptr(ipt)+1) = ich
          col_ind(row_ptr(ipt)+2) = ir
        enddo
      enddo
C$OMP END PARALLEL DO      
      print *, "done computing row_ptr and col_ind"


c
c  Compute diagonal interpolation matrices 
c
      

      allocate(tslf(2*nslf0,kg),wslf(2*nslf0,kg))
c
c   xslfmat(i,j,l) i - basis function number, j
c   self quadrature node number, l - target number
c

      allocate(xslfmat(k,2*nslf0,kg))
      allocate(pslfmat(k,2*nslf0,kg))
      allocate(zpslfmat(k,2*nslf0,kg))
      allocate(zpslfmat2(ngk,2*nslf0,kg))
      print *, "k=",k
      print *, "kg=",kg
      print *, "nslf0=",nslf0
      tl = -1.0d0
      tr = 1.0d0
      print *, "here1"
      do inode=1,kg
        tm = tsg(inode)
        do l=1,nslf0
          tslf(l,inode) = (tslf0(l)+1.0d0)*(tm-tl)/2 + tl
          wslf(l,inode) = wslf0(l)*(tm-tl)/2
          tslf(l+nslf0,inode) = (tslf0(l)+1.0d0)*(tr-tm)/2 + tm
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
      zpslfmat2(1:ngk,1:2*nslf0,1:kg) = zpslfmat(1:ngk,1:2*nslf0,1:kg)
      print *, "done computing self mat"

c
c  Compute off diagonal interpolation matrices 
c
      m = 32
      iref = min(3,ceiling(2*log(k+0.0d0)/log(2.0d0)))
      print *, "iref=",iref
      nmid = 3
      nchquadadj = 2*(iref+1) + nmid
      
      mquad = m*nchquadadj

      allocate(tadj(mquad),wadj(mquad))
      call get_lege_adj_quad(m,nmid,iref,mquad,tadj,wadj)
      print *, "done getting adj quad"
      allocate(xintermat(k,mquad),pmat(k,mquad))
      allocate(zpmat(k,mquad),zpmat2(ngk,mquad))
      print *, "m=",m
      print *, "mquad=",mquad

      do j=1,mquad
        call legepols(tadj(j),k-1,pmat(1,j))
        zpmat(1:k,j) = pmat(1:k,j)
        zpmat2(1:ngk,j) = pmat(1:ngk,j)
      enddo
      call dgemm('t','n',k,mquad,k,alpha,umat,k,pmat,k,beta,xintermat,
     1   k)
      print *, "done getting xintermat"
      
c
c
c
c  zints(i,j) are integrals against polynomials P_{i}(t), if
c     j \in[1,kg] then target on self, if j\in [kg+1,2*kg], then
c     target on left panel, and if j \in[2*kg+1,3*kg] then target
c     on right panel
c
      allocate(zints(ngk,3*kg), zints2(3*kg,ngk), ztmp(ngk))

      allocate(srcover(8,mquad),wover(mquad))
      allocate(fkern(mquad,2*kg))
      allocate(fkernslf(2*nslf0))
      allocate(srcoverslf(8,2*nslf0),woverslf(2*nslf0))
      allocate(xint(k,k,kg),rdlpcoefs(k,kg),rspcoefs(k,kg))
      allocate(rdotns_all(2*nslf0),rdotnt_all(2*nslf0))

c
c  can't use legeinmt to generate rdotn stuff, since
c  tsg are not legendre nodes. Manually generate it here
c

cc      call legeinmt_allnodes(k,kg,xint)
      xint = 0
      allocate(pols(k))
      do ipt=1,kg
        tstart = tsg(ipt)
        do inode = 1,k
          tend = ts(inode)
          hh = tend-tstart
          do l=1,k
            tuse = (tend + tstart)/2 + hh/2*ts(l)
            call legepols(tuse,k-1,pols)
            do ipol=1,k
              xint(ipol,inode,ipt) = xint(ipol,inode,ipt) + 
     1              pols(ipol)*wts(l)*hh/2
            enddo
          enddo
        enddo
      enddo



C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ich,istart,istartg,il,ir)
C$OMP$PRIVATE(rdlpcoefs,rspcoefs,inode,itarg,srcoverslf)
C$OMP$PRIVATE(rdotns_all,rdotnt_all,ds,woverslf,fkernslf,mpt)
C$OMP$PRIVATE(srcover,j,wover,itargl,itargr,l,fkern,icind,zints)
C$OMP$PRIVATE(zints2,ztmp)
      do ich=1,nch
        istart = (ich-1)*k+1
        istartg = (ich-1)*kg+1
        il = adjs(1,ich)
        ir = adjs(2,ich)

        call chunk_to_ldlp_sp_xint(k,kg,ts,tsg,srcinfo(1,istart),
     1     srcinfog(1,istartg),umat,xint,rdlpcoefs,rspcoefs)
        

c  start self quadrature now
        do inode=1,kg
          itarg = (ich-1)*kg + inode
          call dgemm('n','n',6,2*nslf0,k,alpha,srcinfo(1,istart),8,
     1       xslfmat(1,1,inode),k,beta,srcoverslf,8)
          call dgemv('t',k,2*nslf0,alpha,pslfmat(1,1,inode),k,
     1       rdlpcoefs(1,inode),1,beta,rdotns_all,1)
          call dgemv('t',k,2*nslf0,alpha,pslfmat(1,1,inode),k,
     1       rspcoefs(1,inode),1,beta,rdotnt_all,1)
          do j=1,2*nslf0
            ds = sqrt(srcoverslf(3,j)**2 + srcoverslf(4,j)**2)
            srcoverslf(7,j) = srcoverslf(4,j)/ds
            srcoverslf(8,j) = -srcoverslf(3,j)/ds
            woverslf(j) = ds*wslf(j,inode)

            call fkerstab(srcoverslf(1,j),8,srcinfog(1,itarg),
     1         rdotns_all(j),rdotnt_all(j),ndd,dpars,ndz,zpars,
     2         ndi,ipars,fkernslf(j))
c            call h2d_comb(srcoverslf(1,j),8,srcinfog(1,itarg),
c     1         ndd,dpars,ndz,zpars,
c     2         ndi,ipars,fkernslf(j)) 
            fkernslf(j) = fkernslf(j)*woverslf(j)
          enddo
          call zgemv('n',ngk,2*nslf0,zalpha,zpslfmat2(1,1,inode),ngk,
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
        do j=1,kg
          itargl = (il-1)*kg + j 
          itargr = (ir-1)*kg + j

          do l=1,mquad
            call fker(srcover(1,l),8,srcinfog(1,itargl),
     1        ndd,dpars,ndz,zpars,ndi,ipars,fkern(l,j))
            fkern(l,j) = fkern(l,j)*wover(l)
            call fker(srcover(1,l),8,srcinfog(1,itargr),
     1        ndd,dpars,ndz,zpars,ndi,ipars,fkern(l,j+kg)) 
            fkern(l,j+kg) = fkern(l,j+kg)*wover(l)
          enddo
        enddo

        call zgemm('n','n',ngk,2*kg,mquad,zalpha,zpmat2,ngk,fkern,mquad,
     1     zbeta,zints(1,kg+1),ngk)
        do j=1,ngk
          do l=1,3*kg
            zints2(l,j) = zints(j,l)
          enddo
        enddo

        do j=1,ngk
          mpt = (ich-1)*kg+1
          ztmp = 0
          call get_galerkin_proj(2,ngk,kg,tsg,zints2(1,j),qwtsg(mpt),
     1       ztmp)
          do l=1,ngk
            itarg = (ich-1)*ngk+l
            icind = (row_ptr(itarg)+1-1)*ngk + j
            wnear(icind) = ztmp(l)
          enddo

          mpt = (il-1)*kg+1
          ztmp = 0
          call get_galerkin_proj(2,ngk,kg,tsg,zints2(kg+1,j),qwtsg(mpt),
     1      ztmp)
          do l=1,ngk
            itarg = (il-1)*ngk+l
            icind = (row_ptr(itarg)+2-1)*ngk + j
            wnear(icind) = ztmp(l)
          enddo

          mpt = (ir-1)*kg+1
          ztmp = 0
          call get_galerkin_proj(2,ngk,kg,tsg,zints2(2*kg+1,j),
     1      qwtsg(mpt),ztmp)
          do l=1,ngk
            itarg = (ir-1)*ngk+l
            icind = (row_ptr(itarg)-1)*ngk + j
            wnear(icind) = ztmp(l)
          enddo
        enddo
      enddo
c
c
c  now perform the integrals on the target patches, and insert
c  the quadrature at the correct location
c
c        do j=1,ngk
c          itarg = (ich-1)*ngk + j
c          icind = (row_ptr(itarg)+1-1)*ngk + 1
c          do l=1,ngk
c            wnear(icind+l-1) = 0
c            do m=1,kg
c              mpt = (ich-1)*kg + m
c              wnear(icind+l-1) = wnear(icind+l-1) + 
c     1           zints(l,m)*pmattmp(j,m)*qwtsg(mpt)
c            enddo
c          enddo
c          itarg = (il-1)*ngk + j
c          icind = (row_ptr(itarg)+2-1)*ngk + 1
c          do l=1,ngk
c            wnear(icind+l-1) = 0
c            do m=1,kg
c              mpt = (il-1)*kg + m
c              wnear(icind+l-1) = wnear(icind+l-1) + 
c     1           zints(l,m+kg)*pmattmp(j,m)*qwtsg(mpt)
c            enddo
c          enddo
c          itarg = (ir-1)*ngk + j
c          icind = (row_ptr(itarg)-1)*ngk + 1
c          do l=1,ngk
c            wnear(icind+l-1) = 0
c            do m=1,kg
c              mpt = (ir-1)*kg + m
c              wnear(icind+l-1) = wnear(icind+l-1) + 
c     1           zints(l,m+2*kg)*pmattmp(j,m)*qwtsg(mpt)
c            enddo
c          enddo
c        enddo
c      enddo
C$OMP END PARALLEL DO      


      return
      end
c
c
c
c
c
c
c
c
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
      
      do i=iref+2,iref+1+nmid
        tchse(i+1) = tchse(i) + hpan
      enddo
      rpan = hpan/2
      do i=iref+2+nmid,nch-1
        tchse(i+1) = tchse(i) + rpan
        rpan = rpan/2
      enddo
      tchse(nch+1) = 1.0d0

      do ich=1,nch
        ipstart = (ich-1)*m
        h = tchse(ich+1) - tchse(ich)
        do i=1,m
          t(ipstart+i) = tchse(ich) + h/2*(ts(i)+1.0d0)
          w(ipstart+i) = ws(i)*h/2
        enddo
      enddo
      
      
      return
      end
c
c
c
c
c
c
      subroutine lpcomp_galerkin_helm2d(nch,k,ngk,ixys,npts,nb,
     1  srcvals,adjs,qwts,eps,zpars,nnz,row_ptr,
     2  col_ind,iquad,nquad,wnear,sigmacoefs,nover,
     3  nptso,ixyso,srcover,whtsover,potcoefs)
c
c  This subroutine evaluates the helmholtz combined field
c  layer potential in the Galerkin formulation
c
c
      implicit real *8 (a-h,o-z)
      integer nch,k,ixys(nch+1),npts,ngk
      integer adjs(2,nch)
      real *8 srcvals(8,npts),eps, qwts(npts)
      complex *16 zpars(3)
      integer nnz,row_ptr(nb+1),col_ind(nnz),iquad(nnz+1)
      integer nquad
      complex *16 wnear(nquad),sigmacoefs(nb)
      integer nover,nptso,ixyso(nch+1)
      real *8 srcover(8,nptso),whtsover(nptso)
      complex *16 potcoefs(nb)

      integer norder,npols,npolso
      complex *16, allocatable :: potsort(:)


      real *8, allocatable :: sources(:,:),targvals(:,:)
      complex *16, allocatable :: charges(:),dipstr(:),sigmaover(:)
      complex *16, allocatable :: pot(:)
      real *8, allocatable :: dipvec(:,:)
      integer ns,nt
      real *8 dalpha,dbeta
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart


      integer ifaddsub

      integer ntj
      
      complex *16 zdotu,pottmp
      real *8 radexp,epsfmm
      real *8, allocatable :: tso(:),wso(:)
      real *8, allocatable :: ts(:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: pmat(:,:)
      real *8 umato,vmato

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: ctmp2(:),dtmp2(:)
      real *8, allocatable :: dipvec2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      real *8 over4pi
      integer nss,ii,l,npover
      integer nmax,ier,iper

      integer nd,ntarg0

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4

      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      itype = 1
      
      allocate(tso(nover),wso(nover),pot(npts))
      call legeexps(itype,nover,tso,umato,vmato,wso)


c
c    estimate max number of sources in neear field of 
c    any target
c
      nmax = 6*(k + nover)
      allocate(srctmp2(2,nmax),ctmp2(nmax),dtmp2(nmax))
      allocate(dipvec2(2,nmax))
           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(2,ns),targvals(2,npts))
      allocate(charges(ns),dipstr(ns),dipvec(2,ns))
      allocate(sigmaover(ns))

c 
c   evaluate density at oversampled nodes: assumes all
c   oversampled orders are identical
c
      allocate(pmat(ngk,nover))
      do i=1,nover
        call legepols(tso(i),ngk-1,pmat(1,i))
      enddo
      
      dalpha = 1.0d0
      dbeta = 0.0d0
      do i=1,nch
        istart = (i-1)*ngk + 1 
        istarto = ixyso(i)
        call dgemm('n','n',2,nover,ngk,dalpha,sigmacoefs(istart),2,pmat,
     1    ngk,dbeta,sigmaover(istarto),2)
      enddo
c
c       set relevatn parameters for the fmm
c
      alpha = zpars(2)
      beta = zpars(3)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)

        charges(i) = sigmaover(i)*whtsover(i)*alpha
        dipstr(i) = sigmaover(i)*whtsover(i)*beta
        dipvec(1,i) = srcover(7,i)
        dipvec(2,i) = srcover(8,i)
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        targvals(1,i) = srcvals(1,i)
        targvals(2,i) = srcvals(2,i)
      enddo
C$OMP END PARALLEL DO      

      ifcharge = 1
      ifdipole = 1

      if(alpha.eq.0) ifcharge = 0
      if(beta.eq.0) ifdipole = 0

c
c
c       call the fmm
c
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call hfmm2d(nd,eps,zpars(1),ns,sources,ifcharge,charges,
     1  ifdipole,dipstr,dipvec,iper,ifpgh,tmp,tmp,tmp,npts,
     1  targvals,ifpghtarg,pot,tmp,tmp,ier)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

            
      timeinfo(1) = t2-t1


c
c        compute threshold for ignoring local computation
c
      call get_fmm2d_thresh(2,ns,sources,2,npts,targvals,thresh)

c

ccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ich,ipt,i,j,jpatch,srctmp2)
ccC$OMP$PRIVATE(ctmp2,dtmp2,dipvec2,nss,l,jstart,ii,val)
      do ich=1,nch
        do i=1,k
          ipt = (ich-1)*k + i
          nss = 0
          do j=1,3
            if (j.eq.1) jpatch = ich
            if (j.eq.2) jpatch = adjs(1,ich)
            if (j.eq.3) jpatch = adjs(2,ich)

            do l=ixyso(jpatch),ixyso(jpatch+1)-1
              nss = nss+1
              srctmp2(1,nss) = srcover(1,l)
              srctmp2(2,nss) = srcover(2,l)

              if(ifcharge.eq.1) ctmp2(nss) = charges(l)
              if(ifdipole.eq.1) then
                dtmp2(nss) = dipstr(l)
                dipvec2(1,nss) = dipvec(1,l)
                dipvec2(2,nss) = dipvec(2,l)
              endif
            enddo
          enddo

          val = 0
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
            call h2d_directcp(nd,zpars(1),srctmp2,nss,ctmp2,
     1        targvals(1,ipt),1,val,thresh)
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
            call h2d_directdp(nd,zpars(1),srctmp2,nss,dtmp2,
     1          dipvec2,targvals(1,ipt),1,val,thresh)
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
            call h2d_directcdp(nd,zpars(1),srctmp2,nss,ctmp2,dtmp2,
     1          dipvec2,targvals(1,ipt),1,val,thresh)
          endif
          pot(ipt) = pot(ipt) - val
        enddo
      enddo
ccC$OMP END PARALLEL DO      
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1

      potcoefs = 0
c
c  project pot onto potcoefs now
c
      do ich=1,nch
        istart = (ich-1)*k + 1
        jstart = (ich-1)*ngk + 1
        call get_galerkin_proj(2,ngk,k,ts,pot(istart),qwts(istart),
     1        potcoefs(jstart))
      enddo

c
c       add in precomputed quadrature

      call cpu_time(t1)
C$      t1 = omp_get_wtime()
ccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
ccC$OMP$PRIVATE(jstart,l)
      do i=1,nb
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jquadstart = iquad(j)
          jstart = (jpatch-1)*ngk+1
          do l=1,ngk
             potcoefs(i) = potcoefs(i)+wnear(jquadstart+l-1)*
     1         sigmacoefs(jstart+l-1)
          enddo
        enddo
      enddo
ccC$OMP END PARALLEL DO


      
      return
      end
c
c
c
c
c
c

      subroutine helm_comb_dir_galerkin_solver2d(nch,k,ngk,ixys,
     1    npts,nb,srcvals,adjs,qwts,eps,zpars,nnz,row_ptr,col_ind,iquad,
     2    nquad,wnear,nover,npts_over,ixyso,srcover,wover,
     2    numit,ifinout,rhscoefs,eps_gmres,niter,errs,rres,solncoefs)
c
c  Solve the Helmholtz boundary value problem using the combined 
c  field integral equation
c
c  .. math ::
c  
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  Same ordered chunks, and same order oversampling
c
c  Input arguments:
c    - nch: integer
c        number of chunks
c    - k: integer
c        order of discretization
c    - ngk: integer
c        order of galerkin projection
c    - ixys: integer(nch+1)
c        starting location of data on patch i
c    - npts: integer
c        total number of discretization points on the boundary
c    - nb: integer
c        total number of basis functions
c    - srcvals: real *8 (8,npts)
c        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
c    - adjs: integer(2,nch)
c        adjacency info
c    - qwts: real *8 (npts)
c        quadrature weights for integrating smooth functions on the
c        curve
c    - eps: real *8
c        precision requested
c    - zpars: complex *16 (3)
c        kernel parameters (Referring to formula (1))
c        zpars(1) = k 
c        zpars(2) = alpha
c        zpars(3) = beta
c    - nnz: integer
c        number of non-zero entries in quadrature correction array
c    - row_ptr: integer(nb+1)
c        row_ptr(i) is the pointer to col_ind array where list of 
c        relevant source patches for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - nquad: integer
c        number of entries in wnear
c    - wnear: complex *16 (nquad)
c        quadrature corrections in coefficient land
c    - nover: integer
c        oversampling parameter
c    - npts_over: integer
c        total number of oversampled points
c    - ixyso: integer(nch+1)
c        starting location for oversampled data
c    - srcover: real *8 (8,npts_over)
c        total number of oversampled points
c    - wover: real *8(npts_over)
c        oversampled smooth quadrature weights
c    - numit: integer
c        max number of gmres iterations
c    - ifinout: integer
c        flag for interior or exterior problems (normals assumed to 
c        be pointing in exterior of region)
c        * ifinout = 0, interior problem
c        * ifinout = 1, exterior problem
c    - rhscoefs: complex *16(nb)
c        right hand side
c    - eps_gmres: real *8
c        gmres tolerance requested
c 
c  Output arguments:
c    - niter: integer
c        number of gmres iterations required for relative residual
c    - errs: real *8 (1:niter)
c        relative residual as a function of iteration number
c    - rres: real *8
c        relative residual for computed solution
c    - solncoefs: complex *16(nb)
c        density which solves the dirichlet problem
c-----------------------------------
c
      implicit none
      integer nch,k,npts,ixys(nch+1)
      integer ngk, nb
      real *8 srcvals(8,npts),eps,eps_gmres,qwts(npts)
      integer adjs(2,nch)
      complex *16 zpars(3)
      complex *16 rhscoefs(nb)
      complex *16 solncoefs(nb)
      integer nnz,row_ptr(nb+1),col_ind(nnz),iquad(nnz+1)
      integer nquad
      complex *16 wnear(nquad)
      integer nover,npts_over,ixyso(nch+1)
      real *8 srcover(8,npts_over),wover(npts_over)
      integer numit,ifinout


      real *8 errs(numit+1)
      real *8 rres,eps2
      integer niter


      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime
      real *8 ttot,done,pi

c
c
c       gmres variables
c
      complex *16 zid,ztmp
      real *8 rb,wnrm2
      integer it,iind,it1,l
      real *8 rmyerr
      complex *16 temp
      complex *16, allocatable :: vmat(:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)


      allocate(vmat(nb,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(nb),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4

c
c
c     start gmres code here
c
c     NOTE: matrix equation should be of the form (z*I + K)x = y
c       the identity scaling (z) is defined via zid below,
c       and K represents the action of the principal value 
c       part of the matvec
c
      zid = -(-1)**(ifinout)*zpars(3)/2


      niter=0

c
c      compute norm of right hand side and initialize v
c 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo


c
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
      do i=1,nb
        rb = rb + abs(rhscoefs(i))**2
      enddo
C$OMP END PARALLEL DO      
      rb = sqrt(rb)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nb
        vmat(i,1) = rhscoefs(i)/rb
      enddo
C$OMP END PARALLEL DO      

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c
        call lpcomp_galerkin_helm2d(nch,k,ngk,ixys,npts,nb,
     1    srcvals,adjs,qwts,eps,zpars,nnz,row_ptr,col_ind,iquad,
     2    nquad,wnear,vmat(1,it),nover,npts_over,ixyso,
     3    srcover,wover,wtmp)

        do l=1,it
          ztmp = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:ztmp)          
          do j=1,nb
            ztmp = ztmp + wtmp(j)*conjg(vmat(j,l))
          enddo
C$OMP END PARALLEL DO          
          hmat(l,it) = ztmp

C$OMP PARALLEL DO DEFAULT(SHARED) 
          do j=1,nb
            wtmp(j) = wtmp(j)-hmat(l,it)*vmat(j,l)
          enddo
C$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)        
        do j=1,nb
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
C$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

C$OMP PARALLEL DO DEFAULT(SHARED) 
        do j=1,nb
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo
C$OMP END PARALLEL DO        

        do l=1,it-1
          temp = cs(l)*hmat(l,it)+conjg(sn(l))*hmat(l+1,it)
          hmat(l+1,it) = -sn(l)*hmat(l,it)+cs(l)*hmat(l+1,it)
          hmat(l,it) = temp
        enddo

        ztmp = wnrm2

        call zrotmat_gmres2d(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+conjg(sn(it))*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

c
c            solve the linear system corresponding to
c            upper triangular part of hmat to obtain yvec
c
c            y = triu(H(1:it,1:it))\s(1:it);
c
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo



c
c          estimate x
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
          do j=1,nb
            solncoefs(j) = 0
            do i=1,it
              solncoefs(j) = solncoefs(j) + yvec(i)*vmat(j,i)
            enddo
          enddo
C$OMP END PARALLEL DO          


          rres = 0
C$OMP PARALLEL DO DEFAULT(SHARED)          
          do i=1,nb 
            wtmp(i) = 0
          enddo
C$OMP END PARALLEL DO          
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c

        call lpcomp_galerkin_helm2d(nch,k,ngk,ixys,npts,nb,
     1    srcvals,adjs,qwts,eps,zpars,nnz,row_ptr,col_ind,iquad,
     2    nquad,wnear,solncoefs,nover,npts_over,ixyso,
     3    srcover,wover,wtmp)

C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)            
          do i=1,nb
            rres = rres + abs(zid*solncoefs(i) + wtmp(i)-rhscoefs(i))**2
          enddo
C$OMP END PARALLEL DO          
          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo
c
      return
      end
c
c
c
c
c
c        

      subroutine get_galerkin_proj(nd,ngk,k,ts,vals,qwts,vcoefs)
c
c  this subroutine evaluates the constant/linear projection of a function
c  on an interval not parameterized by unit speed
c
c  Input arguments:
c    - k: integer
c        number of points on the panel
c    - ngk: integer
c        order of galerkin projection 
c    - ts: real *8 (k)
c        values in parameter space where function is sampled 
c    - vals: real *8 (nd,k)
c        values of the functions at the legendre nodes
c    - qwts: real *8 (k)
c        quadrature weights for integrating smooth functions
c
c  Output arguments:
c    - vcoefs: real *8 (nd,2)
c        coeffs of the linear projection
c       
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nd,k,ngk
      real *8, intent(in) :: ts(k),vals(nd,k),qwts(k)
      real *8, intent(out) :: vcoefs(nd,ngk)

      integer i,idim,j,l
      real *8 ra
      real *8 rhs(nd,ngk)
      real *8, allocatable :: amat(:,:),amatinv(:,:),pols(:)

      allocate(amat(ngk,ngk),amatinv(ngk,ngk), pols(ngk))

      do j=1,ngk
        do l=1,ngk
          amat(l,j) = 0
          amatinv(l,j) = 0
        enddo
        do idim=1,nd
          rhs(idim,j) = 0
        enddo
      enddo

      do i=1,k
        call legepols(ts(i), ngk-1, pols)
        do j=1,ngk
          do l=1,ngk
            amat(l,j) = amat(l,j) + pols(l)*pols(j)*qwts(i)
          enddo
          do idim=1,nd
            rhs(idim,j) = rhs(idim,j) + vals(idim,i)*pols(j)*qwts(i)
          enddo
        enddo
      enddo

      info = 0

      call dinverse(ngk, amat, info, amatinv)

      alpha = 1.0d0
      beta = 0.0d0

      call dgemm('n','t',nd,ngk,ngk,alpha,rhs,nd,
     1     amatinv,ngk,beta,vcoefs,nd)
      
      return
      end
c        
c
c
c
c
c
