      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcinfo(:,:),qwts(:),srccoefs(:,:)
      real *8, allocatable :: srcinfog(:,:),qwtsg(:)
      real *8, allocatable :: srcover(:,:),wover(:)
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
      integer, allocatable :: novers(:),ixyso(:)
      complex *16 zk,zpars(3),ima,z1,z2
      complex *16, allocatable :: sigmacoefsg(:),potcoefsg(:)
      complex *16, allocatable :: sigmacoefs(:),potcoefs(:)
      complex *16, allocatable :: potexcoefsg(:),potexcoefs(:)
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
c
c  k is the boundary discretization order
c
c  kg is the galerkin polynomial order
c
c  nover is the oversampled source order
c

      k = 20
      kg = 12
      nover = 24
      itype = 2
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      allocate(tsg(kg),umatg(kg,kg),vmatg(kg,kg),wtsg(kg))
      call legeexps(itype,kg,tsg,umatg,vmatg,wtsg)

      allocate(tover(nover),wtsover(nover))
      itype = 1
      call legeexps(itype,nover,tover,umo,vmo,wtsover)


      nch = 12
      allocate(srcinfo(8,k*nch),qwts(k*nch),qwtsg(kg*nch))
      allocate(srcinfog(8,kg*nch))
      allocate(srccoefs(6,k*nch))
      allocate(srcover(8,nover*nch),wover(nover*nch))
      allocate(norders(nch),iptype(nch),ixys(nch+1))
      allocate(novers(nch),ixyso(nch+1),ixysg(nch+1))
      h = 2*pi/(nch+0.0d0)
      npts = nch*k
      npts_over = nch*nover
      nptsg = nch*kg

      print *, "nptsg=",nptsg
      allocate(sigma(npts),sigmacoefsg(nptsg),pot(npts))
      allocate(sigmag(nptsg))
      allocate(potcoefsg(nptsg))
      allocate(sigmacoefs(npts),potcoefs(npts))
      allocate(potexcoefs(npts))

      allocate(potex(npts),potexg(nptsg),potexcoefsg(nptsg))
      do ich=1,nch
        tstart = (ich-1.0d0)*h
        tend = (ich+0.0d0)*h
c
c  get source info
c

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
        call dgemm('n','t',2,k,k,alpha,sigma(istart),
     1     2,umat,k,beta,sigmacoefs(istart),2)
        call dgemm('n','t',2,k,k,alpha,potex(istart),
     1     2,umat,k,beta,potexcoefs(istart),2)
c
c  get density info
c
        do j=1,kg
          ipt = (ich-1)*kg + j
          tuse = tstart + (tend-tstart)*(tsg(j)+1)/2
          srcinfog(1,ipt) = cos(tuse)
          srcinfog(2,ipt) = sin(tuse)
          srcinfog(3,ipt) = -sin(tuse)*h/2
          srcinfog(4,ipt) = cos(tuse)*h/2
          srcinfog(5,ipt) = -cos(tuse)*h*h/2/2
          srcinfog(6,ipt) = -sin(tuse)*h*h/2/2
          srcinfog(7,ipt) = cos(tuse)
          srcinfog(8,ipt) = sin(tuse)
          sigmag(ipt) = exp(ima*(imode+0.0d0)*tuse)
          potexg(ipt) = zfac*sigmag(ipt)
          qwtsg(ipt) = h/2*wtsg(j)
        enddo
        istart = (ich-1)*kg + 1
        call dgemm('n','t',2,kg,kg,alpha,sigmag(istart),
     1     2,umatg,kg,beta,sigmacoefsg(istart),2)
        call dgemm('n','t',2,kg,kg,alpha,potexg(istart),
     1     2,umatg,kg,beta,potexcoefsg(istart),2)

c 
c  get oversampled fun info
c    
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
        ixysg(ich) = (ich-1)*kg+1
      enddo
      ixys(nch+1) = npts+1
      ixyso(nch+1) = npts_over+1
      ixysg(nch+1) = nptsg+1
      ra = sum(wover)
      ra2 = sum(qwts)
      call prin2('Error in perimeter of circle=*',abs(ra-2*pi),1)
      call prin2('Error in perimeter of circle=*',abs(ra2-2*pi),1)
      call prin2('sigma=*',sigma,32)
      call prin2('sigmacoefs=*',sigmacoefs,32)


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

      nnzg = 3*kg*nch
      nquadg = nnz*kg

      allocate(row_ptr(npts+1),col_ind(nnz),wnear(nquad))
      allocate(wnearcoefs(nquad))
      allocate(row_ptrg(nptsg+1),col_indg(nnzg),wnearg(nquadg))
      allocate(wnearcoefsg(nquadg))



      print *, "starting trid quad"
      print *, "nch=",nch
      call get_helm_dir_trid_quad_corr(zk,nch,k,k,npts,npts,srcinfo,
     1   srcinfo,ndz,zpars,nnz,row_ptr,col_ind,nquad,wnear,wnearcoefs)

      call get_helm_dir_trid_quad_corr(zk,nch,k,kg,npts,nptsg,srcinfo,
     1   srcinfog,ndz,zpars,nnzg,row_ptrg,col_indg,nquadg,wnearg,
     2   wnearcoefsg)


      allocate(iquad(nnz+1),iquadg(nnzg+1))
      call prinf('nch=*',nch,1)
      call prinf('npts=*',npts,1)
      call prinf('nnz=*',nnz,1)
c      call prinf('ixys=*',ixys,nch+1)
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,iquad)

      call prinf('ixysg=*',ixysg,nch+1)
      stop
      call get_iquad_rsc2d(nch,ixysg,nptsg,nnzg,row_ptrg,col_indg,
     1   iquadg)
      iquadtype = 1
      eps = 0.51d-11

      call prin2('zpars=*',zpars,6)
      call lpcomp_helm_comb_dir_addsub_2d(nch,norders,ixys,
     1  iptype,npts,srccoefs,srcinfo,8,npts,srcinfo,eps,
     2  zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,
     3  novers,npts_over,ixyso,srcover,wover,pot)
      
      erra = 0
      ra = 0
      do i=1,npts
        pot(i) = pot(i) + sigma(i)/2*zpars(3)
        erra = erra + abs(pot(i)-potex(i))**2
        ra = ra + abs(potex(i))**2
        if(i.le.5) print *, i,real(pot(i)),real(potex(i)),
     1     real(pot(i))/real(potex(i))
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in pot=*',erra,1)

      potcoefs = 0
      call lpcomp_galerkin_helm2d(nch,k,ixys,npts,
     1  srcinfo,eps,ndz,zpars,nnz,row_ptr,col_ind,iquad,
     2  nquad,wnearcoefs,sigmacoefs,novers(1),npts_over,ixyso,
     3  srcover,wover,potcoefs)
      erra = 0
      ra = 0
      do i=1,npts
        potcoefs(i) = potcoefs(i) + sigmacoefs(i)/2*zpars(3)
        erra = erra + abs(potcoefs(i)-potexcoefs(i))**2
        ra = ra + abs(potexcoefs(i))**2
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in pot galerkin=*',erra,1)
      stop

      potcoefsg = 0
      call lpcomp_galerkin_helm2d(nch,kg,ixysg,nptsg,
     1  srcinfog,eps,ndz,zpars,nnzg,row_ptrg,col_indg,iquadg,
     2  nquadg,wnearcoefsg,sigmacoefsg,novers(1),npts_over,ixyso,
     3  srcover,wover,potcoefsg)
      erra = 0
      ra = 0
      do i=1,nptsg
        potcoefsg(i) = potcoefsg(i) + sigmacoefsg(i)/2*zpars(3)
        erra = erra + abs(potcoefsg(i)-potexcoefsg(i))**2
        ra = ra + abs(potexcoefsg(i))**2
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in pot galerkin=*',erra,1)

       

      return
      end




      subroutine get_helm_dir_trid_quad_corr(zk,nch,k,kg,npts,nptsg,
     1   srcinfo,srcinfog,ndz,zpars,nnz,row_ptr,
     2   col_ind,nquad,wnear,wnearcoefs)
      implicit real *8 (a-h,o-z)
      complex *16 zk
      integer nch,k,npts,kg,nptsg
      real *8 srcinfo(8,npts),srcinfog(8,nptsg)
      integer nnz,nquad,row_ptr(nptsg+1),col_ind(nnz)
      complex *16 wnear(nquad),wnearcoefs(nquad)
      real *8, allocatable :: tadj(:),wadj(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: srcoverslf(:,:),woverslf(:)
      real *8, allocatable :: xintermat(:,:),umat(:,:),vmat(:,:)
      real *8, allocatable :: ts(:),wts(:),pmat(:,:)
      real *8, allocatable :: tsg(:),wtsg(:),umatg(:,:),vmatg(:,:)
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

c
c  Assumption k>=kg
c
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))

      itype = 2
      call legeexps(itype,k,ts,umat,vmat,wts)
      
      allocate(tsg(kg),umatg(kg,kg),vmatg(kg,kg),wtsg(kg))
      call legeexps(itype,kg,tsg,umatg,vmatg,wtsg)
      

      do i=1,nptsg+1
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
        do j=1,kg
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

      allocate(tslf(2*nslf0,kg),wslf(2*nslf0,kg))
      allocate(xslfmat(k,2*nslf0,kg))
      allocate(pslfmat(k,2*nslf0,kg))
      allocate(zpslfmat(k,2*nslf0,kg))
      alpha = 1.0d0
      beta = 0.0d0
      zalpha = 1.0d0
      zbeta = 0.0d0

      tl = -1.0d0
      tr = 1.0d0
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
      allocate(zints(kg,3*kg),zints2(kg,3*kg))
      ndd = 0
      ndi = 0

      allocate(srcover(8,mquad),wover(mquad))
      allocate(fkern(mquad,2*kg))
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
        do inode=1,kg
          itarg = (ich-1)*kg + inode
          call dgemm('n','n',6,2*nslf0,k,alpha,srcinfo(1,istart),8,
     1       xslfmat(1,1,inode),k,beta,srcoverslf,8)
          do j=1,2*nslf0
            ds = sqrt(srcoverslf(3,j)**2 + srcoverslf(4,j)**2)
            srcoverslf(7,j) = srcoverslf(4,j)/ds
            srcoverslf(8,j) = -srcoverslf(3,j)/ds
            woverslf(j) = ds*wslf(j,inode)
            call h2d_comb_stab(srcoverslf(1,j),8,srcinfog(1,itarg),
     1         rdotns,rdotnt,ndd,dpars,ndz,zpars,ndi,ipars,fkernslf(j)) 
            fkernslf(j) = fkernslf(j)*woverslf(j)
          enddo
          call zgemv('n',kg,2*nslf0,zalpha,zpslfmat(1,1,inode),k,
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
            call h2d_comb_stab(srcover(1,l),8,srcinfog(1,itargl),rdotns,
     1         rdotnt,ndd,dpars,ndz,zpars,ndi,ipars,fkern(l,j))
            fkern(l,j) = fkern(l,j)*wover(l)
            call h2d_comb_stab(srcover(1,l),8,srcinfog(1,itargr),rdotns,
     1         rdotnt,ndd,dpars,ndz,zpars,ndi,ipars,fkern(l,j+kg)) 
            fkern(l,j+kg) = fkern(l,j+kg)*wover(l)
          enddo
        enddo
        call zgemm('n','n',kg,2*kg,mquad,zalpha,zpmat,k,fkern,mquad,
     1     zbeta,zints(1,kg+1),kg)
        call zrmatmatt(3*kg,kg,zints,k,umatg,zints2)
c
c
c   insert quadrature at correct location
c
        do j=1,kg
          itarg = (ich-1)*kg + j
          icind = (row_ptr(itarg)+1-1)*k + 1
          do l=1,kg
            wnearcoefs(icind + l-1) = zints(l,j)
            wnear(icind + l-1) = zints2(l,j)
          enddo
          itarg = (il-1)*kg+j
          icind = (row_ptr(itarg)+2-1)*kg+1
          do l=1,kg
            wnearcoefs(icind + l-1) = zints(l,j+kg)
            wnear(icind + l-1) = zints2(l,j+kg)
          enddo
          itarg = (ir-1)*kg+j
          icind = (row_ptr(itarg)-1)*kg+1
          do l=1,kg
            wnearcoefs(icind + l-1) = zints(l,j+2*kg)
            wnear(icind + l-1) = zints2(l,j+2*kg)
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
      subroutine lpcomp_galerkin_helm2d(nch,k,ixys,npts,
     1  srcvals,eps,ndz,zpars,nnz,row_ptr,
     2  col_ind,iquad,nquad,wnearcoefs,sigmacoefs,nover,
     3  nptso,ixyso,srcover,whtsover,potcoefs)
c
c  This subroutine evaluates the helmholtz combined field
c  layer potential in the Galerkin formulation
c
c
      implicit real *8 (a-h,o-z)
      integer nch,k,ixys(nch+1),npts,ndz
      real *8 srcvals(8,npts),eps
      complex *16 zpars(ndz)
      integer nnz,row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      integer nquad
      complex *16 wnearcoefs(nquad),sigmacoefs(npts)
      integer nover,nptso,ixyso(nch+1)
      real *8 srcover(8,nptso),whtsover(nptso)
      complex *16 potcoefs(npts)

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
      ntarg = npts
      done = 1
      pi = atan(done)*4

      itype = 2
      call prinf('k=*',k,1)
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      itype = 1
      print *, "nover=",nover
      
      allocate(tso(nover),wso(nover),pot(npts))
      call legeexps(itype,nover,tso,umato,vmato,wso)


c
c    estimate max number of sources in neear field of 
c    any target
c
      nmax = 0
      call get_near_corr_max2d(ntarg,row_ptr,nnz,col_ind,nch,
     1  ixyso,nmax)
      allocate(srctmp2(2,nmax),ctmp2(nmax),dtmp2(nmax))
      allocate(dipvec2(2,nmax))
           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(2,ns),targvals(2,ntarg))
      allocate(charges(ns),dipstr(ns),dipvec(2,ns))
      allocate(sigmaover(ns))
      print *, "ns=",ns
      print *, "ntarg=",ntarg

c 
c   evaluate density at oversampled nodes: assumes all
c   oversampled orders are identical
c
      allocate(pmat(k,nover))
      do i=1,nover
        call legepols(tso(i),k-1,pmat(1,i))
      enddo
      
      dalpha = 1.0d0
      dbeta = 0.0d0
      do i=1,nch
        istart = ixys(i)
        istarto = ixyso(i)
        call dgemm('n','n',2,nover,k,dalpha,sigmacoefs(istart),2,pmat,
     1    k,dbeta,sigmaover(istarto),2)
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
      do i=1,ntarg
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
     1  ifdipole,dipstr,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,
     1  targvals,ifpghtarg,pot,tmp,tmp,ier)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

            
      timeinfo(1) = t2-t1


c
c        compute threshold for ignoring local computation
c
      call get_fmm2d_thresh(2,ns,sources,2,ntarg,targvals,thresh)


c
c       add in precomputed quadrature

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixys(jpatch+1)-ixys(jpatch)
          jquadstart = iquad(j)
          jstart = ixys(jpatch) 
          do l=1,npols
             pot(i) = pot(i)+wnearcoefs(jquadstart+l-1)*
     1         sigmacoefs(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO

c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,dipvec2,nss,l,jstart,ii,val,npover)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
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
     1        targvals(1,i),1,val,thresh)
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          call h2d_directdp(nd,zpars(1),srctmp2,nss,dtmp2,
     1          dipvec2,targvals(1,i),1,val,thresh)
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call h2d_directcdp(nd,zpars(1),srctmp2,nss,ctmp2,dtmp2,
     1          dipvec2,targvals(1,i),1,val,thresh)
        endif
        pot(i) = pot(i) - val
      enddo
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


cc      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
cc      call prin2('time in lpcomp=*',ttot,1)
      do i=1,nch
        istart = ixys(i)
        call dgemm('n','t',2,k,k,dalpha,pot(istart),2,umat,k,dbeta,
     1     potcoefs(istart),2)
      enddo

      
      return
      end
