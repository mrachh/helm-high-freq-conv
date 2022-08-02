      subroutine get_helm_dir_trid_quad_corr(zk,nch,k,kg,npts,nptsg,
     1   adjs,srcinfo,srcinfog,ndz,zpars,nnz,row_ptr,col_ind,nquad,
     2   wnear,wnearcoefs)
      implicit real *8 (a-h,o-z)
      complex *16 zk
      integer nch,k,npts,kg,nptsg
      integer adjs(2,nch)
      real *8 srcinfo(8,npts),srcinfog(8,nptsg)
      integer nnz,nquad,row_ptr(nptsg+1),col_ind(nnz)
      complex *16 wnear(nquad),wnearcoefs(nquad)
      real *8 dpars(1)
      complex *16 zpars(ndz)
      integer ipars
      procedure (), pointer :: fker,fkerstab
      external h2d_comb,h2d_comb_stab

      ndd = 0
      ndi = 0
      fker => h2d_comb
      fkerstab => h2d_comb_stab

      call get_helm_guru_trid_quad_corr(zk,nch,k,kg,npts,nptsg,
     1   adjs,srcinfo,srcinfog,fker,fkerstab,ndd,dpars,
     2   ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,nquad,wnear,
     3   wnearcoefs)
      return
      end
c
c
c
c
c
c
      subroutine get_helm_guru_trid_quad_corr(zk,nch,k,kg,npts,nptsg,
     1   adjs,srcinfo,srcinfog,fker,fkerstab,ndd,dpars,ndz,zpars,
     2   ndi,ipars,nnz,row_ptr,col_ind,nquad,wnear,wnearcoefs)
      implicit real *8 (a-h,o-z)
      complex *16 zk
      integer nch,k,npts,kg,nptsg
      integer adjs(2,nch)
      real *8 srcinfo(8,npts),srcinfog(8,nptsg)
      integer nnz,nquad,row_ptr(nptsg+1),col_ind(nnz)
      complex *16 wnear(nquad),wnearcoefs(nquad)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

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
      real *8, allocatable :: xint(:,:,:),rdlpcoefs(:,:),rspcoefs(:,:)
      real *8, allocatable :: rdotns_all(:),rdotnt_all(:)

      complex *16 ima,zalpha,zbeta
      external fker,fkerstab
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
      
C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nptsg+1
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
        do j=1,kg
          ipt = (ich-1)*kg + j
          col_ind(row_ptr(ipt)) = il
          col_ind(row_ptr(ipt)+1) = ich
          col_ind(row_ptr(ipt)+2) = ir
        enddo
      enddo
C$OMP END PARALLEL DO      
      print *, "done computing row_ptr and col_ind"
cc      call prinf('row_ptr=*',row_ptr,npts+1)

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
      print *, "k=",k
      print *, "kg=",kg
      print *, "nslf0=",nslf0
      alpha = 1.0d0
      beta = 0.0d0
      zalpha = 1.0d0
      zbeta = 0.0d0

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
      allocate(zpmat(k,mquad))
      print *, "m=",m
      print *, "mquad=",mquad

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

      allocate(srcover(8,mquad),wover(mquad))
      allocate(fkern(mquad,2*kg))
      allocate(fkernslf(2*nslf0))
      allocate(srcoverslf(8,2*nslf0),woverslf(2*nslf0))
      allocate(xint(k,k,kg),rdlpcoefs(k,kg),rspcoefs(k,kg))
      allocate(rdotns_all(2*nslf0),rdotnt_all(2*nslf0))
      call legeinmt_allnodes(k,kg,xint)


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ich,istart,istartg,il,ir)
C$OMP$PRIVATE(rdlpcoefs,rspcoefs,inode,itarg,srcoverslf)
C$OMP$PRIVATE(rdotns_all,rdotnt_all,ds,woverslf,fkernslf)
C$OMP$PRIVATE(srcover,j,wover,itargl,itargr,l,fkern,icind,zints)
C$OMP$PRIVATE(zints2)
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
            call fker(srcover(1,l),8,srcinfog(1,itargl),
     1        ndd,dpars,ndz,zpars,ndi,ipars,fkern(l,j))
            fkern(l,j) = fkern(l,j)*wover(l)
            call fker(srcover(1,l),8,srcinfog(1,itargr),
     1        ndd,dpars,ndz,zpars,ndi,ipars,fkern(l,j+kg)) 
            fkern(l,j+kg) = fkern(l,j+kg)*wover(l)
          enddo
        enddo

        call zgemm('n','n',kg,2*kg,mquad,zalpha,zpmat,k,fkern,mquad,
     1     zbeta,zints(1,kg+1),kg)
        call zrmatmatt(3*kg,kg,zints,kg,umatg,zints2)
c
c
c   insert quadrature at correct location
c
        do j=1,kg
          itarg = (ich-1)*kg + j
          icind = (row_ptr(itarg)+1-1)*kg + 1
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
C$OMP END PARALLEL DO      


      return
      end
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
      subroutine lpcomp_galerkin_helm2d(nch,k,ixys,npts,
     1  srcvals,eps,zpars,nnz,row_ptr,
     2  col_ind,iquad,nquad,wnearcoefs,sigmacoefs,nover,
     3  nptso,ixyso,srcover,whtsover,potcoefs)
c
c  This subroutine evaluates the helmholtz combined field
c  layer potential in the Galerkin formulation
c
c
      implicit real *8 (a-h,o-z)
      integer nch,k,ixys(nch+1),npts
      real *8 srcvals(8,npts),eps
      complex *16 zpars(3)
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
      allocate(ts(k),umat(k,k),vmat(k,k),wts(k))
      call legeexps(itype,k,ts,umat,vmat,wts)

      itype = 1
      
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
c
c
c
c
c
c
      subroutine helm_comb_dir_galerkin_solver2d(nch,k,ixys,
     1    npts,srcvals,eps,zpars,nnz,row_ptr,col_ind,iquad,
     2    nquad,wnearcoefs,nover,npts_over,ixyso,srcover,wover,
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
c    - ixys: integer(nch+1)
c        starting location of data on patch i
c    - npts: integer
c        total number of discretization points on the boundary
c    - srcvals: real *8 (8,npts)
c        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
c    - eps: real *8
c        precision requested
c    - zpars: complex *16 (3)
c        kernel parameters (Referring to formula (1))
c        zpars(1) = k 
c        zpars(2) = alpha
c        zpars(3) = beta
c    - nnz: integer
c        number of non-zero entries in quadrature correction array
c    - row_ptr: integer(npts+1)
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
c    - wnearcoefs: complex *16 (nquad)
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
c        
c    - numit: integer
c        max number of gmres iterations
c    - ifinout: integer
c        flag for interior or exterior problems (normals assumed to 
c        be pointing in exterior of region)
c        * ifinout = 0, interior problem
c        * ifinout = 1, exterior problem
c    - rhscoefs: complex *16(npts)
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
c    - solncoefs: complex *16(npts)
c        density which solves the dirichlet problem
c-----------------------------------
c
      implicit none
      integer nch,k,npts,ixys(nch+1)
      real *8 srcvals(8,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhscoefs(npts)
      complex *16 solncoefs(npts)
      integer nnz,row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      integer nquad
      complex *16 wnearcoefs(nquad)
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


      allocate(vmat(npts,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(npts),svec(numit+1),yvec(numit+1))


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
      do i=1,npts
        rb = rb + abs(rhscoefs(i))**2
      enddo
C$OMP END PARALLEL DO      
      rb = sqrt(rb)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
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


        call lpcomp_galerkin_helm2d(nch,k,ixys,npts,
     1    srcvals,eps,zpars,nnz,row_ptr,col_ind,iquad,
     2    nquad,wnearcoefs,vmat(1,it),nover,npts_over,ixyso,
     3    srcover,wover,wtmp)

        do l=1,it
          ztmp = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:ztmp)          
          do j=1,npts
            ztmp = ztmp + wtmp(j)*conjg(vmat(j,l))
          enddo
C$OMP END PARALLEL DO          
          hmat(l,it) = ztmp

C$OMP PARALLEL DO DEFAULT(SHARED) 
          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(l,it)*vmat(j,l)
          enddo
C$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)        
        do j=1,npts
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
C$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

C$OMP PARALLEL DO DEFAULT(SHARED) 
        do j=1,npts
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
          do j=1,npts
            solncoefs(j) = 0
            do i=1,it
              solncoefs(j) = solncoefs(j) + yvec(i)*vmat(j,i)
            enddo
          enddo
C$OMP END PARALLEL DO          


          rres = 0
C$OMP PARALLEL DO DEFAULT(SHARED)          
          do i=1,npts
            wtmp(i) = 0
          enddo
C$OMP END PARALLEL DO          
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c

        call lpcomp_galerkin_helm2d(nch,k,ixys,npts,
     1    srcvals,eps,zpars,nnz,row_ptr,col_ind,iquad,
     2    nquad,wnearcoefs,solncoefs,nover,npts_over,ixyso,
     3    srcover,wover,wtmp)

C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)            
          do i=1,npts
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
