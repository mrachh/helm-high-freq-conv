%function [ref_sols] = get_ref_sols(gpars,kh,bd_params,ppw)
% This subroutine computes the reference solutions 
% for plane wave incident data, and singular data
% given by r^(2*bd_params.spow+1), where spow was an integer
% 
%  Input arguments:
%    gpars.igeomtype: geometry type
%        igeomtype = 1, unit circle
%        igeomtype = 2, star shaped geometries (requires gpars.nosc, and
%          gpars.rfac) r(t) = 1 + rfac*cos(nosc*t)
%    kh: helmholtz wavenumber
%    bd_params: parameters for defining the boundary data
%         bd_params.dir - incident direction
%         bd_params.spow - power of singularity
%    ppw: integer
%      min points per wavelength
%
%  Output arguments:
%     ref_sols: reference solution struct
%        ref_sols.ts - t values of chunkie parameterization
%        ref_sols.wts - quadrature weights for integrating 
%                      smooth functions on chunkie discretization
%        ref_sols.pwsol - density corresponding to plane wave data
%        ref_sols.singsol - density corresponding to singular data
%    


   %if(nargin <=3) 
       ppw = 20;
   %end
   %if(nargin <= 2)
       bd_params.inc = 0;
       bd_params.spow = 0;
   %end
   %if(nargin <= 1)
       kh = 20.0;
   %end
   %if(nargin <= 0)
       gpars.igeomtype = 1;
       gpars.nosc = 3;
       gpars.rfac = 0.3;
   %end
   if(gpars.igeomtype == 1)
       nosc = 0;
       rfac = 0;
   end
   if(gpars.igeomtype == 2)
       nosc = gpars.nosc;
       rfac = gpars.rfac;
   end
   
   if(~isfield(bd_params,'inc'))
       bd_params.inc = 0;
   end
   if(~isfield(bd_params,'spow'))
       bd_params.spow = 0;
       
   end
   ctr = [0,0];
   phi = 0;
   scale = 1.0;
   fcurve = @(t) starfish(t,nosc,rfac,ctr,phi,scale);
   
   r = @(t) 1 + rfac*cos(nosc*t);
   drdt = @(t) -rfac*nosc*sin(nosc*t);
   dxdt = @(t) -r(t).*sin(t) + drdt(t).*cos(t);
   dydt = @(t) r(t).*cos(t) + drdt(t).*sin(t);
   dsdt = @(t) sqrt(dxdt(t).^2 + dydt(t).^2);
   L = integral(dsdt,0,2*pi);
   n  = ceil(kh*ppw*L/2/pi);%number of points in the domain
   clear r drdt dxdt dydt dsdt
   k = 16;
   nchmid_half = max(ceil(n/k/2),10);
   
   iref = 10;
   nchhalf = iref+1+nchmid_half;
   nch = 2*nchhalf;
   tchse = zeros(1,nchhalf+1);
   ab = zeros(2,nch);
   
   hpan = pi/(nchmid_half+1.0);
   rpan = hpan*(0.5d0)^(iref);
   tchse(1) = 0.0;
   tchse(2) = tchse(1) + rpan;
   for i=2:iref+1
       tchse(i+1) = tchse(i) + rpan;
       rpan = rpan*2;       
   end
   tchse((iref+3):(nchhalf+1)) = tchse(iref+2) + (1:(nchhalf-iref-1))*hpan;
   
    [xs,~] = lege.exps(k);
   
    ab(1,1:nchhalf) = tchse(1:nchhalf);
    ab(2,1:nchhalf) = tchse(2:nchhalf+1);
   
    ab(1,nchhalf+1:end) = -tchse(nchhalf+1:-1:2);
    ab(2,nchhalf+1:end) = -tchse(nchhalf:-1:1);
   
   
    adjs = zeros(2,nch);
    adjs(1,:) = (1:nch) - 1;
    adjs(2,:) = (1:nch) + 1;
    
    adjs(1,1)=nch;
    adjs(2,nch)=1;

    pref = chunkerpref();
    ta = 0;
    
    chnkr = chunker(pref); % empty chunker
    chnkr = chnkr.addchunk(nch);
    
    
    dim = checkcurveparam(fcurve,ta);
    pref.dim = dim;
    nout = 3;
    out = cell(nout,1);
    
    tvals = zeros(k,nch);



    for i = 1:nch
        a=ab(1,i);
        b=ab(2,i);

        ts = a + (b-a)*(xs+1)/2;
        tvals(:,i) = ts;
        [out{:}] = fcurve(ts);
        chnkr.r(:,:,i) = reshape(out{1},dim,k);
        chnkr.d(:,:,i) = reshape(out{2},dim,k);
        chnkr.d2(:,:,i) = reshape(out{3},dim,k);
        chnkr.h(i) = (b-a)/2;
    end

    chnkr.adj = adjs(:,1:nch);

    % added by Shidong Jiang
    chnkr.n = normals(chnkr);
    chnkr = chnkr.makedatarows(1);
    chnkr.data = tvals;
    
    bdpt = [1,0];

    
    
    fkern_chunkie = @(s,t) chnk.helm2d.kern(kh,s,t,'c',-kh);
    fkern = @(s,t) helm_circ_kern(kh,s,t,-kh);
    opdims(1) = 1; opdims(2) = 1;
    opts = [];
    opts.l2scale = true;
    mat = chunkermat(chnkr,fkern,opts);
    mat = 2*mat + eye(chnkr.k*chnkr.nch);
    mat_chunkie = chunkermat(chnkr,fkern_chunkie,opts);
    mat_chunkie = 2*mat_chunkie + eye(chnkr.k*chnkr.nch);
    
    
    xs = chnkr.r(1,:);
    xs = xs(:);
    xx = xs(:);
    ys = chnkr.r(2,:);
    ys = ys(:);
    yy = ys(:);
    inc = bd_params.inc;
    spow  = bd_params.spow;
    u_pw = exp(1i*kh*(xx*cos(inc) + yy*sin(inc)));
    rtmp = (xx-bdpt(1)).^2 + (yy-bdpt(2)).^2;
    u_sing = rtmp.^spow.*sqrt(rtmp);
    source = [0.01 0.3];
    rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
    rr = rr(:);
    u_test = (1i/4)*besselh(0,1,kh*rr);
    wts = weights(chnkr);
    u_pw = u_pw.*sqrt(wts(:));
    u_sing = u_sing.*sqrt(wts(:));
    u_test = u_test.*sqrt(wts(:));
    
    
    
    mat_inv = inv(mat);
    mat_inv_chunkie = inv(mat_chunkie);
    sol_test = mat_inv*2*u_test;
    sol_sing = mat_inv*u_sing;
    sol_pw = mat_inv*u_pw;
    sol_pw_chunkie = mat_inv_chunkie*u_pw;
    sol_test_chunkie = mat_inv_chunkie*2*u_test;
    
    sol_test = sol_test./sqrt(wts(:));
    sol_sing = sol_sing./sqrt(wts(:));
    sol_pw = sol_pw./sqrt(wts(:));
    sol_pw_chunkie = sol_pw_chunkie./sqrt(wts(:));
    sol_test_chunkie = sol_test_chunkie./sqrt(wts(:));
    tgt = [3.1 2.8];
    tgt = tgt(:);
    
    pot_test = chunkerkerneval(chnkr,fkern_chunkie,sol_test,tgt);
    pot_test_chunkie = chunkerkerneval(chnkr,fkern_chunkie,sol_test_chunkie,tgt);
    rr = sqrt((tgt(1)-source(1)).^2+(tgt(2)-source(2)).^2);
    rr = rr(:);
    pot_ex = (1i/4)*besselh(0,1,kh*rr);
    err_ex = norm(pot_ex-pot_test)/norm(sol_test.*sqrt(wts(:)));
    err_ex_chunkie = norm(pot_ex-pot_test_chunkie)/norm(sol_test.*sqrt(wts(:)));
    
    fprintf('Error in analytic solution test:%d\n',err_ex);
    fprintf('Error in analytic solution chunkie test:%d\n',err_ex_chunkie);
    
    figure(1)
    clf
    plot(chnkr,'k.'); hold on; quiver(chnkr);
    
%     figure(2)
%     clf
%     plot(tvals(:),real(sol_sing),'k.')
%     
    figure(3)
    clf
    plot(tvals(:),real(sol_pw),'r.'); hold on;
    plot(tvals(:),real(sol_pw_chunkie),'k.');
%     
%     figure(4)
%     clf
%     plot(tvals(:),real(u_sing./sqrt(wts(:))),'k.')
%     
%     figure(5)
%     clf
%     plot(tvals(:),real(u_pw./sqrt(wts(:))),'k.')
%     
    u_sing = u_sing./sqrt(wts(:));
    sol_use = reshape(sol_pw,[chnkr.k,chnkr.nch]);
    sol_use2 = reshape(sol_pw_chunkie,[chnkr.k,chnkr.nch]);
    rhs_use = reshape(u_sing,[chnkr.k,chnkr.nch]);
    [~,~,u,v] = lege.exps(chnkr.k);
    rhs_coefs = complex(zeros(size(rhs_use)));
    sol_coefs = complex(zeros(size(sol_use)));
    sol_coefs2 = complex(zeros(size(sol_use2)));
    for i=1:chnkr.nch
        sol_coefs(:,i) = u*sol_use(:,i);
        sol_coefs2(:,i) = u*sol_use2(:,i);
        rhs_coefs(:,i) = u*rhs_use(:,i);
    end
    
    
    rr = abs(sol_coefs(15,1:nch)) + abs(sol_coefs(16,1:nch));
    rr2 = abs(sol_coefs2(15,1:nch)) + abs(sol_coefs2(16,1:nch));
    rr3 = rr2.*(ab(2,1:nch)-ab(1,1:nch));
    
    figure(6)
    clf
    semilogy(1:nch,rr,'r.'); hold on;
    semilogy(1:nch,rr3,'k--'); hold on;
    semilogy(1:nch,rr2,'k.'); hold on;
    
    rr = abs(rhs_coefs(15,1:nch)) + abs(rhs_coefs(16,1:nch));
%     rr = rr.*(ab(2,1:nch)-ab(1,1:nch));
   
    semilogy(1:nch,rr,'b.') 
    
%     n_dir = 200;
%     t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
%     
%     sig_use = repmat((sol_sing.*wts(:)).', [n_dir,1]);
%     ffsig = exp(1j*kh*(cos(t_dir.')*xs.' + sin(t_dir.')*ys.')).*sig_use;
%     usig = sum(ffsig,2);
%     
%     
%     figure(8)
%     plot(t_dir,real(usig),'k.');
%     
    % run density and data interpolation test
   
%end
