function [err_out,err_geom,err_ex,rnorm,rnorm_inv] = chunkie_err_starn(kh,ppw,rfac,...
  nosc,tref,sol_c_ref,dsdt_ref,wref)
    r = @(t) 1 + rfac*cos(nosc*t);
    drdt = @(t) -rfac*nosc*sin(nosc*t);
    dxdt = @(t) -r(t).*sin(t) + drdt(t).*cos(t);
    dydt = @(t) r(t).*cos(t) + drdt(t).*sin(t);
    dsdt = @(t) sqrt(dxdt(t).^2 + dydt(t).^2);
    L = integral(dsdt,0,2*pi);
    n  = ceil(kh*ppw*L/2/pi);%number of points in the domain
    
    
    clear r drdt dxdt dydt dsdt
    k = 16;
    nch = ceil(n/k);
    chnkr = chunkerfuncuni(@(t) starfish(t,nosc,rfac),nch);
    
    fkern = @(s,t) chnk.helm2d.kern(kh,s,t,'c',-1j*kh);
    opdims(1) = 1; opdims(2) = 1;
    opts = [];
    mat = chunkermat(chnkr,fkern,opts);
    mat = 2*mat + eye(chnkr.k*chnkr.nch);
    

    
    xs = chnkr.r(1,:);
    xs = xs(:);
    xx = xs(:);
    ys = chnkr.r(2,:);
    ys = ys(:);
    yy = ys(:);
    u_bdry = exp(1i*kh*xx);
    source = [0.01 0.3];
    rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
    rr = rr(:);
    u_test = (1i/4)*besselh(0,1,kh*rr);
    
    wts = weights(chnkr);
    
    mat_inv = inv(mat);
    sol_c = mat_inv*2*u_bdry;
    sol_test = mat_inv*2*u_test;
    
    tgt = [3.1 2.8];
    tgt = tgt(:);
    
    pot_test = chunkerkerneval(chnkr,fkern,sol_test,tgt);
    rr = sqrt((tgt(1)-source(1)).^2+(tgt(2)-source(2)).^2);
    rr = rr(:);
    pot_ex = (1i/4)*besselh(0,1,kh*rr);
    err_ex = norm(pot_ex-pot_test)/norm(sol_test.*sqrt(wts(:)));
   
    
    [~,~,u,~] = lege.exps(chnkr.k);
    sol_c = reshape(sol_c,[chnkr.k chnkr.nch]);
    dx = chnkr.d(1,:);
    dy = chnkr.d(2,:);
    dsdt = sqrt(dx.^2 + dy.^2);
    dsdt = reshape(dsdt, [chnkr.k chnkr.nch]);
    
    sol_c_coeff = u*sol_c;
    dsdt_coeff = u*dsdt;
    
    
    h = 2*pi/nch;
    ichuse = floor(tref/h) + 1;
    ichuse(ichuse>nch) = nch;
    ttuse = (tref - (ichuse-1)*h)/h*2-1;
    ttuse = ttuse(:);
   
    ichuse = ichuse(:);
    sol_c_coeffuse = sol_c_coeff(:,ichuse);
    dsdt_coeffuse = dsdt_coeff(:,ichuse);
    
    pols = lege.pols(ttuse,chnkr.k-1);
    pval = sol_c_coeffuse.*pols;
    sol_c_interp = sum(pval,1);
    sol_c_interp = sol_c_interp(:);
    dsdt_interp = sum(dsdt_coeffuse.*pols,1);
    dsdt_interp = dsdt_interp(:);
    
    err1 = (dsdt_interp(:)-dsdt_ref(:)).*sqrt(wref(:));
    r1 = dsdt_ref(:).*sqrt(wref(:));
    err_geom = norm(err1)/norm(r1);
    
    err1 = (sol_c_interp(:) - sol_c_ref(:)).*sqrt(wref(:));
    r1 = sol_c_ref(:).*sqrt(wref(:));
    err_out = norm(err1)/norm(r1);
    
    
    D = diag(sqrt(wts(:)));
    Dinv = diag(1.0/sqrt(wts(:)));
    rnorm = norm(D*mat*Dinv);
    rnorm_inv = norm(D*mat_inv*Dinv);
    
end