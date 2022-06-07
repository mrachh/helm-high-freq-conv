function [err_out,err_geom,err_ex,rnorm,rnorm_inv] = alpert_err_starn(kh,ppw,rfac,...
  nosc,norder,tref,sol_c_ref,dsdt_ref)

    %variables
    
    r = @(t) 1 + rfac*cos(nosc*t);
    drdt = @(t) -rfac*nosc*sin(nosc*t);
    dxdt = @(t) -r(t).*sin(t) + drdt(t).*cos(t);
    dydt = @(t) r(t).*cos(t) + drdt(t).*sin(t);
    dsdt = @(t) sqrt(dxdt(t).^2 + dydt(t).^2);
    L = integral(dsdt,0,2*pi);
    n  = ceil(kh*ppw*L/2/pi);%half the number of points in the domain
    
    clear r drdt dxdt dydt dsdt
    
    t  = 0:2*pi/n:2*pi-pi/n; %discretization of boundary
    
    rt = 1 + rfac*cos(nosc*t);
    drdt = -rfac*nosc*sin(nosc*t);
    d2rdt2 = -rfac*nosc*nosc*cos(nosc*t);
    d3rdt3 = rfac*nosc*nosc*nosc*sin(nosc*t);
    
    %boundary
    xs   = rt.*cos(t);
    ys   = rt.*sin(t);
    dxs  = -rt.*sin(t) + drdt.*cos(t);
    dys  = rt.*cos(t) + drdt.*sin(t);
    ddxs = -rt.*cos(t) - 2*drdt.*sin(t) + d2rdt2.*cos(t);
    ddys = -rt.*sin(t) + 2*drdt.*cos(t) + d2rdt2.*sin(t);
    d3xs = d3rdt3.*cos(t) - d2rdt2.*sin(t) - 2*(d2rdt2.*sin(t) + ... 
       drdt.*cos(t)) - (drdt.*cos(t) -rt.*sin(t));
    d3ys = d3rdt3.*sin(t) + d2rdt2.*cos(t) + 2*(d2rdt2.*cos(t) -drdt.*sin(t)) ...
           - (drdt.*sin(t) + rt.*cos(t));
    
    ds   = sqrt(dxs.^2 + dys.^2);
    
    
    

    xx = xs(:);
    u_bdry = exp(1i*kh*xx);
    source = [0.01 0.3];
    rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
    rr = rr(:);
    u_test = (1i/4)*besselh(0,1,kh*rr);

    %setting domain
    
    src_info = zeros(6,n);
    src_info(1,:) = xs;
    src_info(2,:) = ys;
    src_info(3,:) = dys./ds;
    src_info(4,:) = -dxs./ds;
    src_info(5,:) = ds;
    h = 2*pi/n;
    
    zpars = [kh -1j*kh 1];
    
    
    
    mat = 2*comb_ext_mat(zpars,norder,h,src_info);
    
    %mat = 2*(eye(2*n)/2 + D + 1i*eta*S);
    mat_inv = inv(mat);
    sol_c = mat_inv*2*u_bdry;
    sol_test = mat_inv*2*u_test;
    
    tgt = [3.1 2.8];
    tgt = tgt(:);
    
    zs = helm_c_p(kh,src_info,tgt);
    zd = helm_d_p(kh,src_info,tgt);
    z = zpars(2)*zs + zpars(3)*zd;
    pot_test = (z.*src_info(5,:))*sol_test*h;
    %pot_test = (dlmat_out(kh,src,tgt)+1i*eta*slmat_out(kh,src,tgt))*sol_test;
    
    rr = sqrt((tgt(1)-source(1)).^2+(tgt(2)-source(2)).^2);
    rr = rr(:);
    pot_ex = (1i/4)*besselh(0,1,kh*rr);
    err_ex = norm(pot_ex-pot_test)/norm(sol_test.*sqrt(ds(:)));
    
    
    
    
    ms = length(sol_c);
    eps = 1e-15;
    isign = -1;
    sol_c_hat = finufft1d1(t,sol_c,isign,eps,ms)/ms;
    
    dsdt_hat = finufft1d1(t,complex(ds),isign,eps,ms)/ms;
    
    
    isign = 1;
    sol_c_interp = finufft1d2(tref,isign,eps,sol_c_hat);
    dsdt_interp = finufft1d2(tref,isign,eps,dsdt_hat);
    
    err1 = (sol_c_interp(:) - sol_c_ref(:)).*sqrt(dsdt_ref(:));
    r1 = sol_c_ref(:).*sqrt(dsdt_ref(:));
    err_out = norm(err1)/norm(r1);
    
    err1 = (dsdt_interp(:) - dsdt_ref(:)).*sqrt(dsdt_ref(:));
    r1 = dsdt_ref(:).*sqrt(dsdt_ref(:));
    err_geom = norm(err1)/norm(r1);
    
    rnorm = norm(mat);
    rnorm_inv = norm(mat_inv);
end
    