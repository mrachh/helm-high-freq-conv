function [err_out,err_geom,err_ex,rnorm,rnorm_inv] = kress_err_larrycup(kh,ppw,a,...
  b,tref,sol_c_ref,dsdt_ref)

    n0 = 600;
    src_info = larrycup(a,b,n0,n0);
    ds = sqrt(src_info(3,:).^2 + src_info(4,:).^2);
    L = sum(ds)*2*pi/n0;
    n  = ceil(kh*ppw*L/2/pi);
    
    if(mod(n,2))
        n = n+1;
    end
    src = larrycup(a,b,n,n0);
    xs = src(1,:);
    ys = src(2,:);
    dxs = src(3,:);
    dys = src(4,:);
    ds = sqrt(dxs.^2 + dys.^2);
    
    nhalf = n/2;
    t  = 0:pi/nhalf:2*pi-pi/nhalf; %discretization of boundary
    
    
    xx = xs(:);
    yy = ys(:);
    thet = pi/4;
    u_bdry = exp(1i*kh*(xx*cos(thet) + yy*sin(thet)));
    source = [0.01 -1.01];
    rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
    rr = rr(:);
    u_test = (1i/4)*besselh(0,1,kh*rr);

   

    %creating potential
    S  = slmat(kh,src,t);
    D  = dlmat(kh,src,t);

    %Combined layer
    %solving for the field 
    eta   = -kh;
    
    mat = 2*(eye(n)/2 + D + 1i*eta*S);
    mat_inv = inv(mat);
    sol_c = mat_inv*2*u_bdry;
    sol_test = mat_inv*2*u_test;
    
    tgt = [3.1 2.8];
    tgt = tgt(:);
    pot_test = (dlmat_out(kh,src,tgt)+1i*eta*slmat_out(kh,src,tgt))*sol_test;
    
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
    
    h = 2*pi/n;
    D = diag(sqrt(ds(:)*h));
    Dinv = diag(1.0/sqrt(ds(:)*h));
    rnorm = norm(D*mat*Dinv);
    rnorm_inv = norm(D*mat_inv*Dinv);
    

end