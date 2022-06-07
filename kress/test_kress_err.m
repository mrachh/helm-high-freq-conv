%function [err_out,pot_c,pot_ex,rnorm,rnorm_inv] = kress_err(kh,ppw,nmax)

    kh = 3.0;
    ppw = 100;
    
    rfac = 0.3;
    nosc = 3;
    %variables
   
    
    r = @(t) 1 + rfac*cos(nosc*t);
    drdt = @(t) -rfac*nosc*sin(nosc*t);
    dxdt = @(t) -r(t).*sin(t) + drdt(t).*cos(t);
    dydt = @(t) r(t).*cos(t) + drdt(t).*sin(t);
    dsdt = @(t) sqrt(dxdt(t).^2 + dydt(t).^2);
    L = integral(dsdt,0,2*pi);
    n  = ceil(kh*ppw*L/2/pi);%half the number of points in the domain
    
    clear r drdt dxdt dydt dsdt
    
    t  = 0:pi/n:2*pi-pi/n; %discretization of boundary
    
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
    src = zeros(8,2*n);
    src(1,:) = xs;
    src(2,:) = ys;
    src(3,:) = dxs;
    src(4,:) = dys;
    src(5,:) = ddxs;
    src(6,:) = ddys;
    src(7,:) = d3xs;
    src(8,:) = d3ys;


    %creating potential
    S  = slmat(kh,src,t);
    D  = dlmat(kh,src,t);

    %Combined layer
    %solving for the field 
    eta   = kh;
    
    mat = 2*(eye(2*n)/2 + D + 1i*eta*S);
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
    
    rr1 = fftshift(fft(ds));
    nmax = ceil(max(n/4,20));
    err_geom = norm(rr1(1:nmax))/norm(rr1);
    
    ms = length(sol_c);
    eps = 1e-15;
    isign = -1;
    sol_c_hat = finufft1d1(t,sol_c,isign,eps,ms)/ms;
    
    nfine = 1000;
    tfine  = 0:pi/nfine:2*pi-pi/nfine; %discretization of boundary
    isign = 1;
    sol_c_interp = finufft1d2(tfine,isign,eps,sol_c_hat);
    
    
    
    
    figure
    semilogy(abs(sol_c_hat)); hold on;
    semilogy(abs(rr1))

    
    rnorm = norm(mat);
    rnorm_inv = norm(mat_inv);
    
