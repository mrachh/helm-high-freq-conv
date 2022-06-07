function [sol_c_ref,tref,dsdt_ref,theta_ref,ffsig_ref] = kress_ref_starn(kh,rfac,nosc,nmin,ppw)
    r = @(t) 1 + rfac*cos(nosc*t);
    drdt = @(t) -rfac*nosc*sin(nosc*t);
    dxdt = @(t) -r(t).*sin(t) + drdt(t).*cos(t);
    dydt = @(t) r(t).*cos(t) + drdt(t).*sin(t);
    dsdt = @(t) sqrt(dxdt(t).^2 + dydt(t).^2);
    L = integral(dsdt,0,2*pi);
    n  = max(nmin,ceil(kh*ppw*L/2/pi));%half the number of points in the domain
    
    clear r drdt dxdt dydt dsdt
    
    t  = 0:pi/n:2*pi-pi/n; %discretization of boundary
    tref = t;
    
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
    dsdt_ref = ds(:);

    xx = xs(:);
    u_bdry = exp(1i*kh*xx);
    source = [0.01 0.3];
    rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
    rr = rr(:);
    

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
    eta   = -kh;
    
    mat = 2*(eye(2*n)/2 + D + 1i*eta*S);
    sol_c_ref = mat\u_bdry;
    sol_c_ref = 2*sol_c_ref;  
    
    theta_ref = t(:);
    sk = kh*cos(theta_ref);
    tk = kh*sin(theta_ref);
    yy = ys(:);
    
    
    h = pi/n;
    rnx = dys(:)./ds(:);
    rny = -dxs(:)./ds(:);
    c1 = sol_c_ref.*dsdt_ref*h*1i/4;
    c2 = sol_c_ref.*dsdt_ref.*rnx*h*1i/4;
    c3 = sol_c_ref.*dsdt_ref.*rny*h*1i/4;
    ffsig_ref = zeros(size(sk));
    
%     isign = -1;
%     eps = 1e-15;
%     f1 = finufft2d3(xx,yy,c1,isign,eps,sk,tk);
%     f2 = finufft2d3(xx,yy,c2,isign,eps,sk,tk);
%     f3 = finufft2d3(xx,yy,c3,isign,eps,sk,tk);
%     
%     ffsig_ref = - f2.*sk - f3.*tk + eta*f1;
    
end