function [sol_c_ref,tref,dsdt_ref] = kress_ref_larrycup(kh,a,b,nmin,ppw)
%     a = 0.3;
%     b = pi/3;
%     nmin = 50;
%     ppw = 15;
%     kh = 1.3;
    n0 = 600;
    src_info = larrycup(a,b,n0,n0);
    ds = sqrt(src_info(3,:).^2 + src_info(4,:).^2);
    L = sum(ds)*2*pi/n0;
    n  = max(nmin,ceil(kh*ppw*L/2/pi));%half the number of points in the domain
    
    if(mod(n,2))
        n = n+1;
    end
    src = larrycup(a,b,n,n0);
    
    nhalf = n/2;
    t  = 0:pi/nhalf:2*pi-pi/nhalf; %discretization of boundary
    tref = t;
    
    
    xs = src(1,:);
    ys = src(2,:);
    dxs = src(3,:);
    dys = src(4,:);
    ds = sqrt(dxs.^2 + dys.^2);
    dsdt_ref = ds(:);
    
    xx = xs(:);
    yy = ys(:);
    thet = pi/4;
    %u_bdry = exp(1i*kh*xx);
    u_bdry = exp(1i*kh*(xx*cos(thet) + yy*sin(thet)));
    source = [0.01 -1.01];
    rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
    rr = rr(:);
    %u_test = (1i/4)*besselh(0,1,kh*rr);
    

    %creating potential
    S  = slmat(kh,src,t);
    D  = dlmat(kh,src,t);

    %Combined layer
    %solving for the field 
    eta   = -kh;
    
    mat = 2*(eye(n)/2 + D + 1i*eta*S);
    sol_c_ref = mat\u_bdry;
    sol_c_ref = 2*sol_c_ref;
%     sol_test = mat\u_test;
%     sol_test = 2*sol_test;
%     
%     tgt = [3.1 2.8];
%     tgt = tgt(:);
%     pot_test = (dlmat_out(kh,src,tgt)+1i*eta*slmat_out(kh,src,tgt))*sol_test;
%     
%     
%     
%     rr = sqrt((tgt(1)-source(1)).^2+(tgt(2)-source(2)).^2);
%     rr = rr(:);
%     pot_ex = (1i/4)*besselh(0,1,kh*rr);
%     err_ex = norm(pot_ex-pot_test)/norm(sol_test.*sqrt(ds(:)));
%     
    
end