function [err_out,pot_c,pot_ex,rnorm,rnorm_inv] = kress_circ_err(kh,ppw,nmax)

    %variables
    n  = ceil(kh*ppw);%half the number of points in the domain
    t  = 0:pi/n:2*pi-pi/n; %discretization of boundary

    %boundary
    xs   = cos(t);
    ys   = sin(t);
    dxs  = -sin(t);
    dys  = cos(t);
    ddxs = -cos(t);
    d3xs = sin(t);
    ddys = -sin(t);
    d3ys = -cos(t);
    ds   = sqrt(dxs.^2 + dys.^2);

    u_bdry = exp(1i*kh*xs);

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
    mat_inv = inv(mat);
    pot_c = mat_inv*transpose(2*u_bdry);
    
    rnorm = norm(mat);
    rnorm_inv = norm(mat_inv);
    

    i1 = (-nmax):1:(nmax);
    i1 = i1';
    expvals = exp(1j*i1*t);

    besseljs = besselj(i1,kh);
    besselhs = besselh(i1,kh);

    besselj_ders = (besselj(i1-1,kh) - besselj(i1+1,kh))/2;

    coeff = ((1j).^(i1-1)*2.*besseljs)/pi./(besselj_ders - 1j*besseljs)./besselhs;
    pot_ex = coeff.'*expvals/kh;
    pot_ex = pot_ex.';
    err_out = norm(pot_ex - pot_c)/norm(pot_ex);
end