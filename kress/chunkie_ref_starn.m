function [sol_c_ref,tref,dsdt_ref,wref,theta_ref,ffsig_ref] = chunkie_ref_starn(kh,rfac,nosc,nmin,ppw)
    r = @(t) 1 + rfac*cos(nosc*t);
    drdt = @(t) -rfac*nosc*sin(nosc*t);
    dxdt = @(t) -r(t).*sin(t) + drdt(t).*cos(t);
    dydt = @(t) r(t).*cos(t) + drdt(t).*sin(t);
    dsdt = @(t) sqrt(dxdt(t).^2 + dydt(t).^2);
    L = integral(dsdt,0,2*pi);
    n  = max(nmin,ceil(kh*ppw*L/2/pi));%number of points in the domain
    
    
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
    
    
    wref = weights(chnkr);
    
    mat_inv = inv(mat);
    sol_c_ref = mat_inv*2*u_bdry;
    
    
    ts = linspace(0,2*pi,nch+1);
    
    ab = zeros(2,nch);
    ab(1,:) = ts(1:end-1);
    ab(2,:) = ts(2:end);
    
    npts = chnkr.k*chnkr.nch;
    tref = zeros(npts,1);
    
    [xs,~] = lege.exps(chnkr.k);
    
    
    for i = 1:nch
        a=ab(1,i);
        b=ab(2,i);

        istart = (i-1)*chnkr.k + 1;
        iend = i*chnkr.k;
        tref(istart:iend) = a + (b-a)*(xs+1)/2;
    end
    
    dx = chnkr.d(1,:);
    dy = chnkr.d(2,:);
    dsdt_ref = sqrt(dx.^2 + dy.^2);
    
    theta_ref = 0;
    ffsig_ref = 0;
    
    
end