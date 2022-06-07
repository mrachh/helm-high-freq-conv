
    kh = 20.0;
    ppw = 40.0;
    
    rfac = 0.3;
    nosc = 3;
    
    addpath ~/git/finufft/matlab/
    
   
   
    %variables
    
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
    
    
    
    ts = linspace(0,2*pi,nch+1);
    
    ab = zeros(2,nch);
    ab(1,:) = ts(1:end-1);
    ab(2,:) = ts(2:end);
    
    npts = chnkr.k*chnkr.nch;
    tvals = zeros(npts,1);
    
    [xs,ws] = lege.exps(chnkr.k);
    
    
    for i = 1:nch
        a=ab(1,i);
        b=ab(2,i);

        istart = (i-1)*chnkr.k + 1;
        iend = i*chnkr.k;
        tvals(istart:iend) = a + (b-a)*(xs+1)/2;
    end
    
    
    
    
    wts = weights(chnkr);
    
    [x,w,u,v] = lege.exps(chnkr.k);
    sol_c = reshape(sol_c,[chnkr.k chnkr.nch]);
    sol_c_coeff = u*sol_c;
    
    nfine = 1000;
    tfine  = 0:pi/nfine:2*pi-pi/nfine; %discretization of boundary
    tfine = tfine(:);
    
    
    h = ab(2,1) - ab(1,1);
    ichuse = floor(tfine/h) + 1;
    ichuse(ichuse>nch) = nch;
    ttuse = (tfine - (ichuse-1)*h)/h*2-1;
    ttuse = ttuse(:);
   
    ichuse = ichuse(:);
    sol_c_coeffuse = sol_c_coeff(:,ichuse);
    
    pols = lege.pols(ttuse,chnkr.k-1);
    pval = sol_c_coeffuse.*pols;
    sol_c_interp = sum(pval,1);
    sol_c_interp = sol_c_interp(:);
    sol_c = sol_c(:);
    
    
    
    
    figure(1)
    clf
    plot(tvals,real(sol_c),'k-'); hold on
    plot(tfine,real(sol_c_interp),'r.')
    
    
    return
    
    
    
    figure
    semilogy(abs(sol_c_hat)); hold on;
    semilogy(abs(rr1))

    
    rnorm = norm(mat);
    rnorm_inv = norm(mat_inv);
    
