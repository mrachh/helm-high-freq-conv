function [qdiffq,qlnq] = get_circ_q_diff_ln_q(kh,ppw)
    k = kh;
    rfac = 6;
    mmax = ceil(20*rfac*k);
    muse = 2*mmax;

    nover = 6;
    n = (mmax + 2*k)*2*nover;
    ts = 0:2*pi/n:2*pi*(1-1/n);

    slp = -1j*k;
    dlp = 1;

    stt = sin(ts/2);
    fval = -slp*besselj(0,2*k*stt)+dlp*k*stt.*besselj(1,2*k*stt);
    % rescale by 2 pi to make identity term 1;
    fval = fval/2/pi;

    isign = -1;
    eps = 1e-15;
    fcoef = finufft1d1(ts,fval,isign,eps,muse);
    fcoef = fcoef/length(ts);


    gval = 1i*slp*besselh(0,1,2*k*stt)/2 - ...
        1i*dlp*k.*stt.*besselh(1,1,2*k*stt)/2 - ...
        log(4*stt.*stt).*fval;
    C  = 0.57721566490153286060651209008240243104215933593992; %euler-mascheroni's constant
    gval(1) = -dlp/2/pi + slp*(1i*pi-2*C-2*log(k/2))/2/pi;

    isign = -1;
    eps = 1e-15;
    gcoef = finufft1d1(ts,gval,isign,eps,muse);
    gcoef = gcoef/length(ts);
    
    n  = ceil(k*ppw);

    qdiffq = zeros(2*n,1);

    qind = -(n-1):1:n;
    qindinv = 1.0./abs(qind);
    qindinv(qindinv == Inf) = 0;
    qindinv_rep = repmat(qindinv,[1,4]);
    qinduse = [(qind-4*n) (qind-2*n) (qind+2*n) (qind+4*n)];

    pinduse = [-4*n -2*n 2*n 4*n];

    gammas_est = zeros(2*n,1);

    
    

    qgind = -(4*n-1):1:4*n;
    qgindinv = 1.0./abs(qgind);
    qgindinv(qgindinv == Inf) = 0;

    qqind = [-(4*n-1):-n n+1:4*n];
    qqindinv = 1.0./abs(qqind);


    for i=-(n-1):n
        p = pinduse + i + mmax + 1;
        q = qinduse + i + mmax + 1;
        aa = sum(gcoef(p))*2*pi;
        fq = fcoef(q).*qindinv_rep.';
        bb = sum(fq)*2*pi;

        qq = qqind + i + mmax + 1;
        fqq = fcoef(qq).*qqindinv.';
        bbq = sum(fqq)*2*pi;



        qdiffq(i+n) = aa - bb + bbq;

        p = i + mmax + 1;
        aa = gcoef(p)*2*pi;
        q = qgind + i + mmax + 1; 
        fq = fcoef(q).*qgindinv.';
        bb = sum(fq)*2*pi;
        gammas_est(i+n) = aa - bb+dlp;
    end
    qlnq = qdiffq + gammas_est;


end