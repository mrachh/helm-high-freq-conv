clear;
clc;

addpath('./src')
addpath('~/git/helm_transmission_2d/matlab/src');

%% Generate geometry and right hand side for impendence problem


nk = 30;
k0 = 10;
dk = 5;
nppw = 3;



khs = k0:dk:k0+(nk-1)*dk;
err_sol = zeros(size(khs));
err_data = zeros(size(khs));
err_diff = zeros(size(khs));






for iii = 1:nk


    src = [0.01;-0.07];
    targ = [12.1;5.2];

    % Generate geometry
    a = 1.0; b=1.0;


    spars = [];
    spars.ifsplit = true;
    spars.rfac = 16;

    
    gpars = [];
    gpars.igeomtype = 6;
    gpars.xscale = a;
    gpars.yscale = b;
    
    nhalf = ceil(nppw*khs(iii));
    n = 2*nhalf;
    
    [srcinfo,h] = ellipse_new(a,b,n);


    zk = khs(iii);
    uex = helm_c_p(zk,src,targ);
    uin = helm_c_p(zk,src,srcinfo);
    dudnin = helm_c_gn(zk,src,srcinfo);

    xx = srcinfo(1,:).';
    yy = srcinfo(2,:).';

    rnx = srcinfo(3,:).';
    rny = srcinfo(4,:).';

    dir = 0;
    u_pw = exp(1i*zk*(xx*cos(dir) + yy*sin(dir)));
    dudn_pw = 1i*zk*(rnx*cos(dir) + rny*sin(dir)).*u_pw;


    alpha = zk;
    zk2 = 1j*abs(zk);
    norder = 16;
    sikmat = slp_mat(zk2,norder,h,srcinfo);
    smat = slp_mat(zk,norder,h,srcinfo);
    dmat = dlp_ext_mat(zk,norder,h,srcinfo);

    %% Neumann test

    eta = zk;
    alpha = zk;
    zpars = complex(zeros(3,1));
    zpars(1) = zk;
    zpars(2) = alpha;
    zpars(3) = eta;

    rhs = dudnin;

    xmat = rpcomb_neu_ext_mat(zpars,norder,h,srcinfo);
    soln = xmat\rhs;
    soln2 = sikmat*soln;

    soln_pw = -xmat\dudn_pw;
    soln2_pw = sikmat*soln_pw;

    % Test resulting solution
    zs = helm_c_p(zk,srcinfo,targ);
    zd = helm_d_p(zk,srcinfo,targ);
    ucomp = (zs.*srcinfo(5,:))*soln*h + 1j*alpha*(zd.*srcinfo(5,:))*soln2*h;

    err = norm(ucomp(:) - uex(:));
    fprintf("Error in Neumann test = %d\n\n\n",err);

    % Test dirichlet data on the boundary
    ubd = 1j*alpha*dmat*soln2 + smat*soln;
    err = norm(ubd(:)-uin(:));
    fprintf("Error in potential on the boundary = %d \n\n\n",err);


    uscat_pw = 1j*alpha*dmat*soln2_pw + smat*soln_pw;
    % Now test neumann solution without split and using the Neumann formulation
    % in the paper
    
    [src,ts,intpt,~] = get_geom(gpars,nhalf);

%     figure(1)
%     plot(src(1,:),src(2,:),'k.'); hold on; plot(srcinfo(1,:),srcinfo(2,:),'r.');
    S_kress  = slmat(zk,src,ts);
    D_kress  = dlmat(zk,src,ts);
    Sik_kress = slmat(zk2,src,ts,spars);
    T_kress = dprimelmat(zk,src,ts);

    Bk = 1j*eta*(eye(n)/2-D_kress) + Sik_kress*T_kress;

    rhs_kress = 1j*eta*u_pw - Sik_kress*dudn_pw;
    sol_kress = Bk\rhs_kress;
    utot = 2*(u_pw + D_kress*sol_kress);
    uscat = utot-u_pw;

    err_diff(iii) = norm(uscat-uscat_pw);
    fprintf('error between two modes of computing uscat=%d\n',err_diff(iii))

    % Check jacobi anger expansion of Neumann data

    nn = -300:1:300;
    jnm1 = besselj(nn-1,zk).';
    jnp1 = besselj(nn+1,zk).';
    imas = (1j).^(nn).';
    zexp = exp(1j*ts.'*nn);


    bn =  ((jnm1 - jnp1).*imas)*zk/2;

    dudn_series = zexp*bn;

    % Test accuracy with analytic solution
    hn = besselh(nn,1,zk).';
    hnp = (besselh(nn-1,1,zk).' - besselh(nn+1,1,zk).')/2;

    cn = -bn.*hn/zk./hnp;
    uscat_series = zexp*cn;

    err_data(iii) = norm(dudn_series-dudn_pw);
    err_sol(iii) = norm(uscat_series-uscat);
    fprintf('error in data=%d\n',err_data(iii));
    fprintf('error in sol=%d\n',err_sol(iii));
end

