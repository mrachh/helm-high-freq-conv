%% Test Neumann problem with data focussed on modes close to k

nterms = 5;
nn = 2*nterms+1;
coefs = (rand(nn,1)-0.5) + 1j*(rand(nn,1));


nk = 30;
k0 = 10;
dk = 5;
nppw = 5;



khs = k0:dk:k0+(nk-1)*dk;
err_sol = zeros(size(khs));
err_sol2 = zeros(size(khs));

spars = [];
spars.ifsplit = true;
spars.rfac = 8;



for iii=1:length(khs)
    
    kh = khs(iii);
    eta = kh;
        
    gpars = [];
    gpars.igeomtype = 6;
    gpars.xscale = a;
    gpars.yscale = b;
    
    nhalf = ceil(nppw*khs(iii));
    n = 2*nhalf;
    [src,ts,intpt,~] = get_geom(gpars,nhalf);
    
    nn = ceil(kh) - nterms:ceil(kh)+nterms;
    zexp = exp(1j*ts.'*nn);
    bn = coefs(:);
    dudnin = zexp*bn;
    
    nn2 = ceil(kh);
    dudnin2 = exp(1j*ts.'*nn2);
    [mat,Sik] = get_neu_mat_adj(kh,src,ts,eta,spars);
    S  = slmat(kh,src,ts,spars);
    D  = dlmat(kh,src,ts,spars);
    
    rhsuse = -dudnin;
    sol = mat\rhsuse;
    uscat = ((D + eye(n)/2)*Sik - 1j*eta*S)*sol;
    
    rhsuse = -dudnin2;
    sol2 = mat\rhsuse;
    uscat2 = ((D + eye(n)/2)*Sik - 1j*eta*S)*sol2;
    
    
    hn = besselh(nn,1,kh).';
    hnp = (besselh(nn-1,1,kh).' - besselh(nn+1,1,kh).')/2;

    cn = -bn.*hn/kh./hnp;
    uscat_series = zexp*cn;
    err_sol(iii) = norm(uscat-uscat_series)/norm(uscat_series);
    
    
    hn2 = besselh(nn2,1,kh).';
    hnp2 = (besselh(nn2-1,1,kh).' - besselh(nn2+1,1,kh).')/2;

    cn2 = -hn2/kh./hnp2;
    uscat_series2 = dudnin2*cn2;
    err_sol2(iii) = norm(uscat2-uscat_series2)/norm(uscat_series2);
    
    fprintf('kh=%d   err_sol = %d    err_sol2=%d\n',kh,err_sol(iii),err_sol2(iii));
end
figure(1)
clf
semilogy(khs,err_sol,'k.'); hold on; semilogy(khs,err_sol2,'r.');
xlabel('k','FontSize',15);
ylabel('\epsilon','FontSize',15);
saveas(gcf,'./results/neu_bdry_data.pdf');

save('neu_sol_data_jul12022.mat','khs','coefs','err_sol','err_sol2');
