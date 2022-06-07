khs = 20:5:60;
nkh = length(khs);

rnorms = zeros(1,nkh);
rnorms_err = zeros(1,nkh);


gpars.igeomtype = 1;

ppw = 6;

for i=1:nkh
    kh = khs(i);
    bd_params = [];
    bd_params.dir = 0;
    nuse = 1000;
    [ref_sols] = get_ref_sols(gpars,kh,bd_params,20,nuse);

    n = ceil(kh*ppw);
    m1 = nuse;
    qn = get_qmat(n);
    qm1 = get_qmat(m1);

    enm1 = get_emat(m1,n);
    em1n = get_emat(n,m1);
    qnkqn_m1 = enm1*qm1*ref_sols.mat*em1n*qn;
    rn1 = norm(qnkqn_m1);


    m_ind = 1:n;
    js = besselj(m_ind,kh);
    jders = 0.5*(besselj(m_ind-1,kh) - besselj(m_ind+1,kh));
    hs = besselh(m_ind,1,kh);

    gammas = kh*(jders - 1j*js).*hs*pi;
    gmax = max(abs(gammas));

    err1 = norm(rn1-gmax);
    
    rnorms(i) = gmax;
    rnorms_err(i) = err1;
    fprintf('difference between estimated norms = %d\n',err1);
end


figure(1)
clf
plot(khs,rnorms,'k.','MarkerSize',20);
xlabel('k','FontSize',15);
ylabel('|Q_{N}(I+ L)Q_{N}|','FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

saveas(gcf,'./results/qnlqn_norm.pdf');


figure(2)
clf
semilogy(khs,rnorms_err,'k.','MarkerSize',20);
xlabel('k','FontSize',15);
ylabel('\epsilon_{L}','FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

saveas(gcf,'./results/qnlqn_norm_err.pdf');