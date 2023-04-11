A = load('data_all_aug11_2022_arclengthparam.mat');

ik = 4;
ig = 5;


zk = A.kh(ik);
eps = 1e-6;
B = A.ref_sols_all{ik,ig};
[~,nn] = size(B.src);
ds = sqrt(B.src(3,:).^2 + B.src(4,:).^2);
wts = ds*2*pi/nn;

srcinfo = [];
srcinfo.sources = B.src(1:2,:);

srcinfo.dipvec = zeros(size(srcinfo.sources));
srcinfo.dipvec(1,:) = B.src(4,:)./ds;
srcinfo.dipvec(2,:) = -B.src(3,:)./ds;
sigma_use = B.sol_pw.*wts(:);
sigma_use = sigma_use.';
srcinfo.charges = -1j*zk*sigma_use;
srcinfo.dipstr = sigma_use;
srcinfo.nd = 1;

xvals = -4:0.003:4;
xt = -3:0.003:3;
[xx,yy]=meshgrid(xvals,xt);
targs = zeros(2,length(xx(:)));
targs(1,:) = xx(:);
targs(2,:) = yy(:);
pg = 0;
pgt = 1;
U = hfmm2d(eps,zk,srcinfo,pg,targs,pgt);

Uinc = exp(1i*zk*(targs(1,:)*cos(B.dir) + targs(2,:)*sin(B.dir)));
Uplot = Uinc - U.pottarg;
Uplot = reshape(Uplot,size(xx));
xtarg = reshape(targs(1,:),size(xx));
ytarg = reshape(targs(2,:),size(yy));
figure(1)
clf
h = pcolor(xtarg,ytarg,abs(Uplot)); hold on;
set(h,'EdgeColor','none')
fill(B.src(1,:),B.src(2,:),'w');
axis equal
xlim([-4,4])
ylim([-3,3])
cb = colorbar()
ax = gca;
caxis([0,5])
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,join(['./results/fields_ig' num2str(ig) '_kh' num2str(zk) '_aug11_2022.pdf']));




   