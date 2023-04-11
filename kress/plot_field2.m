kh = 60;
gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;

ifsplit = [false true true];
rfacs = [16 16 32];

dirs = [0 0 pi/3 0 pi/3];

gpars = gpars0;
gpars.igeomtype = 5;
nuse = 1000;
bd_params = [];
bd_params.dir = -pi/2 + 0.3;

[ref_sols] = get_ref_sols(gpars,kh,bd_params,20,nuse);
zk = kh;
eps = 1e-6;
B = ref_sols;
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
caxis([0,3])
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,join(['./results/fields_ig' num2str(ig) '_kh' num2str(zk) '_aug11_2022.pdf']));




   