addpath './src';
addpath '~/git/inverse-obstacle-scattering2d/src/';
kh = 42.5;
gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;

ifsplit = [false true true];
rfacs = [16 16 32];

dirs = [0 0 -pi/2+0.2 0 -pi/2+0.2];

nuse = 1000;
ig = 3;
gpars = gpars0;
gpars.igeomtype = 3;
bd_params = [];
bd_params.dir = dirs(ig);

[ref_sols] = get_ref_sols(gpars,kh,bd_params,10,nuse);

xs = ref_sols.src(1,:);
xx = xs(:);
ys = ref_sols.src(2,:);
yy = ys(:);
    
th = dirs(5);
x0 = -1;
y0 = 0.1;
s0 = 1.0/sqrt(kh);
rr = sqrt((xs-x0).^2+(ys-y0).^2);
rr = rr(:);

B = ref_sols;
ttt = B.ts(:)*sqrt(B.src(3,1).^2 + B.src(4,1).^2);
ttt = B.ts(:);
ntt = 1200;
ttuse = ttt(ntt);
ttuse = 3.91322070806951;

u_wp = exp(-(ttt-ttuse).^2/s0.^2).*exp(1j*kh*cos(0.34*pi)*ttt);
sol_wp = 2*ref_sols.mat\u_wp;

%sol_wp = ref_sols.sol_pw;


zk = kh;

eps = 1e-6;

[~,nn] = size(B.src);
ds = sqrt(B.src(3,:).^2 + B.src(4,:).^2);
wts = ds*2*pi/nn;

srcinfo = [];
srcinfo.sources = B.src(1:2,:);

srcinfo.dipvec = zeros(size(srcinfo.sources));
srcinfo.dipvec(1,:) = B.src(4,:)./ds;
srcinfo.dipvec(2,:) = -B.src(3,:)./ds;
sigma_use = sol_wp.*wts(:);
sigma_use = sigma_use.';
srcinfo.charges = -1j*zk*sigma_use;
srcinfo.dipstr = sigma_use;
srcinfo.nd = 1;

xvals = -4:0.002:4;
xt = -3:0.002:3;
[xx,yy]=meshgrid(xvals,xt);
targs = zeros(2,length(xx(:)));
targs(1,:) = xx(:);
targs(2,:) = yy(:);
pg = 0;
pgt = 1;
U = hfmm2d(eps,zk,srcinfo,pg,targs,pgt);

xt = xx(:);
yt = yy(:);
rr = sqrt((xt-x0).^2 + (yt-y0).^2);
Uinc = exp(-rr.^2/s0.^2).*exp(1j*kh*(xt*cos(th) + yt*sin(th)));
Uinc = exp(1i*zk*(xx*cos(th) + yy*sin(th)));
Uinc = Uinc;
Uplot = 0*Uinc(:).' - U.pottarg;
Uplot = reshape(Uplot,size(xx));
xtarg = reshape(targs(1,:),size(xx));
ytarg = reshape(targs(2,:),size(yy));
figure(1)
clf
h = pcolor(xtarg,ytarg,abs(Uplot)); hold on;
set(h,'EdgeColor','none')
fill(B.src(1,:),B.src(2,:),'w'); hold on;
plot(B.src(1,ntt),B.src(2,ntt),'r.','MarkerSize',10);
axis equal
xlim([-4,4])
ylim([-3,3])
cb = colorbar();
caxis([0,5]);
ax = gca;
fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
