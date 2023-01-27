% verify circle interior resonance

addpath './src/'
gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;

gpars0.igeomtype = 1;









kh = fzero(@(z) besselj(0,z),2);

ppw = 10;

% Estimate number of points

n0 = 300;
[src,ts,~,~] = get_geom(gpars0,n0);
ds = sqrt(src(3,:).^2  +src(4,:).^2);
[~,n] = size(src);
L = sum(ds)*2*pi/n;
nmin = 50;

D  = dlmat(kh,src,ts);
mat = eye(n) - 2*D;
fprintf('Determinant of D on the circle = %d \n',det(mat))


gpars0.igeomtype = 3;
n = 900;
[src2,ts,intpt,~] = get_geom(gpars0,n);
ds = sqrt(src2(3,:).^2  +src2(4,:).^2);
[~,n] = size(src2);
L = sum(ds)*2*pi/n;
wts = sqrt(src2(3,:).^2 + src2(4,:).^2)*2*pi/n;
 
   
D = diag(sqrt(wts(:)));
Dinv = diag(1./sqrt(wts(:)));

eps = 1e-7;
p = chebfunpref; p.chebfuneps = eps;
p.splitting = 0; p.maxLength=257;

chebabs = [40,45];

spars = [];
spars.ifsplit = false;

tic,
ncheb = 64;
khvec = (chebpts(ncheb)+1)/2*(chebabs(2)-chebabs(1)) + chebabs(1);
detchebs = zeros(ncheb,1);
detchebs2 = zeros(ncheb,1);
for i=1:ncheb
    Ainv = inv(get_kress_mat(khvec(i),src2,ts));
    detchebs(i) = norm(Ainv);
    detchebs2(i) = norm(D*Ainv*Dinv);
end
toc;
figure 
clf
plot(khvec,detchebs,'r.'); hold on;
plot(khvec,detchebs2,'k.');


