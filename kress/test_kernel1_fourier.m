addpath './src/'
gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;

gpars0.igeomtype = 2;

kh = 20;
ppw = 10;

% Estimate number of points

n0 = 300;
[src,~,~,~] = get_geom(gpars0,n0);
ds = sqrt(src(3,:).^2  +src(4,:).^2);
[~,n] = size(src);
L = sum(ds)*2*pi/n;
nmin = 50;
n  = max(nmin,ceil(kh*ppw*L/2/pi));
[src,ts,intpt,~] = get_geom(gpars0,n);


[kernel_1] = get_kernel1(kh,src);




figure(1)
clf
imagesc(abs(kernel_1))


[tt2,tt1]= meshgrid(ts);

tt1 = tt1(:);
tt2 = tt2(:);
kk1 = kernel_1(:);
%m = ceil(20*kh);
m = 2*n;
isign = 1;
eps = 1e-15;
f = finufft2d1(tt1,tt2,kk1,isign,eps,m,m)/4/n^2/(2*pi)^2;

figure(3)
clf
imagesc(ff)



s = 1.5;
vec1 = -n:n-1;
vec1 = vec1(:)/kh;
vec1 = abs(vec1).^s;


vval = (vec1'*abs(f)*vec1)/kh/sqrt(2*pi);


