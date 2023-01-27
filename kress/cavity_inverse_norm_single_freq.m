% compute norm of inverse for cavity at a single frequency
% kh = 37.212;
% kh = 37.212733389;
kh = 64.5386641434;
% kh = 30.635;

addpath './src/'
gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;

gpars0.igeomtype = 1;
gpars0.igeomtype = 3;
n = 900;
[src2,ts,intpt,~] = get_geom(gpars0,n);
ds = sqrt(src2(3,:).^2  +src2(4,:).^2);
[~,n] = size(src2);
L = sum(ds)*2*pi/n;
wts = sqrt(src2(3,:).^2 + src2(4,:).^2)*2*pi/n;
 
   
D = diag(sqrt(wts(:)));
Dinv = diag(1./sqrt(wts(:)));
Ainv = inv(get_kress_mat(kh,src2,ts));
norm1 = norm(Ainv);
norm2 = norm(D*Ainv*Dinv);

fprintf('norm1=%d\n',norm1);
fprintf('norm2=%d\n',norm2);