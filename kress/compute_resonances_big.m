clear;
clc;
addpath './src'
gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;


gpars0.igeomtype = 3;
n = 300;
[src2,ts,intpt,~] = get_geom(gpars0,n);
ds = sqrt(src2(3,:).^2  +src2(4,:).^2);
[~,n] = size(src2);
L0 = sum(ds)*2*pi/n;




kh0 = 30;
panlen = 5;
npan = 9;
ppw = 5;
ncheb = 1024;

detchebs = zeros(ncheb,npan);
detchebs2 = zeros(ncheb,npan);
khvec = zeros(ncheb,npan);

khi = 0;

for ipan=1:npan
    fprintf('Starting panel=%d\n',ipan);
    kh_start = kh0 + (ipan-1)*panlen;
    kh_end = kh0 + ipan*panlen;
    
    khs = (chebpts(ncheb)+1)/2*panlen + kh_start;
    khvec(:,ipan) = khs;
    for i=1:ncheb
        n = max(ceil(ppw*khs(i)*L0/2/pi),300);
        [src2,ts,intpt,~] = get_geom(gpars0,n);
        ds = sqrt(src2(3,:).^2  +src2(4,:).^2);
        [~,n] = size(src2);
        wts = sqrt(src2(3,:).^2 + src2(4,:).^2)*2*pi/n;
 
   
        D = diag(sqrt(wts(:)));
        Dinv = diag(1./sqrt(wts(:)));  
        khuse = khs(i) + 1j*khi;
        Ainv = inv(get_kress_mat(khuse,src2,ts));
        detchebs(i,ipan) = norm(Ainv);
        detchebs2(i,ipan) = norm(D*Ainv*Dinv);
    end
end
save(['norm_inv_cavity_data_jun13_2022_ima' num2str(khi) '_ncheb' int2str(ncheb) '.mat'],'L0','detchebs','detchebs2',...
    'gpars0','kh0','khvec','ncheb','npan','panlen','ppw');