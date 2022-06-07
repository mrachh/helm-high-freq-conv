kh_res = 40.43440306421;
kh_nonres = 40.247577830244;


addpath './src/'
gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;

gpars0.igeomtype = 5;

ppw = [2.0 2.2 2.5 3.0 4.0 6.0 10.0 12.0 15.0 20.0];
nppw = length(ppw);


n0 = 300;
[src,~,~,~] = get_geom(gpars0,n0);
ds = sqrt(src(3,:).^2  +src(4,:).^2);
[~,n] = size(src);
L = sum(ds)*2*pi/n;
nmin = 50;
n  = max(nmin,ceil(kh_res*ppw*L/2/pi));



n0_four = n(1);
n1_four = n(end-2);
n2_four = n(end-1);
kdecay = 5;


% set up fourier data
tol = 1e-15;

beta = 1/(n2_four-n1_four)*(log(1/tol) + kdecay*log(n0_four/n1_four));

kvec1 = -n0_four:n0_four;
fvec1 = ones(size(kvec1));
kvec2 = -n1_four:(-n0_four-1);
fvec2 = abs(n0_four./kvec2).^kdecay;
kvec3 = (n0_four+1):n1_four;
fvec3 = abs(n0_four./kvec3).^kdecay;
kvec4 = -n2_four:(-n1_four-1);
fvec4 = (n0_four/n1_four)^kdecay*exp(-beta*(abs(kvec4)-n1_four));
kvec5 = (n1_four+1):n2_four;
fvec5 = (n0_four/n1_four)^kdecay*exp(-beta*(abs(kvec5)-n1_four));

kvec = [kvec1 kvec2 kvec3 kvec4 kvec5];
fvec = [fvec1 fvec2 fvec3 fvec4 fvec5];


src_all = cell(nppw,1);
mat_res = cell(nppw,1);
mat_nonres = cell(nppw,1);
ts_all = cell(nppw,1);
bd_data_all = cell(nppw,1);
sigma_res = cell(nppw,1);
sigma_nonres = cell(nppw,1);


for ippw = 1:nppw
    [src_all{ippw},ts_all{ippw},intpt,~] = get_geom(gpars0,n(ippw));
    mat_res{ippw} = get_kress_mat(kh_res,src_all{ippw},ts_all{ippw});
    mat_nonres{ippw} = get_kress_mat(kh_nonres,src_all{ippw},ts_all{ippw});
    isign = -1;
    tol2 = 1e-15;
    bd_data_all{ippw} = real(finufft1d3(kvec,fvec,isign,tol2,ts_all{ippw}))/length(kvec);
    sigma_res{ippw} = mat_res{ippw}\bd_data_all{ippw};
    sigma_nonres{ippw} = mat_nonres{ippw}\bd_data_all{ippw};    
end


err_res = zeros(nppw-1,1);
err_nonres = zeros(nppw-1,1);


ref_ts = ts_all{nppw};
sol_ref_res = sigma_res{nppw};
sol_ref_nonres = sigma_nonres{nppw};

[~,nref] = size(src_all{nppw});
wref = sqrt(src_all{nppw}(3,:).^2 + src_all{nppw}(4,:).^2)*2*pi/nref;
for ippw=1:nppw-1
   sol_res = sigma_res{ippw};
   sol_nonres = sigma_nonres{ippw};
   ts = ts_all{ippw};

   ms = length(sol_res);
   eps = 1e-15;
   isign = -1;
   sol_res_hat = finufft1d1(ts,sol_res,isign,eps,ms)/ms;
   sol_nonres_hat = finufft1d1(ts,sol_nonres,isign,eps,ms)/ms;
   
   
   
   isign = 1;
   sol_res_interp = finufft1d2(ref_ts,isign,eps,sol_res_hat);
   sol_nonres_interp = finufft1d2(ref_ts,isign,eps,sol_nonres_hat);
   
  
   err1 = norm((sol_res_interp(:)-sol_ref_res(:)).*sqrt(wref(:)));
   r1 = norm(sol_ref_res(:).*sqrt(wref(:)));
   
   err_res(ippw) = err1/r1;
   
   
   
   err1 = norm((sol_nonres_interp(:)-sol_ref_nonres(:)).*sqrt(wref(:)));
   r1 = norm(sol_ref_nonres(:).*sqrt(wref(:)));
   
   err_nonres(ippw) = err1/r1;
end

