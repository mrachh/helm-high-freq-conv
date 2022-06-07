function [out] = get_fourier_ker1_hsnorm(gpars,kh,ppw,spars,spow)
% This subroutine computes the scaled hs norm of the kernel multiplying
% the log term. In particular if the kernel is
% f(t,s) = \sum_{mn} f_{mn} e^{imt} e^{ins}, then this 
% subroutine returns |f_{mn}| <m/k>^s <n/k>^s
% where s here is spow
% 
% 
%  Input arguments:
%    gpars.igeomtype: geometry type
%        igeomtype = 1, unit circle
%        igeomtype = 2, star shaped geometries (requires gpars.nosc, and
%          gpars.rfac) r(t) = 1 + rfac*cos(nosc*t)
%        igeomtype = 3, cavity (requires gpars.a (width), gpars.b
%        (opening angle), gpars.n0 (min number of fourier modes to use in
%         constructing radius function
%    kh: helmholtz wavenumber
%    ppw: integer
%      min points per wavelength
%    spars: splitting parameters struct
%        spars.ifsplit: boolean
%           if true, then use partition of unity given by rfac
%        spars.rfac: double
%           partition of unity given by 
%
%  Output arguments:
%    out = output structure
%      out.val: computes hs norm
%      out.fhat: fourier coefs, ordered in [-n,n-1]^2
%        
%    

   addpath ~/git/finufft/matlab
   if(nargin <5)
       spow = 2;
   end
   if(nargin <4) 
       spars = [];
       spars.ifsplit = false;
       spars.rfac = 16;
   end
   if(nargin <3) 
       ppw = 5.0;
   end
   if(nargin < 2)
       kh = 20.0;
   end
   if(nargin < 1)
       gpars.igeomtype = 1;
   end
   n0 = 300;
   [src,~,~,~] = get_geom(gpars,n0);
   ds = sqrt(src(3,:).^2  +src(4,:).^2);
   [~,n] = size(src);
   L = sum(ds)*2*pi/n;
   nmin = 50;
   n  = max(nmin,ceil(kh*ppw*L/2/pi));
   
   [src,ts,~,~] = get_geom(gpars,n);
   
   zpars = [-1j*kh 1];
   

   [kernel_1] = get_kernel1(kh,src,zpars,spars);


    [tt2,tt1]= meshgrid(ts);

    tt1 = tt1(:);
    tt2 = tt2(:);
    kk1 = kernel_1(:);
    m = 2*n;
    isign = 1;
    eps = 1e-15;
    f = finufft2d1(tt1,tt2,kk1,isign,eps,m,m)/4/n^2/(2*pi)^2;

    
    vec1 = -n:n-1;
    vec1 = vec1(:)/kh;
    vecr = abs(vec1).^(2*spow);
    vecl = abs(vec1).^(spow);


    out.val = (vecl'*sqrt(abs(f.^2)*vecr))/kh/sqrt(2*pi);
    out.fhat = f;
end
