function [ref_sols] = get_ref_sols(gpars,kh,bd_params,ppw,nuse)
% This subroutine computes the reference solutions 
% for plane wave incident data, and singular data
% given by r^(2*bd_params.spow+1), where spow was an integer
% 
%  Input arguments:
%    gpars.igeomtype: geometry type
%        igeomtype = 1, unit circle
%        igeomtype = 2, star shaped geometries (requires gpars.nosc, and
%          gpars.rfac) r(t) = 1 + rfac*cos(nosc*t)
%        igeomtype = 3, cavity (requires gpars.a (width), gpars.b
%        (opening angle), gpars.n0 (min number of fourier modes to use in
%         constructing radius function
%        igeomtype = 4, star shaped geometries, with arc-length
%          parametrization
%        igeomtype = 5, cavity geometry with arc-length parametrization
%        igeomtype = 6, ellipse, (requires gpars.xscale, gpars.yscale,
%          defining the axis length along the x-axis and y-axis)
%        igeomtype = 7, ellipse, with arc-length parametrization
%        igeomtype = 8, ice-cream cone geometry (requires gpars.alpha,
%           gpars.beta, gpars.gamma)
%        igeomtype = 9, arc-length parametrization of ice-cream cone
%           geometry)
%    kh: helmholtz wavenumber
%    bd_params: parameters for defining the boundary data
%         bd_params.dir - incident direction
%    ppw: integer
%      min points per wavelength
%
%  Output arguments:
%     ref_sols: reference solution struct
%        res_sols.src - src struct used to generate reference solution
%        ref_sols.ts - t values of chunkie parameterization
%        ref_sols.dir = bd_params.dir
%        ref_sols.sol_pw - density corresponding to plane wave data
%        
%    


   if(nargin <=3) 
       ppw = 20;
   end
   if(nargin <= 2)
       bd_params.dir = 0;
   end
   if(nargin <= 1)
       kh = 20.0;
   end
   if(nargin <= 0)
       gpars.igeomtype = 1;
   end
   
   
   
   if(nargin == 5)
       n = nuse;
   else
       n0 = 300;
      [src,~,~,~] = get_geom(gpars,n0);
      ds = sqrt(src(3,:).^2  +src(4,:).^2);
      [~,n] = size(src);
      L = sum(ds)*2*pi/n;
      nmin = 50;
      n  = max(nmin,ceil(kh*ppw*L/2/pi));
   end
   fprintf('n=%d\n',n);
   [src,ts,intpt,~] = get_geom(gpars,n);
   
   
   
   if(~isfield(bd_params,'dir'))
       bd_params.dir = 0;
   end
   
   mat = get_kress_mat(kh,src,ts);
   
   [~,nn] = size(src);
   wts = sqrt(src(3,:).^2 + src(4,:).^2)*2*pi/nn;
   
   ref_sols = [];
   ref_sols.src = src;
   ref_sols.ts = ts;
   ref_sols.dir = bd_params.dir;
   
   xs = src(1,:);
   xx = xs(:);
   ys = src(2,:);
   yy = ys(:);
    
   
   dir = bd_params.dir;
   u_pw = exp(1i*kh*(xx*cos(dir) + yy*sin(dir)));
   source = intpt;
   rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
   rr = rr(:);
   u_test = (1i/4)*besselh(0,1,kh*rr); 
    
   mat_inv = inv(mat);
   sol_test = mat_inv*2*u_test;
   sol_pw = mat_inv*2*u_pw;
    
   
   tgt = [3.1 2.8];
   tgt = tgt(:);
   eta = -kh;
   pot_test = (dlmat_out(kh,src,tgt)+1i*eta*slmat_out(kh,src,tgt))*sol_test; 
   
   rr = sqrt((tgt(1)-source(1)).^2+(tgt(2)-source(2)).^2);
   rr = rr(:);
   pot_ex = (1i/4)*besselh(0,1,kh*rr);
   err_ex = norm(pot_ex-pot_test)/norm(sol_test.*sqrt(wts(:)));
   ref_sols.err_ex = err_ex;
   ref_sols.sol_pw = sol_pw;
   ref_sols.mat = mat;

   
end
