function [out] = dens_interp_err(ref_sols,gpars,kh,ppw,spars)
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
%    kh: helmholtz wavenumber
%    ppw: integer
%      min points per wavelength
%    ref_sols: reference solution struct
%        res_sols.src - src struct used to generate reference solution
%        ref_sols.ts - t values of chunkie parameterization
%        ref_sols.dir = bd_params.dir
%        ref_sols.sol_pw - density corresponding to plane wave data
%    spars: splitting parameters struct
%        spars.ifsplit: boolean
%           if true, then use partition of unity given by rfac
%        spars.rfac: double
%           partition of unity given by 
%
%  Output arguments:
%    out = output structure
%      out.err_ex: error in computing analytic solution on given grid
%      out.err_dens: error in interpolated density
%      out.rnorm: norm of matrix
%      out.rnorm_inv = norm of inverse of matrix
%        
%    

   addpath ~/git/finufft/matlab
   if(nargin <=4) 
       spars = [];
       spars.ifsplit = false;
       spars.rfac = 16;
   end
   if(nargin <=3) 
       ppw = 5.0;
   end
   if(nargin <= 2)
       kh = 20.0;
   end
   if(nargin <= 1)
       gpars.igeomtype = 1;
   end
   
   
   ds = sqrt(ref_sols.src(3,:).^2 + ref_sols.src(4,:).^2);
   [~,n] = size(ref_sols.src);
   L = sum(ds)*2*pi/n;
   nmin = 50;
   n  = max(nmin,ceil(kh*ppw*L/2/pi));
   
   [src,ts,intpt,~] = get_geom(gpars,n);
   
   
   mat = get_kress_mat(kh,src,ts,spars);
   
   [~,nn] = size(src);
   wts = sqrt(src(3,:).^2 + src(4,:).^2)*2*pi/nn;
   
   xs = src(1,:);
   xx = xs(:);
   ys = src(2,:);
   yy = ys(:);
    
   
   dir = ref_sols.dir;
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
   
   ms = length(sol_pw);
   eps = 1e-15;
   isign = -1;
   sol_pw_hat = finufft1d1(ts,sol_pw,isign,eps,ms)/ms;
   
   
   
   isign = 1;
   sol_pw_interp = finufft1d2(ref_sols.ts,isign,eps,sol_pw_hat);
   
  
   [~,nref] = size(ref_sols.src);
   wref = sqrt(ref_sols.src(3,:).^2 + ref_sols.src(4,:).^2)*2*pi/nref;
   
   err1 = norm((sol_pw_interp(:)-ref_sols.sol_pw(:)).*sqrt(wref(:)));
   r1 = norm(ref_sols.sol_pw(:).*sqrt(wref(:)));
   
   err_dens = err1/r1;
   
   
   D = diag(sqrt(wts(:)));
   Dinv = diag(1./sqrt(wts(:)));
   rnorm = norm(D*mat*Dinv);
   rnorm_inv = norm(D*mat_inv*Dinv);
   
   out = [];
   out.err_ex = err_ex;
   out.err_dens = err_dens;
   out.rnorm = rnorm;
   out.rnorm_inv = rnorm_inv;   
   fprintf('error in exact solution test=%d\n',out.err_ex);
   fprintf('error in pw solution test=%d\n',out.err_dens);
   
   m1 = nref/2;
   qn = get_qmat(n);
   qm1 = get_qmat(m1);
    
   enm1 = get_emat(m1,n);
   em1n = get_emat(n,m1);
   qnknqn = qn*mat;
   qnkqn_m1 = qn*enm1*qm1*ref_sols.mat*em1n*qn;
   
   err_qnkqn = norm((qnknqn-qnkqn_m1)*Dinv);
   err_qnkqn2 = norm(D*(enm1*qm1*ref_sols.mat*em1n*qn - mat)*Dinv);
   err_qnkqn3 = norm(enm1*qm1*ref_sols.mat*em1n*qn - mat);
   out.qnknqn_norm = norm(qnknqn);
   out.qnkqn_norm = norm(qnkqn_m1);
   out.err_qnkqn = err_qnkqn;
   out.err_qnkqn2 = err_qnkqn2;
   out.err_qnkqn3 = err_qnkqn3;
   fprintf('error in qnkqn =%d\n',out.err_qnkqn);

end
