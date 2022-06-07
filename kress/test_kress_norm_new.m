   
   nppw = 5
   rnorms = zeros(1,nppw);
   rnorms_inv = zeros(1,nppw);
   
   rnorms_unsc = zeros(1,nppw);
   rnorms_inv_unsc = zeros(1,nppw);
   ppws = (1:nppw)*5;
   
   
   
   

   for ippw=1:nppw
   %if(nargin <=3) 
       ppw = ppws(ippw);
   %end
   %if(nargin <= 2)
       bd_params.dir = 0;
   %end
   %if(nargin <= 1)
       kh = 40.0;
   %end
   %if(nargin <= 0)
       gpars.igeomtype = 1;
   %end
   
   
   
   %if(nargin == 5)
       n = nuse;
   %else
   n0 = 300;
   [src,~,~,~] = get_geom(gpars,n0);
   ds = sqrt(src(3,:).^2  +src(4,:).^2);
   [~,n] = size(src);
   L = sum(ds)*2*pi/n;
   nmin = 50;
   n  = max(nmin,ceil(kh*ppw*L/2/pi));
   %end
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
   
   fprintf('error in analytic solution test=%d\n',err_ex);
   
   D = diag(sqrt(wts(:)));
   Dinv = diag(1./sqrt(wts(:)));
   rnorms(ippw) = norm(D*mat*Dinv);
   rnorms_inv(ippw) = norm(D*mat_inv*Dinv);
   rnorms_unsc(ippw) = norm(mat);
   rnorms_inv_unsc(ippw) = norm(mat_inv);
   end
   
   m = 1:100;
   hm = besselh(m,1,kh);
   jmp = (besselj(m-1,kh) - besselj(m+1,kh))/2;
   jm = besselj(m,kh);
   evals = pi*hm.*(1j*kh*jmp + kh*jm);
   

   figure(1)
   clf()
  
   plot(ppws,rnorms,'k.','MarkerSize',20); hold on;
   plot(ppws,rnorms_unsc,'k^','MarkerSize',20);
   plot(m,abs(evals),'ks','MarkerSize',20);
   
   figure(2)
   clf()
   plot(ppws,rnorms_inv,'b.','MarkerSize',20); hold on;
   plot(ppws,rnorms_inv_unsc,'b^','MarkerSize',20);
   
   