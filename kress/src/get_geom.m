function [src,t,intpt,bdpt] = get_geom(gpars,n)
%  This subroutine constructs an equispaced discretization of
%  three types of geometries with 2n equispaced nodes
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
%    n: half number of discretization points
%  
%  Output arguments:
%    src: (8,2*n)
%      boundary info needed for generating kress matrices
%    t: (1,2*n)
%      t values along the curve
%
   
   src = zeros(8,2*n);
   t = 0:pi/n:2*pi-pi/n;
   
   xsc = 1;
   ysc = 1;
   if(gpars.igeomtype == 6 || gpars.igeomtype == 7)
       xsc = gpars.xscale;
       ysc = gpars.yscale;
   end
   if(gpars.igeomtype == 1 || gpars.igeomtype == 6)
     xs   = xsc*cos(t);
     dxs  = -xsc*sin(t);
     ddxs = -xsc*cos(t);
     d3xs = xsc*sin(t);
     
     ys   = ysc*sin(t);
     dys  = ysc*cos(t);
     ddys = -ysc*sin(t);
     d3ys = -ysc*cos(t);
     intpt = [0.01*xsc 0.3*ysc];  
     bdpt = [xsc,0];
       
   elseif(gpars.igeomtype == 2)
     rfac = gpars.rfac;
     nosc = gpars.nosc;
     rt = 1 + rfac*cos(nosc*t);
     drdt = -rfac*nosc*sin(nosc*t);
     d2rdt2 = -rfac*nosc*nosc*cos(nosc*t);
     d3rdt3 = rfac*nosc*nosc*nosc*sin(nosc*t);
    
     %boundary
     xs   = rt.*cos(t);
     ys   = rt.*sin(t);
     dxs  = -rt.*sin(t) + drdt.*cos(t);
     dys  = rt.*cos(t) + drdt.*sin(t);
     ddxs = -rt.*cos(t) - 2*drdt.*sin(t) + d2rdt2.*cos(t);
     ddys = -rt.*sin(t) + 2*drdt.*cos(t) + d2rdt2.*sin(t);
     d3xs = d3rdt3.*cos(t) - d2rdt2.*sin(t) - 2*(d2rdt2.*sin(t) + ... 
        drdt.*cos(t)) - (drdt.*cos(t) -rt.*sin(t));
     d3ys = d3rdt3.*sin(t) + d2rdt2.*cos(t) + 2*(d2rdt2.*cos(t) -drdt.*sin(t)) ...
           - (drdt.*sin(t) + rt.*cos(t));
     intpt = [0.01 0.3];
     bdpt = [xs(1), ys(1)];
     
   end
   
   if(gpars.igeomtype <=2 || gpars.igeomtype == 6)
       %setting domain
       src(1,:) = xs;
       src(2,:) = ys;
       src(3,:) = dxs;
       src(4,:) = dys;
       src(5,:) = ddxs;
       src(6,:) = ddys;
       src(7,:) = d3xs;
       src(8,:) = d3ys;
   end
   
   if(gpars.igeomtype == 3)
       a = gpars.a;
       b = gpars.b;
       n0 = gpars.n0;
       src = larrycup(a,b,2*n,n0);
       intpt = [0.01 -1.01];
       bdpt = [src(1,1), src(2,1)];
   end
   
   
   if(gpars.igeomtype == 4)
     rfac = gpars.rfac;
     nosc = gpars.nosc;
     
     coefs = zeros(2*nosc+1,1);
     coefs(1) = 1;
     coefs(nosc+1) = rfac;
     
     src_use = geometries.starn(coefs,nosc,2*n);
     intpt = [0.01 0.3];
     
     
     
     nh = 1;
     h_upd = zeros(2*nh+1,1);
     src_out = rla.update_geom(src_use,nh,h_upd);
   end

   
   if(gpars.igeomtype == 5)
       a = gpars.a;
       b = gpars.b;
       n0 = gpars.n0;
       src_use = geometries.larrycup(a,b,2*n,n0);
       intpt = [0.01 -1.01];
       
       nh = 1;
       h_upd = zeros(2*nh+1,1);
       src_out = rla.update_geom(src_use,nh,h_upd);
   end
   
   
   if(gpars.igeomtype == 7)
     src_use = geometries.ellipse(xsc,ysc,2*n);
       
     intpt = [0.01*xsc 0.3*ysc];
       
     nh = 1;
     h_upd = zeros(2*nh+1,1);
     src_out = rla.update_geom(src_use,nh,h_upd);
   end
   
   if(gpars.igeomtype == 8 || gpars.igeomtype == 9)
       
     alpha = gpars.alpha;
     beta = gpars.beta;
     gamma = gpars.gamma;
     src_out = geometries.cone(alpha,beta,gamma,2*n);
     intpt = [0.7 0.1];
     
     
     
     nh = 1;
     h_upd = zeros(2*nh+1,1);
     if(gpars.igeomtype == 9)
         src_out = rla.update_geom(src_out,nh,h_upd);
     end  
   end
   
   if(gpars.igeomtype == 4 || gpars.igeomtype == 5 || ...
       gpars.igeomtype == 7 || gpars.igeomtype == 8 || ...
       gpars.igeomtype == 9)
   
       xs = src_out.xs;
       ys = src_out.ys;
       Z = xs + 1j*ys;


       zhat = fft(Z(:))/n/2;
       t1 = (0:(2*n-1))/n/2;
       xy = fourierZ(zhat,t1);
       dxydt = fourierZp(zhat,t1);
       d2xydt = fourierZpp(zhat,t1);
       d3xydt = fourierZppp(zhat,t1);
       src(1,:) = real(xy);
       src(2,:) = imag(xy);
       src(3,:) = real(dxydt)/2/pi;
       src(4,:) = imag(dxydt)/2/pi;
       src(5,:) = real(d2xydt)/4/pi/pi;
       src(6,:) = imag(d2xydt)/4/pi/pi;
       src(7,:) = real(d3xydt)/8/pi/pi/pi;
       src(8,:) = imag(d3xydt)/8/pi/pi/pi;  
       bdpt = [src(1,1), src(2,1)];
   end
end