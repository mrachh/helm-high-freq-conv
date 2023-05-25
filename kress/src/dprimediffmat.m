function T = dprimediffmat(kh,src,t,spars)%N,L)
% function to generate the matrix for the derivative of the double layer potential
% D'_{k} - D'_{ik}
% Input:
% N -> number of boundary points, should be even
% L -> arclegnth
% src -> domain
% src(1,N) = x(t)
% src(2,N) = y(t)
% src(3,N) = x'(t)
% src(4,N) = y'(t)
% src(5,N) = x''(t)
% src(6,N) = y''(t)
% Output:
% T is the normal derivative of the double layer potential

if(nargin == 3) 
    spars = [];
    spars.ifsplit = false;
    spars.rfac = 8;
end

kh2 = 1j*kh;
[~,kernel_k1_kh1,kernel_k2_kh1] = dprimelmat(kh,src,t,spars);
[~,kernel_k1_kh2,kernel_k2_kh2] = dprimelmat(kh2,src,t,spars);


x  = src(1,:);
y  = src(2,:);
dx  = src(3,:);
dy  = src(4,:);
%dx2 = src(5,:);
%dy2 = src(6,:);
%dx3 = src(7,:);
%dy3 = src(8,:);

N  = length(dx);
n  = N/2;
C  = 0.57721566490153286060651209008240243104215933593992; %euler-mascheroni's constant


drr = dx.^2+dy.^2;

quad_r = layer_quad(n);


T = 0.5 ./ repmat(sqrt(drr)',1,N) .* (gallery('circul',quad_r) .* (kernel_k1_kh1-kernel_k1_kh2) + ...
    pi/n * (kernel_k2_kh1-kernel_k2_kh2));


% S1= 0.5 ./ repmat(sqrt(drr)',1,N) .* (gallery('circul',quad_r).* kernel_m1 +...
%     pi/n * kernel_m2) .* ( bsxfun(@times,dx',dx)  + bsxfun(@times,dy',dy) );

