function S = slmat(kh,src,t,spars)%N,L)
% function to generate the matrix for the single layer potential
% Input:
% N -> number of boundary points, should be even
% L -> arclegnth
% src -> domain
% src(1,N) = x(t)
% src(2,N) = y(t)
% src(3,N) = x'(t)
% src(4,N) = y'(t)
% Output:
% S is the single layer potential


if(nargin == 3) 
    spars = [];
    spars.ifsplit = false;
    spars.rfac = 16;
end
x  = src(1,:);
y  = src(2,:);
dx = src(3,:);
dy = src(4,:);

N  = length(x);
n  = N/2;
C  = 0.57721566490153286060651209008240243104215933593992; %euler-mascheroni's constant

rr = sqrt(bsxfun(@minus,x',x).^2+bsxfun(@minus,y',y).^2);
drr = dx.^2+dy.^2;
sdrr = sqrt(drr);




kernel = (1i/4)*besselh(0,1,kh*rr).*repmat(sdrr,N,1);
chi_kr = ones(size(kernel));
if(spars.ifsplit) 
    chi_kr = exp(-36*(kh*rr/spars.rfac/pi).^8);
end
%chi_kr = ones(size(kernel));
kernel_1 = -1/(4*pi)*besselj(0,kh*rr).*repmat(sdrr,N,1).*chi_kr;


kernel_2 = kernel - kernel_1.*log(4*sin(bsxfun(@minus,t',t)/2).^2);

kernel_2_diag = (1i/4-C/(2*pi)-1/(4*pi)*log(kh^2/4*(drr))).*sdrr;

kernel_2(1:N+1:end) = kernel_2_diag;

quad = layer_quad(n);

S = gallery('circul',quad).*kernel_1 + pi/n*kernel_2;
