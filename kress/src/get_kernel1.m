function [kernel_1] = get_kernel1(kh,src,zpars,spars)
%
%  This function returns the kernel multiplying 
%  the log term in the helmholtz layer potentials
%
%  Input arguments:
%    kh: double
%       wave number
%    src: double (8,n)
%       source information, only the first 4 components of src
%       are used
%    t: double(n,1)
%       parameter space discretization
%    zpars: complex(2,1)
%      zpars(1): single layer strength
%      zpars(2): double layer strength
%    spars: splitting params


    if(nargin <= 2)
        zpars = [-1j*kh 1];
        %zpars = [1 0];
    end
    if(nargin <= 3) 
        spars = [];
        spars.ifsplit = false;
        spars.rfac = 16;
    end

    xs  = src(1,:);
    ys  = src(2,:);
    dx  = src(3,:);
    dy  = src(4,:);

    N  = length(xs);

    rr = sqrt(bsxfun(@minus,xs',xs).^2+bsxfun(@minus,ys',ys).^2);
    drr = dx.^2+dy.^2;
    sdrr = sqrt(drr);

    kernel = (1i*kh/4)*(bsxfun(@minus,xs',xs).*repmat(dy,N,1)-bsxfun(@minus,ys',ys).*repmat(dx,N,1)).* ...
        besselh(1,1,kh*rr)./rr;
    chi_kr = ones(size(kernel));
    if(spars.ifsplit) 
        chi_kr = exp(-36*(kh*rr/spars.rfac/pi).^8);
    end
    
    j1s = besselj(1,kh*rr)./rr;
    j1s(1:N+1:end) = ones(N,1)*kh/2; 
    kernel_1_d = -kh/(4*pi)*(bsxfun(@minus,xs',xs).*repmat(dy,N,1)-bsxfun(@minus,ys',ys).*repmat(dx,N,1)).* ...
        j1s;
    kernel_1_s = -1/(4*pi)*besselj(0,kh*rr).*repmat(sdrr,N,1);
    kernel_1 = zpars(1)*kernel_1_s + zpars(2)*kernel_1_d;
    kernel_1 = kernel_1.*chi_kr;

end