function [mat] = get_neu_mat(kh,src,ts,eta,spars)
%
%  this function returns the matrix corresponding to a neumann 
%  boundary value problem with kress quadratures and calderon
%  identites, and where the self part of the log is scaled 
%  by a bump function around the identity
%

if(nargin == 3)
    eta = kh;
    spars = [];
    spars.ifsplit = true;
    spars.rfac = 8;
elseif(nargin == 4)
    spars = [];
    spars.ifsplit = true;
    spars.rfac = 8;
end

kh2 = 1j*kh;
[~,n] = size(src);
D  = dlmat(kh,src,ts,spars);
Sik = slmat(kh2,src,ts,spars);
Dik = dlmat(kh2,src,ts,spars);
Tdiff = dprimediffmat(kh,src,ts,spars);

mat = 1j*eta*(eye(n)/2-D) + Sik*Tdiff + -eye(n)/4 + Dik*Dik;


end