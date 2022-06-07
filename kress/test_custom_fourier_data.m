

n0 = 3;
n1 = 6;
n2 = 13;
kdecay = 5;
tol = 1e-8;

beta = 1/(n2-n1)*(log(1/tol) + kdecay*log(n0/n1));

kvec1 = -n0:n0;
fvec1 = ones(size(kvec1));
kvec2 = -n1:(-n0-1);
fvec2 = abs(n0./kvec2).^kdecay;
kvec3 = (n0+1):n1;
fvec3 = abs(n0./kvec3).^kdecay;
kvec4 = -n2:(-n1-1);
fvec4 = (n0/n1)^kdecay*exp(-beta*(abs(kvec4)-n1));
kvec5 = (n1+1):n2;
fvec5 = (n0/n1)^kdecay*exp(-beta*(abs(kvec5)-n1));

kvec = [kvec1 kvec2 kvec3 kvec4 kvec5];
fvec = [fvec1 fvec2 fvec3 fvec4 fvec5];

nnuse = 200;
ts = 0:2*pi/nnuse:2*pi*(1-1/nnuse);
isign = -1;
tol2 = 1e-15;
f = real(finufft1d3(kvec,fvec,isign,tol2,ts))/length(kvec);

figure(1)
clf
plot(ts,f,'k.');
ff = fft(f);

figure(2)
clf
semilogy(abs(ff),'r.')

