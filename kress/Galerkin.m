function [ks,norms]=Galerkin()
scale=2
k=140
ppw=3;
ks=[];
norms=[];
for j=1:3
    ks=[ks,k];
    norms=[norms,GalerkinRun(k,ppw)];
    k=k+10;
end
plot(ks,norms)

end

function [n]=GalerkinRun(k,ppw)
projection=true;



N=ppw*k;%num points
M=round(3*ppw*k);%maxFreq

lamt=@(m)lam(m,k);
A=lamt(-M:M);
if(projection)
P_N=zeros(2*M+1,2*M+1);
for j=1:2*M+1
    P_N(j,:)=PN(j-M-1,M,N)/sqrt(2*pi);
end

pert=diag(A);
op=eye(2*M+1)+P_N*pert;

minN=min(svd(op));
n=(norm(op\(eye(2*M+1)-P_N)));
end
minReal=min(abs(1+A));

end

function [colOut]=PN(m,M,N)
colOut=zeros(1,2*M+1);
for j=1:2*M+1
    if m==0
        colOut(M+1)=sqrt(2*pi);
    elseif j==M+1
    elseif mod(m,N)==mod(j-M-1,N)
        l=(j-M-1-m)/N;
        colOut(j)=1/pi/sqrt(2*pi)*N^2/m/(j-M-1)*2*(sin(pi*m/N))^2;
    end
end
end

function [lam]=lam(m,k)
etaN=1;
lam=zeros(size(m));
besselJp=@(m,k)-besselj(abs(m)+1,k)+1/k*m.*besselj(m,k);
besselHp=@(m,k)-besselh(m+1,k)+1/k*m.*besselh(m,k);
for j=1:length(m)
    if(abs(m(j))<1.1*k)
        lam(j)=(1i*etaN*(1-1i*pi/2*k*besselJp(abs(m(j)),k).*besselh(abs(m(j)),k))+1i/8*k^2*pi.*besselJp(abs(m(j)),k).*besselHp(abs(m(j)),k)./sqrt(m(j).^2+k^2))/(1i*etaN/2-1/4)-1;
    else
        lam(j)=0;
    end
end
end

