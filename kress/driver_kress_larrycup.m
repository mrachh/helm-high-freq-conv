%driver 
%global variables for old code
clear
addpath ~/git/finufft/matlab/

kh = 20:4:60;
ppw = (3.0:0.5:4.5)*2;
nkh = length(kh);
nppw = length(ppw);

errs = zeros(nkh,nppw);
err_geom = zeros(nkh,nppw);
rnorms = zeros(nkh,nppw);
rinvs = zeros(nkh,nppw);


nmin = 400;
a = 0.3;
b = pi/12;
ppw0 = 40;

for i=1:nkh
    [sol_c_ref,tref,dsdt_ref] = kress_ref_larrycup(kh(i),a,b,nmin,ppw0);
    for j=1:nppw
        [errs(i,j),err_geom(i,j),~,rnorms(i,j),rinvs(i,j)] = ...
          kress_err_larrycup(kh(i),ppw(j),a,...
          b,tref,sol_c_ref,dsdt_ref);
        
    end
end


figure(4)
subplot(2,2,1)

semilogy(kh,errs(:,1),'k.','MarkerSize',20); hold on;
semilogy(kh,errs(:,2),'b.','MarkerSize',20); hold on;
semilogy(kh,errs(:,3),'r.','MarkerSize',20); hold on;
semilogy(kh,errs(:,4),'m.','MarkerSize',20); hold on;
title('Kress density error');

subplot(2,2,2)
semilogy(kh,err_geom(:,1),'k.','MarkerSize',20); hold on;
semilogy(kh,err_geom(:,2),'b.','MarkerSize',20); hold on;
semilogy(kh,err_geom(:,3),'r.','MarkerSize',20); hold on;
semilogy(kh,err_geom(:,4),'m.','MarkerSize',20); hold on;
title('Kress geometry error');
legend('C=2','C=2.5','C=3','C=3.5');


subplot(2,2,3)
semilogy(kh,rnorms(:,1),'k.','MarkerSize',20); hold on;
semilogy(kh,rnorms(:,2),'b.','MarkerSize',20); hold on;
semilogy(kh,rnorms(:,3),'r.','MarkerSize',20); hold on;
semilogy(kh,rnorms(:,4),'m.','MarkerSize',20); hold on;
legend('C=2','C=2.5','C=3','C=3.5','Location','SouthEast');
title('Kress matrix norm');


subplot(2,2,4)
semilogy(kh,rinvs(:,1),'k.','MarkerSize',20); hold on;
semilogy(kh,rinvs(:,2),'b.','MarkerSize',20); hold on;
semilogy(kh,rinvs(:,3),'r.','MarkerSize',20); hold on;
semilogy(kh,rinvs(:,4),'m.','MarkerSize',20); hold on;
legend('C=2','C=2.5','C=3','C=3.5','Location','SouthEast');
title('Kress inverse norm');






