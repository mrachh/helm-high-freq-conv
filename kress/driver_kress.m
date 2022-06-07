%driver 
%global variables for old code
clear


kh = 1:1:30;
ppw = 3:1:5;
nkh = length(kh);
nppw = length(ppw);

errs = zeros(nkh,nppw);
rnorms = zeros(nkh,nppw);
rinvs = zeros(nkh,nppw);
for j=1:nppw
    for i=1:nkh
        nmax = 10*kh(i);
        [errs(i,j),~,~,rnorms(i,j),rinvs(i,j)] = kress_circ_err(kh(i),ppw(j),nmax);
    end
end

figure(1)
clf
semilogy(kh,errs(:,1),'k.','MarkerSize',20); hold on;
semilogy(kh,errs(:,2),'b.','MarkerSize',20); hold on;
semilogy(kh,errs(:,3),'r.','MarkerSize',20); hold on;
legend('C=3','C=4','C=5');


figure(2)
clf
semilogy(kh,rnorms(:,1),'k.','MarkerSize',20); hold on;
semilogy(kh,rnorms(:,2),'b.','MarkerSize',20); hold on;
semilogy(kh,rnorms(:,3),'r.','MarkerSize',20); hold on;
legend('C=3','C=4','C=5','Location','SouthEast');


figure(3)
clf
semilogy(kh,rinvs(:,1),'k.','MarkerSize',20); hold on;
semilogy(kh,rinvs(:,2),'b.','MarkerSize',20); hold on;
semilogy(kh,rinvs(:,3),'r.','MarkerSize',20); hold on;
legend('C=3','C=4','C=5','Location','SouthEast');



