A = load('norm_inv_cavity_data_jun13_2022_ima0_ncheb1024.mat');
B = load('norm_inv_cavity_data_jun13_2022.mat');


[n,npan] = size(A.detchebs);


f = cell(npan,1);
g = cell(npan,1);
figure
clf

for i=1:npan
    f{i} = chebfun(log10(A.detchebs2(:,i)),[A.kh0 + (i-1)*A.panlen, A.kh0 + i*A.panlen]);
    r = roots(diff(f{i}));
    subplot(3,3,i);
    plot(f{i},'b-'); hold on;
    plot(r,f{i}(r),'.g')
    g{i} = chebfun(log10(B.detchebs2(:,i)),[B.kh0 + (i-1)*B.panlen, B.kh0 + i*B.panlen]);
    r2 = roots(diff(g{i}));
    plot(g{i},'k-'); hold on;
    plot(r2,g{i}(r2),'.r')
    xlabel('Helmholtz wavenumber');
    ylabel('log10(Norm of inverse)');
end
set(gcf,'Position',[10 10 1500 1500]);
print(gcf,'-bestfit','-dpdf','inverse_norm_cavity_new.pdf');
