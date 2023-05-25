addpath ./src/
A = load('data_postproc_aug11_2022_arclengthparam.mat');

% Plot rnorm inv, rnorm, err_pw, err_qnkqn for 
% all geometry types, 

geoms = ["Circle", "3-Star", "Cavity", "3-star equi", "Cavity-equi"];
splits = ["None", "s0=16", "s0=32"];


[nppw,nsplit,nkh,ng] = size(A.err_pw);
geoms = ["Circle", "3-Star", "Cavity", "3-star", "Cavity"];
for ig=1:ng
    figure(ig)
    clf
    
   for ik=1:1:nkh
       noverk = 2*A.ppw;
       dat = A.err_pw(:,1,ik,ig);
       semilogy(noverk,dat,'.','MarkerSize',20); hold on;
   end
   ylim([1e-14,1]);
   yticks([1e-12,1e-9,1e-6,1e-3,1]);
   ax = gca;
   ax.FontSize = 15;
   xlabel('Points per wavelength');
   xlim([3,13])
   ylabel('|\sigma - \sigma_{\rm{N}}|');
   
   title(geoms(ig),'FontSize',15);
   legend(sprintfc('k = %d', A.kh(1:1:end)));
   saveas(gcf,join(['./results/geom' num2str(ig) '_kh' num2str(A.kh(ik)) '_aug11_2022.pdf']));
   
   
end
