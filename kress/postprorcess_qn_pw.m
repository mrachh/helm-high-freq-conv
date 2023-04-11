addpath ./src/
A = load('data_postproc_aug11_2022_arclengthparam.mat');

% Plot rnorm inv, rnorm, err_pw, err_qnkqn for 
% all geometry types, 

geoms = ["Circle", "3-Star", "Cavity", "3-star equi", "Cavity-equi"];
splits = ["None", "s0=16", "s0=32"];


[nppw,nsplit,nkh,ng] = size(A.err_pw);


for ig=1:ng
    for isplit=1:nsplit
        ifig = (ig-1)*nsplit + isplit;
        figure(ifig)
        clf
        sgtitle(join(['Geometry:' geoms(ig) '    Split: ' splits(isplit)]),'FontSize',20);
        subplot(2,3,1); hold on;
        for ippw=1:nppw
            pplot = squeeze(A.err_pw(ippw,isplit,:,ig));
            semilogy(A.kh(:),pplot(:),'.','MarkerSize',20); 
        end
        set(gca,'YScale','log')
        title('Plane wave density error','FontSize',10)
        ylim([1e-14,1])
        
        subplot(2,3,2); hold on;
        for ippw=1:nppw
            %pplot = squeeze(A.err_qnkqn(ippw,isplit,:,ig));
            %pplot2 = squeeze(A.err_qnkqn2(ippw,isplit,:,ig));
            pplot3 = squeeze(A.err_qnkqn3(ippw,isplit,:,ig));
            %semilogy(A.kh(:),pplot(:),'.','MarkerSize',20); 
            %semilogy(A.kh(:),pplot2(:),'s','MarkerSize',5); 
            semilogy(A.kh(:),pplot3(:),'.','MarkerSize',20); 
            
        end
        set(gca,'YScale','log')
        title('|Q_{n}(L-L_{n})Q_{n}|','FontSize',10)
        ylim([1e-6,1])
        
        
        subplot(2,3,4); hold on;
        for ippw=1:nppw
            pplot = squeeze(A.rnorms(ippw,isplit,:,ig));
            plot(A.kh(:),pplot(:),'.','MarkerSize',20); 
        end
        title('|I+L_{n}|','FontSize',10);
        set(gca,'YScale','linear')
        
        subplot(2,3,5); hold on;
        for ippw=1:nppw
            pplot = squeeze(A.rnorms_inv(ippw,isplit,:,ig));
            plot(A.kh(:),pplot(:),'.','MarkerSize',20); 
        end
        title('|I+L_{n}|^{-1}','FontSize',10);
        set(gca,'YScale','linear')
        
        
        subplot (2, 3, [3 6]) % merge remaining subplots and put legend here
        plot(A.kh(:), nan,'.','MarkerSize',15) % plot nans (hack to generate correct legend but plot no data)
        
        

        legend(sprintfc('ppw = %0.2f', A.ppw), 'Location', 'west','FontSize',15); 
        axis off
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        saveas(gcf,join(['./results/geom' num2str(ig) '_isplit' num2str(isplit) '_aug11_2022.pdf']));
        
    end
end

gpars = A.gpars0;
gpars.igeomtype = 1;
src_1 = get_geom(gpars,300);
gpars.igeomtype = 2;
src_2 = get_geom(gpars,300);
gpars.igeomtype = 3;
src_3 = get_geom(gpars,300);
figure(16)
clf
subplot(1,3,1)
plot(src_1(1,:),src_1(2,:),'k.')
axis equal
axis off
title(geoms(1),'FontSize',10)

subplot(1,3,2)
plot(src_2(1,:),src_2(2,:),'k.')
axis equal
axis off
title(geoms(2),'FontSize',10)


subplot(1,3,3)
plot(src_3(1,:),src_3(2,:),'k.')
axis equal
axis off
title(geoms(3),'FontSize',10)
saveas(gcf,'./results/geometries.pdf');
