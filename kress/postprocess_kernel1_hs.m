addpath ./src/
A = load('data_kernel1_hs_jun1_2022_newscaling.mat');

% Plot rnorm inv, rnorm, err_pw, err_qnkqn for 
% all geometry types, 

geoms = ["Circle", "3-Star", "Cavity", "3-star equi", "Cavity-equi","Ellipse","Ellipse-equi","Ice-cream cone","Ice-cream cone equi"];
splits = ["None", "s0=16", "s0=32"];


[nppw,nsplit,nspow,nkh,ng] = size(A.val);

mms = cell(1,nspow);
mms{1}.marker = '.';
mms{1}.ms = 20;
mms{2}.marker = 's';
mms{2}.ms = 10;
mms{3}.marker = '*';
mms{3}.ms = 10;

colors = cell(1,nppw);
colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'm';
colors{4} = 'k';

isubfig = [1 2 4];
for ig=6:ng
    ifig = ig;
    figure(ifig)
    clf
    sgtitle(join(['Geometry:' geoms(ig)]),'FontSize',20);
    for isplit=1:nsplit
        subplot(2,3,isubfig(isplit)); hold on;
        for ippw=1:nppw
            for ispow=1:nspow
               muse = [colors{ippw} mms{ispow}.marker];
               pplot = squeeze(A.val(ippw,isplit,ispow,:,ig));
               semilogy(A.kh(:),pplot(:),muse,'MarkerSize',mms{ispow}.ms); 
            end
        end
        title(join(['Split: ' splits(isplit)]),'FontSize',15)
        set(gca,'YScale','log')
        ylim([1e-8,10])
    end
    subplot (2, 3, [3 6]) % merge remaining subplots and put legend here
    plot(A.kh(:), nan,'.','MarkerSize',15) % plot nans (hack to generate correct legend but plot no data)



    legend(sprintfc('ppw = %0.2f', A.ppw), 'Location', 'west','FontSize',15); 
    axis off
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    saveas(gcf,join(['./results/hs_norms_geom' num2str(ig) '_jun1_2022.pdf']));
end



gpars = A.gpars0;
gpars.igeomtype = 6;
src_6 = get_geom(gpars,300);
gpars.igeomtype = 8;
src_8 = get_geom(gpars,300);
figure
clf
subplot(1,2,1)
plot(src_6(1,:),src_6(2,:),'k.')
axis equal
axis off
title(geoms(6),'FontSize',10)

subplot(1,2,2)
plot(src_8(1,:),src_8(2,:),'k.')
axis equal
axis off
title(geoms(8),'FontSize',10)
saveas(gcf,'./results/geometries_new.pdf');

