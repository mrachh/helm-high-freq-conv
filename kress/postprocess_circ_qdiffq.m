A = load('data_postproc_apr26_2022_circ_only.mat');

% Plot rnorm inv, rnorm, err_pw, err_qnkqn for 
% all geometry types, 

geoms = ["Circle", "3-Star", "Cavity"];
splits = ["None", "s0=16", "s0=32"];


[nppw,nsplit,nkh,ng] = size(A.err_pw);

rnorms = A.rnorms(nppw,1,:,1);
rnorms =  rnorms(:);

rnorms_diff = A.err_qnkqn3(nppw,1,:,1);
rnorms_diff = rnorms_diff(:);
rnorms_err = size(rnorms);
rnorms_diff_err = size(rnorms_diff);

ppw = A.ppw(nppw);

for i=1:nkh
    kh = A.kh(i);
    [qdiffq,qlnq] = get_circ_q_diff_ln_q(kh,ppw);
    rnorms_diff_err(i) = abs(max(abs(qdiffq))-rnorms_diff(i));
    rnorms_err(i) = abs(max(abs(qlnq))-rnorms(i));
end





figure(1)
clf
plot(A.kh,rnorms,'k.','MarkerSize',20);
xlabel('k','FontSize',15);
ylabel('|Q_{N}(I+ L_{N})Q_{N}|','FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,'./results/qnlnqn_norm.pdf');


figure(2)
clf
semilogy(A.kh,rnorms_err,'k.','MarkerSize',20);
xlabel('k','FontSize',15);
ylabel('\epsilon_{L_{N}}','FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,'./results/qnlnqn_norm_err.pdf');



figure(3)
clf
plot(A.kh,rnorms_diff,'k.','MarkerSize',20);
xlabel('k','FontSize',15);
ylabel('|Q_{N}(L- L_{N})Q_{N}|','FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,'./results/qndiffqn_norm.pdf');


figure(4)
clf
semilogy(A.kh,rnorms_diff_err,'k.','MarkerSize',20);
xlabel('k','FontSize',15);
ylabel('\epsilon_{L-L_{N}}','FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf,'./results/qndiffqn_norm_err.pdf');


figure(5)
clf
hold on;
for ippw=1:nppw
        
    pplot3 = squeeze(A.err_qnkqn3(ippw,1,:,1));
    semilogy(A.kh(:),pplot3(:),'.','MarkerSize',20); 
end

xlabel('k','FontSize',15);
ylabel('|Q_{N}(L- L_{N})Q_{N}|','FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
legend(sprintfc('ppw = %0.2f', A.ppw)); 
saveas(gcf,'./results/qndiffqn_norm_ppw.pdf');

