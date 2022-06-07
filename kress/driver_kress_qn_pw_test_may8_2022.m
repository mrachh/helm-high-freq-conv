addpath './src';
addpath '~/git/inverse-obstacle-scattering2d/src/';
kh = 20:10:60;
ppw = [2.0 3.0 4.0 6.0];

nkh = length(kh);
nppw = length(ppw);

ng = 5;
nsplit = 3;

gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;

ifsplit = [false true true];
rfacs = [16 16 32];

dirs = [0 0 pi/3 0 pi/3];

nuse = 1500;

errs = zeros(nppw,nsplit,nkh,ng);
rnorms = zeros(nppw,nsplit,nkh,ng);
rnorms_inv = zeros(nppw,nsplit,nkh,ng);
err_pw = zeros(nppw,nsplit,nkh,ng);
err_qnkqn = zeros(nppw,nsplit,nkh,ng);
err_qnkqn2 = zeros(nppw,nsplit,nkh,ng);
err_qnkqn3 = zeros(nppw,nsplit,nkh,ng);

ref_sols_all = cell(nkh,ng);

for ig=1:ng
    gpars = [];
    gpars = gpars0;
    gpars.igeomtype = ig;
    
    bd_params = [];
    bd_params.dir = dirs(ig);
    for ikh = 1:nkh
        fprintf('\n\n\nStarting reference solution for ikh=%d, ig=%d\n\n\n',ikh,ig);
        [ref_sols] = get_ref_sols(gpars,kh(ikh),bd_params,20,nuse);
        ref_sols_all{ikh,ig} = ref_sols;
    
        for isplit=1:nsplit
            spars = [];
            spars.ifsplit = ifsplit(isplit);
            spars.rfac = rfacs(isplit);
            for ippw=1:nppw
                [out] = dens_interp_err(ref_sols,gpars,kh(ikh),ppw(ippw),spars);
                errs(ippw,isplit,ikh,ig) = out.err_ex;
                err_pw(ippw,isplit,ikh,ig) = out.err_dens;
                rnorms(ippw,isplit,ikh,ig) = out.rnorm;
                rnorms_inv(ippw,isplit,ikh,ig) = out.rnorm_inv;
                err_qnkqn(ippw,isplit,ikh,ig) = out.err_qnkqn;
                err_qnkqn2(ippw,isplit,ikh,ig) = out.err_qnkqn2;
                err_qnkqn3(ippw,isplit,ikh,ig) = out.err_qnkqn3;
                
            end
        end
        
                
    end
end
    
save('data_postproc_may8_2022_arclengthparam.mat','err_pw','err_qnkqn','errs','gpars0',...
    'ifsplit','kh','ppw','rfacs','rnorms','rnorms_inv','err_qnkqn2','err_qnkqn3');


save('data_all_may8_2022_arclengthparam.mat','err_pw','err_qnkqn','errs','gpars0',...
    'ifsplit','kh','ppw','rfacs','rnorms','rnorms_inv','err_qnkqn2','err_qnkqn3','ref_sols_all','-v7.3');
