addpath './src';
kh = 20:5:60;
ppw = [2.0 2.2 2.5 3.0 4.0 6.0];

nkh = length(kh);
nppw = length(ppw);

ng = 3;
nsplit = 3;

gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;

ifsplit = [false true true];
rfacs = [16 16 32];

dirs = [0 0 pi/3];

nuse = 1500;

errs = zeros(nppw,nsplit,nkh,ng);
rnorms = zeros(nppw,nsplit,nkh,ng);
rnorms_inv = zeros(nppw,nsplit,nkh,ng);
err_pw = zeros(nppw,nsplit,nkh,ng);
err_qnkqn = zeros(nppw,nsplit,nkh,ng);

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
                
            end
        end
        
                
    end
end
    
