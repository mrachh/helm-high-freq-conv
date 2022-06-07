addpath './src';
addpath '~/git/inverse-obstacle-scattering2d/src/';
kh = 20:10:60;
ppw = [2.0 3.0 4.0 6.0 20.0];
spows = [1.5 2 4];

nkh = length(kh);
nppw = length(ppw);
nspow = length(spows);

ng = 9;
nsplit = 3;



gpars0 = [];
gpars0.a = 0.2;
gpars0.b = pi/12;
gpars0.n0 = 200;

gpars0.nosc = 3;
gpars0.rfac = 0.3;


gpars0.xscale = 1.2;
gpars0.yscale = 1.1;

gpars0.alpha = 0.2;
gpars0.beta = 0.75;
gpars0.gamma = 0;

ifsplit = [false true true];
rfacs = [16 16 32];

dirs = [0 0 pi/3 0 pi/3 0 0 0 0];

nuse = 1500;
val = zeros(nppw,nsplit,nspow,nkh,ng);

for ig=6:ng
    gpars = gpars0;
    gpars.igeomtype = ig;
    
    bd_params = [];
    bd_params.dir = dirs(ig);
    for ikh = 1:nkh
        fprintf('\n\n\nStarting reference solution for ikh=%d, ig=%d\n\n\n',ikh,ig);
        
        for isplit=1:nsplit
            spars = [];
            spars.ifsplit = ifsplit(isplit);
            spars.rfac = rfacs(isplit);
            for ispow=1:nspow
                spow = spows(ispow); 
                for ippw=1:nppw
                    
                    [out] = get_fourier_ker1_hsnorm(gpars,kh(ikh),ppw(ippw),spars,spow);
                    val(ippw,isplit,ispow,ikh,ig) = out.val;
                end
            end
        end
        
                
    end
end


save('data_kernel1_hs_jun1_2022_newscaling.mat','gpars0','ifsplit','kh','ppw', ...
  'spows','rfacs','val');
