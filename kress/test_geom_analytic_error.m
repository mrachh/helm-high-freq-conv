addpath('./src')
gpars = [];
gpars.a = 0.2;
gpars.b = pi/12;
gpars.n0 = 200;

gpars.nosc = 3;
gpars.rfac = 0.3;

gpars.xscale = 1.2;
gpars.yscale = 1.1;

gpars.alpha = 0.2;
gpars.beta = 0.75;
gpars.gamma = 0;

gpars.igeomtype = 9;

kh = 10;

[ref_sols] = get_ref_sols(gpars,kh);
fprintf('Error in analytic test=%d\n',ref_sols.err_ex);



