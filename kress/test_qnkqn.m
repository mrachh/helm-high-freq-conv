addpath 'src/'

kh = 20.0;
n = 100;

m1 = 800;
m2 = 900;


gpars = [];
gpars.igeomtype = 1;
gpars.a = 0.3;
gpars.b = pi/2;
gpars.n0 = 100;
[srcn,tn] = get_geom(gpars,n);
[srcm1,tm1] = get_geom(gpars,m1);
[srcm2,tm2] = get_geom(gpars,m2);
qn = get_qmat(n);
qm1 = get_qmat(m1);
qm2 = get_qmat(m2);
enm1 = get_emat(m1,n);
em1n = get_emat(n,m1);

enm2 = get_emat(m2,n);
em2n = get_emat(n,m2);

kn = get_kress_mat(kh,srcn,tn);
km1 = get_kress_mat(kh,srcm1,tm1);
km2 = get_kress_mat(kh,srcm2,tm2);

qnknqn = qn*kn;
qnkqn_m1 = qn*enm1*qm1*km1*em1n*qn;
qnkqn_m2 = qn*enm2*qm2*km2*em2n*qn;

err1 = norm(qnkqn_m2-qnkqn_m1)/norm(qnkqn_m2);
err2 = norm(qnknqn-qnkqn_m2);

fprintf('self convergence error: %d\n',err1);
fprintf('norm qnkqn: %d\n',err2);


f1 = sin(10*tm1);
f2 = sin(10*tn);
f3 = em1n*qn*f2.';
err3 = norm(f1-f3.');
fprintf('error in interpolation %d\n',err3);
