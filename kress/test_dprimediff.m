%% Test dprimediff

addpath ./src
clear
a = 1.0;
b = 1.0;

src0 = [0.01;0.13];

gpars = [];
gpars.igeomtype = 6;
gpars.xscale = a;
gpars.yscale = b;

nppw = 6;

kh = 20;
kh2 = 1j*kh;

nhalf = ceil(nppw*kh);
n = 2*nhalf;
    
[srcinfo,h] = ellipse_new(a,b,n);
[src,ts,intpt,~] = get_geom(gpars,nhalf);


spars = [];
spars.ifsplit = true;
spars.rfac = 8;

uin = helm_c_p(kh,src0,srcinfo);
dudnin = helm_c_gn(kh,src0,srcinfo);
T = dprimelmat(kh,src,ts,spars);
T2 = dprimelmat(kh2,src,ts,spars);

Tdiff = T-T2;

Tdiff2 = dprimediffmat(kh,src,ts,spars);
rhs = dudnin;

v1 = Tdiff*rhs;
v2 = Tdiff2*rhs;

err1 = norm(v1-v2)/norm(v1);
fprintf('Error in dprimediff mat=%d\n',err1);

fprintf('Error in matrix norm=%d\n',norm(Tdiff2-Tdiff));


%% Test neumann solution with calderon identities

xx = srcinfo(1,:).';
yy = srcinfo(2,:).';

rnx = srcinfo(3,:).';
rny = srcinfo(4,:).';

dir = 0;
u_pw = exp(1i*kh*(xx*cos(dir) + yy*sin(dir)));
dudn_pw = 1i*kh*(rnx*cos(dir) + rny*sin(dir)).*u_pw;


S  = slmat(kh,src,ts,spars);
D  = dlmat(kh,src,ts,spars);
Sik = slmat(kh2,src,ts,spars);
Dik = dlmat(kh2,src,ts,spars);


eta = kh;
Bk = 1j*eta*(eye(n)/2-D) + Sik*T;


Bk2 = 1j*eta*(eye(n)/2-D) + Sik*Tdiff2 + -eye(n)/4 + Dik*Dik;

Bk3 = get_neu_mat(kh,src,ts,eta,spars);

rhs = 1j*eta*u_pw - Sik*dudn_pw;
sol = Bk\rhs;
utot = 2*(u_pw + D*sol);
uscat = utot-u_pw;


sol2 = Bk2\rhs;
utot2 = 2*(u_pw + D*sol2);
uscat2 = utot2-u_pw;


sol3 = Bk3\rhs;
utot3 = 2*(u_pw + D*sol3);
uscat3 = utot3-u_pw;


Bk4 = get_neu_mat_adj(kh,src,ts,eta,spars);
rhs_use = -dudn_pw;
sol4 = Bk4\rhs_use;
uscat4 = ((D + eye(n)/2)*Sik - 1j*eta*S)*sol4 ;

nn = -ceil(kh)-20:1:ceil(kh)+20;
jnm1 = besselj(nn-1,kh).';
jnp1 = besselj(nn+1,kh).';
imas = (1j).^(nn).';
zexp = exp(1j*ts.'*nn);


bn =  ((jnm1 - jnp1).*imas)*kh/2;

dudn_series = zexp*bn;

% Test accuracy with analytic solution
hn = besselh(nn,1,kh).';
hnp = (besselh(nn-1,1,kh).' - besselh(nn+1,1,kh).')/2;

cn = -bn.*hn/kh./hnp;
uscat_series = zexp*cn;

fprintf('Error in dudn =%d\n',norm(dudn_series-dudn_pw)/norm(dudn_series));
fprintf('Error in uscat no calderon=%d\n',norm(uscat-uscat_series)/norm(uscat_series));
fprintf('Error in uscat calderon=%d\n',norm(uscat2-uscat_series)/norm(uscat_series));
fprintf('Error in uscat get_neu_mat=%d\n',norm(uscat3-uscat_series)/norm(uscat_series));
fprintf('Error in uscat get_neu_mat_adj=%d\n',norm(uscat4-uscat_series)/norm(uscat_series));







