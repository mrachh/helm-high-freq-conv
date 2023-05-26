zk = 37.212733389; % Wavenumber

% zk = 3.12;

nch = 200; % number of panels, each panel is 16th order

% setup target locations
x = -4:0.01:4;
[xx1,yy1] = meshgrid(x);
xx = xx1(:);
yy = yy1(:);

ntarg = length(xx);
targs = zeros(2,ntarg);
targs(1,:) = xx;
targs(2,:) = yy;

thet = -pi/2 + 0.2;  % plane wave angle
is_analytic = 1;     % flag for performing analytic test

n = nch*16;
m = 200;

src_info = geometries.larrycup(0.2,pi/12,n,m);

nh = 0;
hcoefs = zeros(1,1);
[src_out] = rla.update_geom(src_info,nh,hcoefs);






fid = fopen('input_data.dat','w');
fprintf(fid,'%12.8f\n',zk);
fprintf(fid,'%d\n',nch);
fprintf(fid,'%d\n',ntarg);
fprintf(fid,'%12.8f\n',thet);
fprintf(fid,'%d\n',is_analytic);
fclose(fid);

fid = fopen('target_data.dat','w');
fprintf(fid,'%12.8f  %12.8f \n',targs);
fclose(fid);

!./run_cavity_reference_sol.sh

p1 = load('potential.dat');

pot = p1(:,1) + 1j*p1(:,2);
uinc = p1(:,3) + 1j*p1(:,4);

if(is_analytic)

    p2 = load('potential_analytic.dat');


    pot_analytic = p2(:,1) + 1j*p2(:,2);
    pot_ex = p2(:,3) + 1j*p2(:,4);


    pplot = reshape(log10(abs(pot_ex-pot_analytic)),size(xx1));
    figure
    clf
    pcolor(xx1,yy1,pplot); shading interp; colorbar(); 
end

pplot = reshape(real(uinc),size(xx1));
figure
clf
pcolor(xx1,yy1,pplot); shading interp; colorbar();
hold on; fill(src_out.xs,src_out.ys,'w-');





