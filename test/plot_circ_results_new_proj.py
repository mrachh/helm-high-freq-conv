from numpy import *
from pylab import *

xdir=loadtxt('circ_data/circ_res_newproj2.txt')

zks = xdir[:,0]
noverk = xdir[:,2]
errd = xdir[:,7]
errq = xdir[:,8]

cc = ['darkgreen','orangered','navy','magenta','black','cyan']

zkvals = [10.0,20.0,40.0,80.0,160.0]
nzk = size(zkvals)

figure(1)
clf()

figure(2)
clf()
for i in range(nzk):
    mask = abs(zks-zkvals[i])<1e-2
    noverkuse = noverk[mask]
    errduse = errd[mask]
    errquse = errq[mask]
    figure(1)
    plot(noverkuse,errduse,'.',color=cc[i])
    figure(2)
    plot(noverkuse,1.0/errquse,'.',color=cc[i])
    
figure(1)    
legend([r'$k=10$',r'$k=20$',r'$k=40$',r'$k=80$',r'$k=160$'])
xlabel('N/k')
ylabel(r'$|\sigma-\sigma_{n}|$')
savefig('err_dens_circ_dir_new_proj_thet5.pdf',bbox_inches='tight')
figure(2)    
legend([r'$k=10$',r'$k=20$',r'$k=40$',r'$k=80$',r'$k=160$'])
ylabel(r'$|\sigma-\sigma_{n}|/|(I-P_{N})\sigma|$')
savefig('err_cq_circ_dir_new_proj_thet5.pdf',bbox_inches='tight')
show()
