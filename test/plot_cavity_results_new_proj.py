from numpy import *
from pylab import *

fname = 'cavity_data/cavity_res_oct5_pw_mpi2p0p2.txt'
xdir=loadtxt(fname)

zks = xdir[:,0]
noverk = xdir[:,2]
ppw = xdir[:,3]
errd = xdir[:,7]
errq = xdir[:,8]

cc = ['darkgreen','orangered','navy','magenta','black','cyan','green','yellow','pink','teal']

zkvals = [20.0,80.0,160.0,37.213,30.635,47.020,56.690,64.539,37.212]
nzk = size(zkvals)

figure(1)
clf()

figure(2)
clf()
for i in range(nzk):
    mask = abs(zks-zkvals[i])<1e-5
    noverkuse = noverk[mask]
    ppwuse = ppw[mask]
    errduse = errd[mask]
    errquse = errq[mask]
    figure(1)
    semilogy(ppwuse,errduse,'.',color=cc[i])
    figure(2)
    semilogy(ppwuse,1.0/errquse,'.',color=cc[i])
    
figure(1)    
legend([r'$k=20$',r'$k=80$',r'$k=160$',r'$k=37.213$',r'$k=30.635$',r'$k=47.020$',r'$k=56.690$',r'$k=64.539$',r'$k=37.212$'])
xlabel('Points per wavelength')
ylabel(r'$|\sigma-\sigma_{n}|$')
savefig('err_dens_cavity_dir_pw_mpi2p0p2.pdf',bbox_inches='tight')
figure(2)    
legend([r'$k=20$',r'$k=80$',r'$k=160$',r'$k=37.213$',r'$k=30.635$',r'$k=47.020$',r'$k=56.690$',r'$k=64.539$',r'$k=37.212$'])
xlabel('Points per wavelength')
ylabel(r'$|\sigma-\sigma_{n}|/|(I-P_{N})\sigma|$')
savefig('err_cq_cavity_dir_pw_mpi2p0p2.pdf',bbox_inches='tight')
show()
