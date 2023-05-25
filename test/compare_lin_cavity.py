from numpy import *
from pylab import *

fname = 'cavity_data/cavity_res_oct5_pw_mpi2p0p2.txt'
xdir=loadtxt(fname)

fname2 = 'cavity_data/cavity_res_nov15_pw_mpi2p0p2_lin.txt'
xdir2=loadtxt(fname2)


zks = xdir[:,0]
noverk = xdir[:,2]
ppw = xdir[:,3]
errd = xdir[:,7]
errq = xdir[:,8]


zks2 = xdir2[:,0]
noverk2 = xdir2[:,2]
ppw2 = xdir2[:,3]
errd2 = xdir2[:,7]
errq2 = xdir2[:,8]

cc = ['darkgreen','orangered','navy','magenta','black','cyan','green','yellow','pink','teal']

zkvals = [37.213]
nzk = size(zkvals)

rsc = 40.61   # norm of density, computed separately
figure(1)
clf()

figure(2)
clf()
for i in range(nzk):
    mask = abs(zks-zkvals[i])<1e-5
    noverkuse = noverk[mask]
    ppwuse = ppw[mask]
    errduse = errd[mask]/rsc
    errquse = errq[mask]

    mask2 = abs(zks2-zkvals[i])<1e-5
    noverkuse2 = noverk2[mask2]
    ppwuse2 = ppw2[mask2]
    errduse2 = errd2[mask2]/rsc
    errquse2 = errq2[mask2]
    figure(1)
    semilogy(ppwuse,errduse,'.',color=cc[2*i])
    semilogy(ppwuse2,errduse2,'.',color=cc[2*i+1])
    figure(2)
    semilogy(ppwuse,1.0/errquse,'.',color=cc[2*i])
    semilogy(ppwuse2,1.0/errquse2,'.',color=cc[2*i+1])
    
figure(1)    
legend(['Const','Linear'])
xlabel('Points per wavelength')
ylabel(r'$|\sigma-\sigma_{n}|$/$|\sigma|$')
xlim([0,50])
savefig('err_dens_cavity_dir_pw_mpi2p0p2_lincomp.pdf',bbox_inches='tight')
figure(2)    
legend(['Const','Linear'])
xlabel('Points per wavelength')
ylabel(r'$|\sigma-\sigma_{n}|/|(I-P_{N})\sigma|$')
xlim([0,50])
savefig('err_cq_cavity_dir_pw_mpi2p0p2_lincomp.pdf',bbox_inches='tight')
show()
