from numpy import *
from pylab import *

fname = 'diamond_data/diamond_res_aug11_thet005.txt'
xdir=loadtxt(fname)

zks = xdir[:,0]
noverk = xdir[:,2]
ppw = xdir[:,3]
errd = xdir[:,7]
errq = xdir[:,8]

cc = ['darkgreen','orangered','navy','magenta','black','cyan']

zkvals = [20.0*sqrt(2.0),40.0*sqrt(2.0),80.0*sqrt(2.0),160*sqrt(2.0),320*sqrt(2.0)]
nzk = size(zkvals)

figure(1)
clf()

figure(2)
clf()
for i in range(nzk):
    mask = abs(zks-zkvals[i])<1e-2
    noverkuse = noverk[mask]
    ppwuse = ppw[mask]
    errduse = errd[mask]
    errquse = errq[mask]
    figure(1)
    plot(ppwuse,errduse,'.',color=cc[i])
    figure(2)
    plot(ppwuse,1.0/errquse,'.',color=cc[i])
    
figure(1)    
legend([r'$k=20\sqrt{2}$',r'$k=40\sqrt{2}$',r'$k=80\sqrt{2}$',r'$k=160\sqrt{2}$',r'$k=320\sqrt{2}$'])
xlabel('Points per wavelength')
ylabel(r'$|\sigma-\sigma_{n}|$')
savefig('err_dens_diamond_dir_' + fname[25:-4] +'.pdf',bbox_inches='tight')
figure(2)    
legend([r'$k=20\sqrt{2}$',r'$k=40\sqrt{2}$',r'$k=80\sqrt{2}$',r'$k=160\sqrt{2}$',r'$k=320\sqrt{2}$'])
xlabel('Points per wavelength')
ylabel(r'$|\sigma-\sigma_{n}|/|(I-P_{N})\sigma|$')
savefig('err_cq_diamond_dir_' + fname[25:-4] + '.pdf',bbox_inches='tight')
show()
