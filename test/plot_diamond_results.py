from numpy import *
from pylab import *

xdir=loadtxt('diamond_data/diamond_res3.txt')

zks = xdir[:,0]
noverk = xdir[:,2]
errd = xdir[:,7]
errq = xdir[:,8]

cc = ['darkgreen','orangered','navy','magenta','black','cyan']

nzk = 1
nppw = 20
figure(1)
clf()

figure(2)
clf()
for i in range(nzk):
    istart = i*nppw
    iend = (i+1)*nppw
    figure(1)
    plot(noverk[istart:iend],errd[istart:iend],'.',color=cc[i])
    figure(2)
    plot(noverk[istart:iend],1.0/errq[istart:iend],'.',color=cc[i])
    
figure(1)    
legend([r'$k=10\sqrt{2}$',r'$k=20\sqrt{2}$',r'$k=40\sqrt{2}$',r'$k=80\sqrt{2}$',r'$k=160\sqrt{2}$',r'$k=320\sqrt{2}$'])
legend([r'$k=10\sqrt{2}$',r'$k=20\sqrt{2}$'])
legend([r'$k=10\sqrt{2}$'])
xlabel('N/k')
ylabel(r'$|\sigma-\sigma_{n}|$')
savefig('err_dens_diamond_dir_wp.pdf',bbox_inches='tight')
figure(2)    
legend([r'$k=10\sqrt{2}$',r'$k=20\sqrt{2}$',r'$k=40\sqrt{2}$',r'$k=80\sqrt{2}$',r'$k=160\sqrt{2}$',r'$k=320\sqrt{2}$'])
legend([r'$k=10\sqrt{2}$'])
ylabel(r'$|\sigma-\sigma_{n}|/|(I-P_{N})\sigma|$')
savefig('err_cq_diamond_dir_wp.pdf',bbox_inches='tight')
show()
