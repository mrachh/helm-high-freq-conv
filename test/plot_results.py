from numpy import *
from pylab import *

xdir=loadtxt('circ_data/circ_res.txt')
xneu=loadtxt('circ_neu_data/circ_res.txt')
xneu2=loadtxt('circ_neu_data/circ_res_pn.txt')

zks = xdir[:,0]
noverk = xdir[:,2]
errd = xdir[:,7]
errq = xdir[:,8]

cc = ['darkgreen','orangered','navy','magenta','black']

nzk = 5
nppw = 20
figure(1)
clf()

figure(2)
clf()
for i in range(nzk-1):
    istart = i*nppw
    iend = (i+1)*nppw
    figure(1)
    plot(noverk[istart+1:iend],errd[istart+1:iend],'.',color=cc[i])
    figure(2)
    plot(noverk[istart+1:iend],1.0/errq[istart+1:iend],'.',color=cc[i])
    
figure(1)    
legend(['k=10','k=20','k=40','k=80'])
xlabel('N/k')
ylabel(r'$|\sigma-\sigma_{n}|$')
savefig('err_dens_dir.pdf',bbox_inches='tight')
title('Dirichlet')
figure(2)    
legend(['k=10','k=20','k=40','k=80'])
ylabel(r'$|\sigma-\sigma_{n}|/|(I-P_{N})\sigma|$')
title('Dirichlet')
savefig('err_cq_dir.pdf',bbox_inches='tight')



zks = xneu[:,0]
noverk = xneu[:,2]
errd = xneu[:,7]
errq = xneu[:,8]
print(shape(errd))
print(shape(errq))
print(shape(noverk))

cc = ['darkgreen','orangered','navy','magenta','black','cyan','green']

nzk = 4
nppw = 20
figure(3)
clf()

figure(4)
clf()
for i in arange(0,7):
    istart = i*nppw
    iend = (i+1)*nppw
    print(istart,iend)
    figure(3)
    plot(noverk[istart+1:iend],errd[istart+1:iend],'.',color=cc[i])
    figure(4)
    plot(noverk[istart+1:iend],1.0/errq[istart+1:iend],'.',color=cc[i])
    
    
figure(3)    
#legend(['k=10','k=20','k=40','k=80'])
legend(['k=10','k=20','k=40','k=80','k=1280','k=2560','k=5120'])
xlabel('N/k')
ylabel(r'$|\sigma-\sigma_{n}|$')
savefig('err_dens_neu_5120.pdf',bbox_inches='tight')
title('Neumann')
figure(4)    
xlabel('N/k')
legend(['k=10','k=20','k=40','k=80','k=1280','k=2560','k=5120'])
#legend(['k=10','k=20','k=40','k=80'])
ylabel(r'$|\sigma-\sigma_{n}|/|(I-P_{N})\sigma|$')
title('Neumann')
savefig('err_cq_neu_5120.pdf',bbox_inches='tight')

show()



zks = xneu2[:,0]
noverk = xneu2[:,1]/xneu2[:,0]
errd = xneu2[:,7]
errq = xneu2[:,8]
print(shape(errd))
print(shape(errq))
print(shape(noverk))

cc = ['darkgreen','orangered','navy','magenta','black','cyan','green']

nzk = 4
nppw = 5
figure(5)
clf()

figure(6)
clf()
for i in arange(nzk):
    istart = i*nppw
    iend = (i+1)*nppw
    print(istart,iend)
    figure(5)
    semilogy(noverk[istart:iend],errd[istart:iend],'.',color=cc[i])
    figure(6)
    plot(noverk[istart:iend],1.0/errq[istart:iend],'.',color=cc[i])
    
    
figure(5)    
#legend(['k=10','k=20','k=40','k=80'])
legend(['k=10','k=40','k=160','k=640'])
xlabel('N/k')
ylabel(r'$|\sigma|$')
savefig('err_dens_neu_pndat.pdf',bbox_inches='tight')
title('Neumann')
figure(6)    
xlabel('N/k')
legend(['k=10','k=40','k=160','k=640'])
#legend(['k=10','k=20','k=40','k=80'])
ylabel(r'$|\sigma-\sigma_{n}|/|(I-P_{N})\sigma|$')
title('Neumann')
savefig('err_cq_neu_pndat.pdf',bbox_inches='tight')

show()
