from numpy import *
from pylab import *

x = loadtxt('fort.38')
y = loadtxt('fort.39')
z = loadtxt('fort.40')

figure(1)
plot(x[:,0],x[:,1],'k.')
plot(y[:,0],y[:,1],'r.')

figure(2)
plot(z[:,0],z[:,1],'k.')
plot(z[:,0],z[:,3],'r.')
show()
