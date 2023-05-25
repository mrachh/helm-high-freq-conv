from numpy import *
from pylab import *
x = loadtxt('fort.37')
plot(x[:,0],x[:,1],'k.')
axis('equal')
show()
