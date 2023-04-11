from numpy import *
from pylab import *
x = loadtxt('fort.37')
y = loadtxt('fort.38')
(m,n) = shape(x)
muse = int(m/4)
figure(1)
colors = ['navy','orangered','darkgreen','cyan']
for i in range(4):
    istart = i*muse
    iend = (i+1)*muse
    plot(x[istart:iend,0],x[istart:iend,1],'.',color=colors[i])
    h = 0.1
    xuse = x[istart:iend,0] + h*x[istart:iend,2]
    yuse = x[istart:iend,1] + h*x[istart:iend,3]
    plot(xuse,yuse,'*',color=colors[i])
plot(y[0],y[1],'r.')
axis('equal')
show()
