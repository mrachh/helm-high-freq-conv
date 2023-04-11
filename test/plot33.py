#!/usr/bin/python
import matplotlib.pyplot as pt
import numpy as np
import colormaps as cmaps

pt.rc("font", size=16)
x = np.loadtxt("plot33.dat1")
x = x.reshape(  601,  601, order="F")
y = np.loadtxt("plot33.dat2")
f=open("plot33.dat3","r")
iss = 0
nd =     4
#for i in range(nd):
#    nn = int(f.readline())
#    pt.fill(y[iss:iss+nn,0],y[iss:iss+nn,1],facecolor="w")
#    iss=iss+nn

pt.imshow(x, extent=[-.10000E+02,0.10000E+02,-.10000E+02,0.10000E+02],cmap = cmaps.parula)
cb=pt.colorbar(shrink=.9)
pt.savefig("plot33.pdf")
