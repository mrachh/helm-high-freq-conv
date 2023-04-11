#!/usr/bin/python
import matplotlib.pyplot as pt
import numpy as np
import colormaps as cmaps

pt.rc("font", size=16)
x = np.loadtxt("plot35.dat1")
x = x.reshape( 3001, 4001, order="F")
y = np.loadtxt("plot35.dat2")
f=open("plot35.dat3","r")
iss = 0
nd =     1
for i in range(nd):
    nn = int(f.readline())
    pt.fill(y[iss:iss+nn,0],y[iss:iss+nn,1],facecolor="w")
    iss=iss+nn

pt.imshow(x, extent=[-.40000E+01,0.40000E+01,-.30000E+01,0.30000E+01],cmap = cmaps.parula)
cb=pt.colorbar(shrink=.9)
pt.savefig("plot35.pdf")
