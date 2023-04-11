from numpy import *
from pylab import *
from scipy.io import FortranFile
from scipy.special import *

ik = 3
nch = 1130
fname = 'cavity_data/sol_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'.bin'
fref_tot = 'cavity_data/fields_ref_utot_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'.pdf'
fref_sc = 'cavity_data/fields_ref_usc_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'.pdf'
fnord1_tot = 'cavity_data/fields_nord1_utot_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'.pdf'
fnord1_sc = 'cavity_data/fields_nord1_usc_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'.pdf'
f = FortranFile(fname,'r')
zk = f.read_record(complex)
nch = int(f.read_ints(int32))
k = int(f.read_ints(int32))
kg = int(f.read_ints(int32))
nover = int(f.read_ints(int32))
sigmacoefs = f.read_record(complex)
sigmacoefs_full = f.read_record(complex)
solncoefs = f.read_record(complex)
solncoefs_full = f.read_record(complex)
erra = f.read_reals(float)
errq = f.read_reals(float)
ntlat = int(f.read_ints(int32))
ntarg = int(f.read_ints(int32))
targs = f.read_reals(float).reshape((2,ntarg),order="F")
pottarg_plot_nord1 = f.read_record(complex)
soln_use = f.read_record(complex)
pottarg_plot_ref = f.read_record(complex)

x = loadtxt('fort.37')

s0 = 1.0/sqrt(abs(zk))
thet = 0.45*pi
ttuse = 3.91322070806951

uinc = zeros(size(targs[0,:]))

xt = targs[0,:]
yt = targs[1,:]


usc = - pottarg_plot_ref
uplot = abs(usc)
uplot = uplot.reshape(ntlat,ntlat,order="F")

figure(1)
clf
imshow(uplot,extent=[-4,4,-4,4],vmin=0,vmax=2)
plot(x[:,0],x[:,1],'w.')
cb = colorbar(shrink=.9)
savefig(fref_sc,bbox_inches='tight')


usc = - pottarg_plot_nord1
uplot = abs(usc)
uplot = uplot.reshape(ntlat,ntlat,order="F")

figure(2)
clf
imshow(uplot,extent=[-4,4,-4,4],vmin=0,vmax=2)
plot(x[:,0],x[:,1],'w.')
cb = colorbar(shrink=.9)
savefig(fnord1_sc,bbox_inches='tight')

show()



