from numpy import *
from pylab import *
from scipy.io import FortranFile

ik = 5
nch = 22640
fname = 'diamond_data/sol_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_ithet005.bin'
fref_tot = 'diamond_data/fields_ref_utot_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_ithet005.pdf'
fref_sc = 'diamond_data/fields_ref_usc_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_ithet005.pdf'
fnord1_tot = 'diamond_data/fields_nord1_utot_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_ithet005.pdf'
fnord1_sc = 'diamond_data/fields_nord1_usc_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_ithet005.pdf'
f = FortranFile(fname,'r')
ncomp = int(f.read_ints(int32))
rsc = f.read_reals(float)
shifts = f.read_reals(float)
zk = f.read_record(complex)
nch0 = int(f.read_ints(int32))
nch = int(f.read_ints(int32))
ithet = f.read_ints(int32)
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
isin = f.read_ints(int32)
pottarg_plot_nord1 = f.read_record(complex)
soln_use = f.read_record(complex)
pottarg_plot_ref = f.read_record(complex)
ct = cos((ithet+0.0)*2*pi/360.0)
st = sin((ithet+0.0)*2*pi/360.0)
uinc = exp(1j*zk[0]*(targs[0,:]*ct + targs[1,:]*st))

xt = targs[0,:]
yt = targs[1,:]

utot = uinc - pottarg_plot_ref
utot[isin > -4] = NaN 
utot[(xt>3.03) & (xt<3.92) & (yt>3.5) & (yt<4.6)] = NaN
uplot = abs(utot)
uplot = uplot.reshape(901,901,order="F")

figure(1)
clf
imshow(uplot,extent=[-10,10,-10,10],vmin=0,vmax=5)
cb = colorbar(shrink=.9)
savefig(fref_tot,bbox_inches='tight')


usc = - pottarg_plot_ref
usc[isin > -4] = NaN 
usc[(xt>3.03) & (xt<3.92) & (yt>3.5) & (yt<4.6)] = NaN
uplot = abs(usc)
uplot = uplot.reshape(901,901,order="F")

figure(2)
clf
imshow(uplot,extent=[-10,10,-10,10],vmin=0,vmax=5)
cb = colorbar(shrink=.9)
savefig(fref_sc,bbox_inches='tight')


utot = uinc - pottarg_plot_nord1
utot[isin > -4] = NaN 
utot[(xt>3.03) & (xt<3.92) & (yt>3.5) & (yt<4.6)] = NaN
uplot = abs(utot)
uplot = uplot.reshape(901,901,order="F")

figure(3)
clf
imshow(uplot,extent=[-10,10,-10,10],vmin=0,vmax=5)
cb = colorbar(shrink=.9)
savefig(fnord1_tot,bbox_inches='tight')


usc = - pottarg_plot_nord1
usc[isin > -4] = NaN 
usc[(xt>3.03) & (xt<3.92) & (yt>3.5) & (yt<4.6)] = NaN
uplot = abs(usc)
uplot = uplot.reshape(901,901,order="F")

figure(4)
clf
imshow(uplot,extent=[-10,10,-10,10],vmin=0,vmax=5)
cb = colorbar(shrink=.9)
savefig(fnord1_sc,bbox_inches='tight')


show()




show()



