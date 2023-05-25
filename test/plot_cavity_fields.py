from numpy import *
from pylab import *
from scipy.io import FortranFile
from scipy.special import *

ik = 6
nch = 3269
fname = 'cavity_data/sol_ik'+str(ik).zfill(1)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.bin'
nch = 1635
fname = 'cavity_data/sol_ik'+str(ik).zfill(2)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2_lin.bin'
fref_tot = 'cavity_data/fields_ref_utot_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
fref_sc = 'cavity_data/fields_ref_usc_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
fnord1_tot = 'cavity_data/fields_nord1_utot_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
fnord1_sc = 'cavity_data/fields_nord1_usc_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
ferr = 'cavity_data/fields_error_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
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

fname2 = 'cavity_data/cavity_interior.bin'
g = FortranFile(fname2,'r')
ntlat = int(g.read_ints(int32))
ntarg = int(g.read_ints(int32))
targs = g.read_reals(float).reshape((2,ntarg),order="F")
isin = g.read_ints(int32)
qwts = g.read_reals(float)

print(ntlat)
x = loadtxt('fort.37')

s0 = 1.0/sqrt(abs(zk))
thet = 0.45*pi
ttuse = 3.91322070806951

thet = -pi/2 + 0.2 

uinc = zeros(size(targs[0,:]))
ct = cos(thet)
st = sin(thet)
uinc = exp(1j*zk[0]*(targs[0,:]*ct + targs[1,:]*st))

xt = targs[0,:]
yt = targs[1,:]


usc = uinc - pottarg_plot_ref
usc0 = usc
uplot = abs(usc)
uplot[isin>-0.5] = NaN

umax = max(uplot)
uplot = uplot.reshape(ntlat,ntlat,order="F")

figure(1)
clf
imshow(uplot,extent=[-4,4,-4,4],vmin=0,vmax=80)
#plot(x[:,0],x[:,1],'w.')
cb = colorbar(shrink=.9)
savefig(fref_sc,bbox_inches='tight')


usc = uinc - pottarg_plot_nord1
usc1 = usc
uplot = abs(usc)
uplot[isin>-0.5] = NaN
uplot = uplot.reshape(ntlat,ntlat,order="F")

figure(2)
clf
imshow(uplot,extent=[-4,4,-4,4],vmin=0,vmax=80)
#plot(x[:,0],x[:,1],'w.')
cb = colorbar(shrink=.9)
savefig(fnord1_sc,bbox_inches='tight')


figure(3)
clf
uplot = log10(abs(usc0-usc1)/umax)
uerrmax = max(uplot)
print(uerrmax)
uplot[isin>-0.5] = NaN
uplot = uplot.reshape(ntlat,ntlat,order="F")

imshow(uplot,extent=[-4,4,-4,4],vmin=-10,vmax=0)
#plot(x[:,0],x[:,1],'w.')
cb = colorbar(shrink=.9)
savefig(fnord1_sc,bbox_inches='tight')


show()



