from numpy import *
from pylab import *
from scipy.io import FortranFile
from scipy.special import *

ik = 6
nch = 2024
fname = 'cavity_data/sol_ik'+str(ik).zfill(2)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2_lin.bin'
fref_tot = 'cavity_data/fields_ref_utot_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
fref_sc = 'cavity_data/fields_ref_usc_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
fnord1_tot = 'cavity_data/fields_nord1_utot_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
fnord1_sc = 'cavity_data/fields_nord1_usc_ik'+str(ik)+'_nch'+str(nch).zfill(5)+'_pw_mpi2p0p2.pdf'
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

snorm = sqrt(dot(abs(soln_use)**2,qwts))
slen = sum(qwts)
print(snorm)
print(slen)
