from astropy.io import ascii
from astropy.table import join
import matplotlib.pylab as pl
import os
import numpy as np

t1 = ascii.read('../../results/Ceph_res_B.csv')
t2 = ascii.read('../../results/Ceph_res_J.csv')
t3 = ascii.read('../../results/Ceph_res_r.csv')

t = join(t1,t3,keys='sn')

pl.plot(t['res_1'],t['res_2'],'ko')
pl.savefig('../../plots/compBr.pdf')

sys.exit()

pl.subplot(2,2,1)
pl.plot(t1['B-V'],t1['res'],'ko',alpha=.3)
pl.xlabel('B-V'), pl.ylabel('HR (B)')
w = np.where(t1['m']<9.0)
pl.plot(t1['B-V'][w],t1['res'][w],'bd')

pl.subplot(2,2,2)
pl.plot(t2['B-V'],t2['res'],'ko',alpha=.3)
pl.xlabel('B-V'), pl.ylabel('HR (J)')
w = np.where(t2['m']<9.0)
pl.plot(t2['B-V'][w],t2['res'][w],'bd')
pl.subplot(2,2,3)
pl.plot(t3['B-V'],t3['res'],'ko',alpha=.3)
pl.xlabel('B-V'), pl.ylabel('HR (H)')
w = np.where(t3['m']<9.0)
pl.plot(t3['B-V'][w],t3['res'][w],'bd')
pl.tight_layout()
pl.savefig('../../plots/compbvBJH_trgb.pdf')
os.system('open ../../plots/compbvBJH_trgb.pdf')

