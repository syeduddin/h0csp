from astropy.io import ascii
import numpy as np
import matplotlib.pylab as pl
import sys
from astropy.table import join

trgb = ascii.read('../../data/working/B_trgb.csv')
ceph = ascii.read('../../data/working/B_ceph.csv')

w1 = np.where(trgb['sample']=='CSPI')
w2 = np.where(trgb['sample']=='CSPII')

print (trgb.colnames)
param = sys.argv[1]

w = np.where ((trgb['dist']>0))



pl.figure(1)
pl.hist(trgb[param][w1],histtype='step',lw=2,label='CSPI')
pl.hist(trgb[param][w2],histtype='step',lw=2,label='CSPII')
#pl.hist(trgb[param][w],histtype='step',lw=2,label='TRGB')
#pl.hist(ceph[param][w],histtype='step',lw=2,label='Ceph')

pl.legend(loc='upper left')
pl.savefig('../../plots/comp'+param+'.pdf')



pl.figure(2)

b = ascii.read('../../data/working/B_trgb.csv')
h = ascii.read('../../data/working/H_trgb.csv')

t = join(b,h,keys='sn')
wb = np.where (b['dist']>0)
wh = np.where ((h['dist']>0))

#wt = np.where ((t['dist_1']<0))

print (h[param][wh])

pl.hist(b[param][wb],histtype='step',lw=2,label='B')
pl.hist(h[param][wh],histtype='step',lw=2,label='H')
pl.xlabel(param)
pl.legend(loc='upper left')
pl.savefig('../../plots/compBH'+param+'.pdf')

