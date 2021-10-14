from astropy.io import ascii
import numpy as np



sbf = ascii.read('data/final/B_sbf.dat')
ceph = ascii.read('data/final/B_ceph.dat')
trgb = ascii.read('data/final/B_trgb.dat')

bla = ascii.read('data/final/calibrators_sbf.csv')
ws = np.where(sbf['dist']>0)
wc = np.where(ceph['dist']<0)
wt = np.where(trgb['dist']<0)

#print sbf['st']

#print (len(sbf['sn'][ws]),len(ceph['sn'][wc]),len(trgb['sn'][wt]))

print len((set(sbf['sn'][ws]).difference(ceph['sn'])))
res = ascii.read('results/B_ceph_result.txt')

#print (res['p0'][0])
