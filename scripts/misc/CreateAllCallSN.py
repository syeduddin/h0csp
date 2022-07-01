from astropy.io import ascii
from astropy.table import join,vstack, Table, unique, hstack
import numpy as np


t1 = ascii.read('../../data/working/H_ceph.csv',format='csv') # distant sample
t2 = ascii.read('../../data/working/H_trgb.csv') # distant sample
t3 = ascii.read('../../data/working/H_sbfcombined.csv') # distant sample


w = np.where(t1['dist']<1.0)
#t1 = unique(t1[w],keys='s')

wc = np.where((t1['dist']>1.0))
#tc = unique(t1[wc],keys='sn')

tc = vstack([t1[w],t1[wc]])
#tc = unique(t1[wc],keys='sn')


wt = np.where((t2['dist']>1.0))
#tt = unique(t2[wt],keys='sn')

ws = np.where((t3['dist']>1.0) &(t3['caltype']=='s'))
#ts = unique(t3[ws],keys='sn')

print (len(t1[w]),len(t1[wc]),len(t2[wt]),len(t3[ws]))
table = vstack([tc,t2[wt],t3[ws]])

#table = unique(table,keys='sn')

table.write('../../data/working/H_all.csv', format='ascii.csv', delimiter=',',overwrite=True)



print (table)
