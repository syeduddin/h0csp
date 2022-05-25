from astropy.io import ascii
from astropy.table import join,vstack, Table, unique, hstack
import numpy as np


t1 = ascii.read('../../data/working/BigTableH.csv',format='csv') # distant sample
#t2 = ascii.read('../../data/working/H_trgb.csv') # distant sample
#t3 = ascii.read('../../data/working/H_sbfcombined.csv') # distant sample


w = np.where(t1['dist']<1.0)
#t1 = unique(t1[w],keys='sn')

wc = np.where((t1['dist']>1.0) &(t1['caltype']=='c'))
#tc = unique(t1[wc],keys='sn')

wt = np.where((t2['dist']>1.0) &(t2['caltype']=='t'))
#tt = unique(t2[wt],keys='sn')

ws = np.where((t3['dist']>1.0) &(t3['caltype']=='s'))
#ts = unique(t3[ws],keys='sn')

print (len(t1[w]),len(t1[wc]),len(t2[wt]),len(t3[ws]))
table = vstack([t1[w],t1[wc],t2[wt],t3[ws]])

table = unique(table,keys='sn')

table.write('../../data/working/H_all.csv', format='ascii.csv', delimiter=',',overwrite=True)



print (table)
