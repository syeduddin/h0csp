from astropy.io import ascii
from astropy.table import join,vstack, Table, unique
import sys
import numpy as np
filter = sys.argv[1]


t0 = ascii.read('../../data/lc/'+filter+'_max_SBF_new.dat') # for SBF  lc data
t0 = unique(t0,keys='name')

t0.remove_column('Ho')
t0.remove_column('eHo')
t0.remove_column('method')


t = ascii.read('../../data/lc/'+filter+'_max.csv') # for SBF  lc data


print (t.colnames)

need = ['SN2008R','PTF13ebh','SN2008ia','SN2006ef']

w1 = np.where(t['name']=='SN2008R')
w2 = np.where(t['name']=='PTF13ebh')
w3 = np.where(t['name']=='SN2008ia')
w4 = np.where(t['name']=='SN2006ef')

t1 = vstack([t[w1],t[w2],t[w3],t[w4]])

t2 = vstack([t0,t1])


del t2['dist','sample', 'subtype', 'quality', 'cosmo', 'phys']

t2.write('../../data/lc/'+filter+'_max_SBFJ21.csv', format='ascii.csv', delimiter=',',overwrite=True)

print (t2)
