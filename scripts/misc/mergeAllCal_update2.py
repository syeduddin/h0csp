from astropy.table import join,vstack,unique
from astropy.io import ascii
import sys
import numpy as np
filter = sys.argv[1]
tab1 = ascii.read('../../data/working/'+filter+'_ceph_update2.csv')

tab2 = ascii.read('../../data/working/'+filter+'_trgb_update2.csv')
tab3 = ascii.read('../../data/working/'+filter+'_sbfcombined_update2.csv')

w2 = np.where(tab2['dist']>0)
w3 = np.where(tab3['dist']>0)

bigt = vstack([tab2[w2],tab3[w3]])

newt=  vstack([tab1,bigt])


newt.write('../../data/'+filter+'_all_update2.csv', format='csv', delimiter=' ',overwrite=True)
