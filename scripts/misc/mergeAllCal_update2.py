from astropy.table import join,vstack,unique
from astropy.io import ascii
import sys
import numpy as np
filter = sys.argv[1]
tab1 = ascii.read('../../data/working/'+filter+'_ceph_update3.csv')

tab2 = ascii.read('../../data/working/'+filter+'_trgb_update3.csv')
tab3 = ascii.read('../../data/working/'+filter+'_sbfj21_update3.csv')

w2 = np.where(tab2['dist']>0)
w3 = np.where(tab3['caltype']=='s')

bigt = vstack([tab2[w2],tab3[w3]])

newt=  vstack([tab1,bigt])

print (len(newt[w2]))

newt.write('../../data/working/'+filter+'_nok21_update3.csv', format='csv', delimiter=' ',overwrite=True)
