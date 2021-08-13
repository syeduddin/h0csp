from astropy.table import join,vstack,unique
from astropy.io import ascii
import matplotlib.pylab as pl
import numpy as np

tab1 = ascii.read('SBF.dat')
tab1.rename_column('SN_name','sn')      # for consistency
tab2 = ascii.read('../final/B_sbf.dat')



# Now we join tab1 to masses using 'sn' as the join key
t1 = join(tab1,tab2, keys='sn')

w = np.where(t1['dist']>0)

print len(t1[w])
pl.plot(t1['Bmax'][w],t1['Mmax'][w],'ro')
pl.show()

