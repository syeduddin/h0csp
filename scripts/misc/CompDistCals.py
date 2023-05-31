from astropy.io import ascii
import matplotlib.pylab as pl
import numpy as np
from scipy.stats import norm
from astropy.table import join, hstack

#tab0 = ascii.read('../../data/calibrators/calibrators_B_sbfj21.csv')
tab1 = ascii.read('../../data/calibrators/calibrators_B_sbfcombined.csv')
tab2 = ascii.read('../../data/calibrators/calibrators_trgb.csv')
tab3 = ascii.read('../../data/calibrators/calibrators_ceph.csv')


t = join(tab2,tab3,keys='sn')

print (hstack([t['sn'], t['dist_1'],t['edist_1']]))
          
#(t['sn', 'dist_1','dist_2','host_1']).pprint_all()




print (np.max(tab1['dist']))
print (np.max(tab2['dist']))
print (np.max(tab3['dist']))



w = np.where(tab1['dist']>32.0)
#print tab1['sn'][w]

xrange=[np.min(tab1['dist']),np.max(tab1['dist'])]



pl.hist(tab3['dist'],bins=10,range=[30,34],histtype='step',label='Cepheid',color='b',lw=2)
pl.hist(tab2['dist'],bins=10,range=[30,34],histtype='step',label='TRGB',color='r',lw=2)
pl.hist(tab1['dist'],bins=10,range=[30,34],histtype='stepfilled',label='SBF',color='k',lw=2, alpha=.3)
#pl.hist(tab0['dist'],bins=10,range=[30,34],histtype='stepfilled',color='k',lw=2, alpha=.3)

pl.xlabel(r'$Distance \ moduli \ (\mu)$ [mag]',fontsize=14), pl.ylabel('$Number$',fontsize=14)
pl.ylim(0,8),pl.xlim(29.9,34.1)
pl.legend(loc='upper left'),pl.grid()
pl.savefig('../../plots/caldist.pdf')
#pl.show()


