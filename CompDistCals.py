from astropy.io import ascii
import matplotlib.pylab as pl
import numpy as np
from scipy.stats import norm

tab1 = ascii.read('data/final/calibrators_sbf.csv')
tab2 = ascii.read('data/final/calibrators_trgb.csv')
tab3 = ascii.read('data/final/calibrators_cepheids.csv')



print (np.mean(tab1['m']))
print (np.mean(tab2['m']))
print (np.mean(tab3['m']))

sys.exit()


w = np.where(tab1['dist']>32.0)
#print tab1['sn'][w]

xrange=[np.min(tab1['dist']),np.max(tab1['dist'])]



pl.hist(tab3['dist'],bins=10,histtype='stepfilled',alpha=.5,label='Cepheid',color='r')

pl.hist(tab2['dist'],bins=10,histtype='stepfilled',alpha=.5,label='TRGB',color='g')
pl.hist(tab1['dist'],bins=10,histtype='stepfilled',alpha=.5,label='SBF',color='b')
pl.xlabel(r'$Distacne \ moduli \ (\mu)$',fontsize=20), pl.ylabel('$Number$',fontsize=20)
pl.ylim(0,10),pl.xlim(30,35)
pl.legend(),pl.grid()
pl.savefig('plots/caldist.pdf')
pl.show()


