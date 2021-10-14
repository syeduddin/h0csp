from astropy.io import ascii
import matplotlib.pylab as pl
import numpy as np
from scipy.stats import norm
from astropy.table import join

tab0 = ascii.read('../data/calibrators/calibrators_sbf.csv')

tab1 = ascii.read('../data/calibrators/calibrators_sbf_new.csv')
tab2 = ascii.read('../data/calibrators/calibrators_trgb.csv')
tab3 = ascii.read('../data/calibrators/calibrators_cepheids.csv')


t = join(tab0,tab1,keys='sn')

(t['sn', 'dist_1','dist_2','host_1']).pprint_all()




print (np.mean(tab1['m']))
print (np.mean(tab2['m']))
print (np.mean(tab3['m']))



w = np.where(tab1['dist']>32.0)
#print tab1['sn'][w]

xrange=[np.min(tab1['dist']),np.max(tab1['dist'])]



pl.hist(tab3['dist'],bins=10,range=[30,34],histtype='step',label='Cepheid',color='b',lw=2)
pl.hist(tab2['dist'],bins=10,range=[30,34],histtype='step',label='TRGB',color='r',lw=2)
pl.hist(tab1['dist'],bins=10,range=[30,34],histtype='step',label='SBF',color='#377eb8',lw=2)

pl.xlabel(r'$Distacne \ moduli \ (\mu)$',fontsize=14), pl.ylabel('$Number$',fontsize=14)
pl.ylim(0,10),pl.xlim(30,35)
pl.legend(),pl.grid()
pl.savefig('../plots/caldist.pdf')
#pl.show()


