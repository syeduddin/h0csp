from astropy.io import ascii
from astropy.table import join,vstack,Table
import sys
import numpy as np
import matplotlib.pylab as pl



ceph = ascii.read('../../data/working/B_ceph_update3.csv')
trgb = ascii.read('../../data/working/B_ceph_update3.csv')
sbfk21 = ascii.read('../../data/working/B_ceph_update3.csv')
sbfj21 = ascii.read('../../data/working/B_ceph_update3.csv')

cals = ['ceph','trgb','sbfk21','sbfj21']

for i in range(len(cals)):
    
    print (cals['sn'][i])
    #pl.errorbar(cals[i]['zcmb'],cals[i]['dist'],yerr=cals[i]['edist'])

