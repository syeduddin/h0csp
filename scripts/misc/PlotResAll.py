import sys
import numpy as np
import emcee
import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import random,os
from astropy.io import ascii
from astropy.table import join


B = ascii.read('../results/Ceph_res_B.csv')

filter = ['u','B','g','V','r','i','Y','J','H']

for i in range(len(filter)):
    data =  ascii.read('../results/Ceph_res_'+filter[i]+'.csv')
    tab = join(B,data,keys='sn')
    w = np.where(tab['zcmb_1']>0.0)
    bres =tab['res_1'][w]
    ares =tab['res_2'][w]
    ebres =tab['eres_1'][w]
    eares =tab['eres_2'][w]
    sn = tab['sn'][w]
    w1= np.where(ares<-1.3)
    print(sn[w1],ares[w1])
    #sys.exit()
    pl.subplot(3,3,i+1)
    pl.errorbar(bres,ares,yerr=eares,xerr=ebres,fmt='.',mfc='white',color='k',label='$'+filter[i]+'$')
    pl.legend(numpoints =1)
pl.tight_layout()
pl.savefig('../plots/compres.pdf')
