import sys
import numpy as np
from numpy import matrix
import emcee
import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import triangle
import random,os
from astropy.cosmology import FlatLambdaCDM
import pandas as pd
from scipy import optimize
from multiprocessing import Pool
from multiprocessing import cpu_count
import time
from astropy.io import ascii



c = 300000. # km/sec
q=-0.59


tab = ascii.read('data/B_ceph.dat')

Ho_dists = tab['dist'] < 0
st = tab['st']
est = tab['est']
zhel = tab['zhel']
zcmb = tab['zcmb']
mmax = tab['Mmax']
emmax = tab['eMmax']
bv = tab['BV']
ebv = tab['eBV']
m_csp = tab['m']
dist = tab['dist']
edist = tab['edist']
c_ms = tab['covMs']
c_mbv = tab['covBV_M']


#initial guess
p = -19.088
p1 = -1.19
p2 = -1.3
rv = 3.05
alpha = -0.02
sig = 0.17
vel = 441.0

# Eqn 9 of Bruns 2018
def distmod(h,z1,z2):
    t1 = (1+z1)/(1+z2)
    t2 = (c*z2)/h
    t3 = 1 + ((1-q)*z1/2)
    return (5*np.log10(t1*t2*t3)) +25

# Liklihood function


st1 = p1*(st - 1)
st2 = p2*((st - 1)**2)
red = rv*(bv)
mu_obs = mmax - p - st1 - st2 - red - alpha*(m_csp-np.median(m_csp))

mu_model = np.where(Ho_dists,distmod(72.89,zhel,zcmb), dist)
fac= (p1+(2*p2*st))

err = (fac*est)**2 +emmax**2 +(rv*ebv)**2+2*fac*c_ms+rv*c_mbv+sig**2+(0.00000723*vel/zcmb)**2
err1 = ((fac*est)**2) +(emmax**2) +((rv*ebv)**2)+(2*fac*c_ms)+(rv*c_mbv)+(edist**2)
mu_stat = np.where(Ho_dists,err,err1)
mu_stat=np.sqrt(mu_stat)
dmu=mu_obs-mu_model

w = (dmu<0) 

print len(dmu[w])
#sys.exit()




pl.errorbar(zcmb[Ho_dists],dmu[Ho_dists],yerr=err[Ho_dists],fmt='o',color='.65',label='$CSP \ SNe  \ Ia$')
Ho_dists = tab['dist'] > 0
pl.errorbar(zcmb[Ho_dists],dmu[Ho_dists],yerr=err[Ho_dists],fmt='o',color='b',label='$Cepheids$')
pl.legend(numpoints=1)
z = np.arange(0,.14,.001)
evel= (((2.17*vel)/(z*c))**2) +(.18**2)

evel = np.sqrt(evel)

from operator import itemgetter
pick_0 = itemgetter(0)
pick_1 = itemgetter(1)
x_decor = sorted(enumerate(z), key=pick_1)
x_idxs = map(pick_0, x_decor)
multi_picker = itemgetter(*x_idxs)
multi_picker(evel)


pl.plot(multi_picker(z),multi_picker(evel),'r-',lw=2)

pl.plot(multi_picker(z),multi_picker(-evel),'r-',lw=2)
pl.ylim(-2,2)
pl.axhline(0,color='k')
pl.ylabel(r'$\Delta \mu \ (mag)$' ,fontsize=18)
pl.xlabel(r'$z_{cmb}$' ,fontsize=18)
#pl.axhline(0.18,color='g')
#pl.axhline(-0.18,color='g')
pl.savefig('plots/hd_Bceph.pdf')
pl.show()






