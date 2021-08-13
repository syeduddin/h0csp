import sys
import numpy as np
from numpy import matrix
import emcee
import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import triangle
import random,os
from astropy.cosmology import FlatLambdaCDM
from scipy import optimize
from multiprocessing import Pool
from multiprocessing import cpu_count
import time
from astropy.io import ascii





c = 300000. # km/sec
q=-0.59


#tab = ascii.read('Bmax+MassCSPI.dat')


#initial guess
pl.figure(figsize=(20,10))

filter = ['u','B','g','V','r','i','Y','J','H']
p0 = [-18.67,-19.09,-19.08,-19.08,-19.03,-18.44,-18.46,-18.6,-18.36]
p1 = [-1.72,-1.19,-1.11,-1.21,-1.13,-0.73,-0.14,-0.64,-0.28]
p2 = [-1.66,-1.3,-1.03,-1.5,-1.14,-0.06,1.58,0.007,1.03]
rv = [3.25,2.06,1.54,1.1,0.65,0.14,-0.53,-0.59,-0.81]
alpha = [-0.06,-0.02,-0.02,-0.02,-0.02,-0.02,-0.03,-0.03,-0.03]
sig = [0.22,.17,0.13,0.17,0.16,0.15,0.15,0.22,0.20]
vel = [471,441,477,457,434,412,342,359,319]
h0 =[74.25,72.89,74.84,72.9,72.2,72.7,75.28,74.19,75.69]          

# Eqn 9 of Bruns 2018
def distmod(h,z1,z2):
    t1 = (1+z1)/(1+z2)
    t2 = (c*z2)/h
    t3 = 1 + ((1-q)*z1/2)
    return (5*np.log10(t1*t2*t3)) +25

# Liklihood function
for i in range(len(filter)):
    tab = ascii.read('data/'+filter[i]+'_ceph.dat')

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
    

    st1 = p1[i]*(st - 1)
    st2 = p2[i]*((st - 1)**2)
    red = rv[i]*(bv)
    mu_obs = mmax - p0[i] - st1 - st2 - red - alpha[i]*(m_csp-np.median(m_csp)) 
    mu_model = distmod(h0[i],zhel,zcmb)
    fac= (p1[i]+(2*p2[i]*st))
    err = ((fac*est)**2) +(emmax**2) +((rv[i]*ebv)**2)+(2*rv[i]*c_ms)+(2*fac*c_mbv)+(sig[i]**2)+(((2.17*vel[i])/(zcmb*c))**2)

    
    dmu=mu_obs-mu_model

   
    err_int=(np.std(dmu)) - (np.mean(err))
    #wt= np.sqrt(1/((ey**2)+(err_int**2))) # weights
    wt= 1/((err**2)+(err_int**2)) # weights

    mean= np.sum(dmu*wt)/np.sum(wt)
    error= np.sqrt((1/np.sum(wt)))
    
    print filter[i], '%0.3f'%mean,'%0.3f'%error, '%0.3f'%np.std(dmu), '%0.3f'%(np.std(dmu)/np.sqrt(len(dmu)))


    
    pl.subplot(3,3,i+1)
   
    pl.grid()
    Ho_dists = tab['dist'] > 0
   
    pl.errorbar(zcmb,dmu,yerr=err,fmt='o',color='0.65',label='$'+filter[i]+'$')
    #pl.errorbar(zcmb[Ho_dists],dmu[Ho_dists],yerr=err[Ho_dists],fmt='o',color='b')
    z = np.arange(0.0001,.15,.001)
    evel= (((2.17*vel[i])/(z*c))**2) +(sig[i]**2)

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
    pl.ylim(-1.5,1.5),pl.xlim(0.0,0.14)
    pl.axhline(0,color='k')
    pl.legend(numpoints =1)
    pl.ylabel(r'$\Delta \mu \ (mag)$' ,fontsize=18)
    pl.xlabel(r'$z_{cmb}$' ,fontsize=18)
    #pl.axhline(sig[i],color='g')
    #pl.axhline(-sig[i],color='g')
    #pl.savefig('plots/hd_burns.pdf')
pl.tight_layout()
pl.savefig('plots/hd_ceph.pdf')
#pl.show()






