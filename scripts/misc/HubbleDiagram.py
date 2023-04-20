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
from astropy.table import Table
c = 300000. # km/sec
q=-0.59 # decelertion parameter


# Eqn 9 of Bruns 2018
def distmod(h,z1,z2):
    t1 = (1+z1)/(1+z2)
    t2 = (c*z2)/h
    t3 = 1 + ((1-q)*z1/2)
    return (5*np.log10(t1*t2*t3)) +25

# Getting results
pl.figure(figsize=(20,10))
filter = ['u','B','g','V','r','i','Y','J','H']

#filter='B'
for i in range(len(filter)):
    
    result = ascii.read('../../results/'+filter[i]+'_ceph_update2_result.txt')
    #result = ascii.read('../../results/B_trgb_result.txt')
    p0=result['p0'][0]
    ep0 = (result['p0'][1]+result['p0'][2])/2
    p1=result['p1'][0]
    ep1 = (result['p1'][1]+result['p1'][2])/2
    p2=result['p2'][0]
    ep2 = (result['p2'][1]+result['p2'][2])/2
    alpha=result['alpha'][0]
    ealpha = (result['alpha'][1]+result['alpha'][2])/2

    beta=result['beta'][0]
    ebeta = (result['beta'][1]+result['beta'][2])/2
    sig=result['sig_int'][0]
    vel=result['vel'][0]
    h0=result['H0'][0]
    eh0 = (result['H0'][1]+result['H0'][2])/2

    #ep2=ep2/10.
    
    
    tab = ascii.read('../../data/working/'+filter[i]+'_ceph_update2.csv')
   
    w = np.where((tab['sn']!='CSP14abk') &  (tab['sn']!='PTF13dyt') &  (tab['sn']!='PTF13dym') & (tab['sn']!='PTF14yw') & (tab['sn']!='PS1-13eao') & (tab['subtype']!='Ia-SC') & (tab['subtype']!='Ia-02cx') & (tab['sn']!='LSQ14fmg')& (tab['sn']!='SN2004dt')& (tab['sn']!='SN2005gj')& (tab['sn']!='SN2005hk')& (tab['sn']!='SN2006bt')& (tab['sn']!='SN2006ot')& (tab['sn']!='SN2007so')& (tab['sn']!='SN2008ae')& (tab['sn']!='SN2008bd')& (tab['sn']!='SN2008ha')& (tab['sn']!='SN2008J')& (tab['sn']!='SN2009dc')& (tab['sn']!='SN2009J')& (tab['sn']!='SN2010ae'))

    st = tab['st'][w]
    est = tab['est'][w]
    zhel = tab['zhel'][w]
    zcmb = tab['zcmb'][w]
    mmax = tab['Mmax'][w]
    emmax = tab['eMmax'][w]
    bv = tab['BV'][w]
    ebv = tab['eBV'][w]
    m_csp = tab['m'][w]
    ml = tab['ml'][w]
    mu = tab['mu'][w]
    eml = (tab['m'][w]-tab['ml'][w])
    emu = (tab['mu'][w]-tab['m'][w])
    em = (eml+emu)/2.
    dist = tab['dist'][w]
    edist = tab['edist'][w]
    c_ms = tab['covMs'][w]
    c_mbv = tab['covBV_M'][w]
    sn = tab['sn'][w]
    sample = tab['sample'][w]
    cal = tab['caltype'][w]
   
    subtype = tab['subtype'][w] 
    Ho_dist = tab['dist'][w]<0

    st1 = p1*(st - 1)
    st2 = p2*((st - 1)**2)
    red = beta*(bv)
    
    mu_obs = mmax - p0 - st1 - st2 - red - alpha*(m_csp-np.median(m_csp)) 
    absmag = p0 + st1 + st2 + red + alpha*(m_csp-np.median(m_csp)) 


    #mu_model = np.where(Ho_dist,distmod(h0,zhel,zcmb), dist)
    fac= (p1+(2*p2*st))

    err1 = ((fac*est)**2) +(emmax**2) +((beta*ebv)**2)+(2*beta*c_ms)+(2*fac*c_mbv)+(sig**2)+(((2.17*vel)/(zcmb*c))**2)+(alpha*(em/np.log(10)*m_csp))**2 

    err2 = ((fac*est)**2) +(emmax**2) +((beta*ebv)**2)+(2*fac*c_ms)+(beta*c_mbv)+(edist**2)

    #err1 = (emmax**2)+ (ep0**2)+ ((1+(st**2))*(ep1**2)) + ((1+(4*st**2)+(st**4))*(ep2**2))                                                                
    #err2 = (emmax**2)+ (ep0**2)+ ((1+(st**2))*(ep1**2))+ ((1+(4*st**2)+(st**4))*(ep2**2))  

    err = np.where(Ho_dist,err1,err2)
    err = np.sqrt(err)
    
    
    pl.subplot(3,3,i+1)
   
    pl.grid()
    Ho_dists = tab['dist'][w] > 0
   
    pl.errorbar(np.log10(zcmb),mu_obs,yerr=err,fmt='o',mfc='white',color='k',ms=8,label='$'+filter[i]+'$')
    
    wt= np.where((subtype=='Ia-91T'))
    wbg= np.where((subtype=='Ia-91bg'))
   
    
    #pl.errorbar(np.log10(zcmb[wbg]),mu_obs[wbg],yerr=err[wbg],fmt='s',color='#f781bf',ms=10,markeredgecolor='k')
    #pl.errorbar(np.log10(zcmb[wt]),mu_obs[wt],yerr=err[wt],fmt='p',color='b',ms=12,markeredgecolor='k')
    pl.errorbar(np.log10(zcmb[Ho_dists]),mu_obs[Ho_dists],yerr=err2[Ho_dists],fmt='o',color='r',ms=12,markeredgecolor='k')

    
    z = np.arange(0.001,.15,.0001)
    evel= (((2.17*vel)/(z*c))**2) +(sig**2)

    evel = np.sqrt(evel)
    mu_model = 5.0*np.log10(c*z/h0) + 25.0
    
    pl.plot(np.log10(z),mu_model,'r-',lw=2)
    from operator import itemgetter
    pick_0 = itemgetter(0)
    pick_1 = itemgetter(1)
    x_decor = sorted(enumerate(z), key=pick_1)
    x_idxs = map(pick_0, x_decor)
    multi_picker = itemgetter(*x_idxs)
    multi_picker(evel)


    pl.plot(multi_picker(np.log10(z)),multi_picker(mu_model+evel),'r--',lw=2,zorder=3)
    pl.plot(multi_picker(np.log10(z)),multi_picker(mu_model-evel),'r--',lw=2,zorder=3)

    #pl.plot(multi_picker(z),multi_picker(-evel),'r-',lw=3,zorder=3)
    #pl.ylim(-1.5,1.5),pl.xlim(0.0,0.15)
    #pl.axhline(0,color='k')
    #pl.legend(numpoints =1,fontsize=16)
    pl.ylabel(r'$ \mu \ (' +filter[i]+') \ (mag)$' ,fontsize=18)
    pl.xlabel(r'$log10(z_{cmb})$' ,fontsize=18)
    #pl.axhline(sig[i],color='g')
    #pl.axhline(-sig[i],color='g')
    #pl.savefig('plots/hd_burns.pdf')
pl.tight_layout()
pl.savefig('../../plots/hd_update2.pdf')
#pl.show()






