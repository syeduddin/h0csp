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
for i in range(len(filter)):
    
    result = ascii.read('../../results/'+filter[i]+'_ceph_result.txt')
    #result = ascii.read('../../results/B_ceph_result.txt')
    p0=result['p0'][0]
    p1=result['p1'][0]
    p2=result['p2'][0]
    alpha=result['alpha'][0]
    beta=result['beta'][0]
    sig=result['sig_int'][0]
    vel=result['vel'][0]
    h0=result['H0'][0]
    
    
    tab = ascii.read('../../data/working/'+filter[i]+'_ceph.csv')
    #w= np.where((tab['subtype']!='Ia-91T') & (tab['subtype']!='Ia-91bg') & (tab['caltype']=='none')& (tab['sample']!='bla'))
   
    w = np.where((tab['sn']!='CSP14abk') &  (tab['sn']!='PTF13dyt') &  (tab['sn']!='PTF13dym') &  (tab['sn']!='PS1-13eao'))

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
    mu_model = np.where(Ho_dist,distmod(h0,zhel,zcmb), dist)
    fac= (p1+(2*p2*st))

    err = ((fac*est)**2) +(emmax**2) +((beta*ebv)**2)+(2*beta*c_ms)+(2*fac*c_mbv)+(sig**2)+(((2.17*vel)/(zcmb*c))**2)

    err2 = ((fac*est)**2) +(emmax**2) +((beta*ebv)**2)+(2*fac*c_ms)+(beta*c_mbv)+(edist**2)
    err = np.where(Ho_dist,err,err2)
    dmu=mu_obs-mu_model
    w1 = np.where(np.abs(dmu)>1.)

    #print (filter[i], sn[w1], dmu[w1])
    
    #sys.exit()
    rms=np.sqrt(np.sum(np.power(dmu,2))/len(dmu))
    #err_int=(np.std(dmu)) - (np.mean(err))
    #wt= np.sqrt(1/((ey**2)+(err_int**2))) # weights
    #wt= 1/((err**2)+(err_int**2)) # weights

    #mean= np.sum(dmu*wt)/np.sum(wt)
    #error= np.sqrt((1/np.sum(wt)))
    
    #print (filter[i], '%0.3f'%mean,'%0.3f'%error, '%0.3f'%np.std(dmu), '%0.3f'%(np.std(dmu)/np.sqrt(len(dmu))))
    print (filter[i],'Full','%0.3f'%rms)
    #sys.exit()
    data=Table()
    data['sn']=sn
    data['res']=dmu
    data['eres']=err
    data['zcmb']=zcmb
    data['st']=st
    data['B-V']=bv
    data['m'] =m_csp
    data['ml'] =ml
    data['mu'] =mu
    data['sample']=sample
    data['cal']=cal
    
    ascii.write(data,'../../results/Ceph_res_'+filter[i]+'.csv',format='csv', delimiter=',',overwrite=True)
    
    pl.subplot(3,3,i+1)
   
    pl.grid()
    Ho_dists = tab['dist'][w] > 0
   
    pl.errorbar(zcmb,dmu,yerr=err,fmt='o',mfc='white',color='k',ms=12,label='$'+filter[i]+'$')
    wt= np.where((subtype=='Ia-91T'))
    wbg= np.where((subtype=='Ia-91bg'))
   
    
    pl.errorbar(zcmb[wbg],dmu[wbg],yerr=err[wbg],fmt='s',color='#f781bf',ms=10,markeredgecolor='k')
    pl.errorbar(zcmb[wt],dmu[wt],yerr=err[wt],fmt='p',color='b',ms=12,markeredgecolor='k')
    pl.errorbar(zcmb[Ho_dists],dmu[Ho_dists],yerr=err2[Ho_dists],fmt='d',color='g',ms=12,markeredgecolor='k')

    rms91t=np.sqrt(np.sum(np.power(dmu[wt],2))/len(dmu[wt]))
    rms91bg=np.sqrt(np.sum(np.power(dmu[wbg],2))/len(dmu[wbg]))
    print (filter[i],'91T','%0.3f'%rms91t, len(dmu[wt]))
    print (filter[i],'91bg','%0.3f'%rms91bg, len(dmu[wt]))

    z = np.arange(0.0001,.15,.001)
    evel= (((2.17*vel)/(z*c))**2) +(sig**2)

    evel = np.sqrt(evel)

    from operator import itemgetter
    pick_0 = itemgetter(0)
    pick_1 = itemgetter(1)
    x_decor = sorted(enumerate(z), key=pick_1)
    x_idxs = map(pick_0, x_decor)
    multi_picker = itemgetter(*x_idxs)
    multi_picker(evel)


    pl.plot(multi_picker(z),multi_picker(evel),'r-',lw=3,zorder=3)

    pl.plot(multi_picker(z),multi_picker(-evel),'r-',lw=3,zorder=3)
    pl.ylim(-1.5,1.5),pl.xlim(0.0,0.15)
    pl.axhline(0,color='k')
    pl.legend(numpoints =1,fontsize=16)
    pl.ylabel(r'$\Delta \mu \ (mag)$' ,fontsize=18)
    pl.xlabel(r'$z_{cmb}$' ,fontsize=18)
    #pl.axhline(sig[i],color='g')
    #pl.axhline(-sig[i],color='g')
    #pl.savefig('plots/hd_burns.pdf')
pl.tight_layout()
pl.savefig('../../plots/hd_ceph.pdf')
#pl.show()






