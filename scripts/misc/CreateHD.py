import sys
import numpy as np
from numpy import matrix
#import emcee
import astropy.io.fits as pyfits
import matplotlib.pylab as pl
#simport triangle
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

#filter=['B','H']
for i in range(len(filter)):

    result = ascii.read('../../results/'+filter[i]+'_sbfj21_update3_result.txt')
    #result = ascii.read('../../results/B_trgb_result.txt')
    p0=result['p0'][0]
    ep0 = (result['p0'][1]+result['p0'][2])/2
    p1=result['p1'][0]
    ep1 = (result['p1'][1]+result['p1'][2])/2
    p2=result['p2'][0]
    ep2 = (result['p2'][1]+result['p2'][2])/2
    alpha=result['alpha'][0]
    ealpha = (result['alpha'][1]+result['alpha'][2])/2

    rv=result['beta'][0]
    erv = (result['beta'][1]+result['beta'][2])/2
    sig=result['sig_int'][0]
    vel=result['vel'][0]
    h0=result['H0'][0]
    eh0 = (result['H0'][1]+result['H0'][2])/2

    #ep2=ep2/10.
    
    
    tab = ascii.read('../../data/working/'+filter[i]+'_sbfj21_update3.csv')
   
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
    c_sbv = tab['covBVs'][w]

    sn = tab['sn'][w]
    sample = tab['sample'][w]
    cal = tab['caltype'][w]
    host = tab['host'][w]
   
    subtype = tab['subtype'][w] 
    Ho_dist = tab['dist'][w]<0

    st1 = p1*(st - 1)
    st2 = p2*((st - 1)**2)
    red = rv*(bv)
    
    mu_obs = mmax - p0 - st1 - st2 - red - alpha*(m_csp-np.median(m_csp)) 
    absmag = p0 + st1 + st2 + red + alpha*(m_csp-np.median(m_csp)) 


    mu_model = np.where(Ho_dist,distmod(h0,zhel,zcmb), dist)
    fac= (p1+(2*p2*st))

    err1 = (emmax**2) + ((fac*est)**2) +((rv*ebv)**2) -(2*fac*c_ms)+(2*rv*fac*c_sbv) -(2*fac*c_mbv)+((alpha*em)**2) + (sig**2) + ((0.00000723*vel/zcmb)**2)
        
    err2 = (emmax**2) + ((fac*est)**2) +((rv*ebv)**2) -(2*fac*c_ms)+(2*rv*fac*c_sbv) -(2*fac*c_mbv)+((alpha*em)**2)  + (edist**2)
            
    

    err = np.where(Ho_dist,err1,err2)
    err = np.sqrt(err)
    dmu=mu_obs-mu_model
    w1 = np.where(np.abs(dmu)>1.)

    rms=np.sqrt(np.sum(np.power(dmu,2))/len(dmu))
    #err_int=(np.std(dmu)) - (np.mean(err))
    #wt= np.sqrt(1/((ey**2)+(err_int**2))) # weights
    #wt= 1/((err**2)+(err_int**2)) # weights

    #mean= np.sum(dmu*wt)/np.sum(wt)
    #error= np.sqrt((1/np.sum(wt)))
    
    #print (filter[i], '%0.3f'%mean,'%0.3f'%error, '%0.3f'%np.std(dmu), '%0.3f'%(np.std(dmu)/np.sqrt(len(dmu))))
    #print (filter[i],'Full','%0.3f'%rms)
    #sys.exit()
    data=Table()
    data['sn']=sn
    #data['res']=dmu
    
    data['mu_obs']=mu_obs.round(3)
    data['err_mu_obs']=err.round(3)
    
    data['zcmb']=zcmb.round(4)
    #data['sBV']=st
    #data['B-V']=bv
    data['mass_best'] =m_csp
    data['mass_low'] =ml
    data['mass_hhigh'] =mu
    data['sample']=sample
    #data['cal']=cal
   #s data['host']=host
    #print(data)
    #ascii.write(data,'../../results/Ceph_res_'+filter[i]+'_update2_vpec.csv',format='csv', delimiter=',',overwrite=True)
    ascii.write(data,'../../results/forBrent/resids_'+filter[i]+'_update3.txt',format='tab',overwrite=True)

    
    pl.subplot(3,3,i+1)
   
    #pl.grid()
    Ho_dists = tab['dist'][w] > 0
   
    #pl.errorbar(zcmb,dmu,yerr=err,fmt='o',mfc='white',color='k',ms=12,label='$'+filter[i]+'$')
    wt= np.where((subtype=='Ia-91T'))
    wbg= np.where((subtype=='Ia-91bg'))
    
    #pl.errorbar(zcmb[wbg],dmu[wbg],yerr=err[wbg],fmt='s',color='#f781bf',ms=10,markeredgecolor='k')
    #pl.errorbar(zcmb[wt],dmu[wt],yerr=err[wt],fmt='p',color='b',ms=12,markeredgecolor='k')
    #pl.errorbar(zcmb[Ho_dists],dmu[Ho_dists],yerr=err2[Ho_dists],fmt='d',color='g',ms=12,markeredgecolor='k')

    rms91t=np.sqrt(np.sum(np.power(dmu[wt],2))/len(dmu[wt]))
    rms91bg=np.sqrt(np.sum(np.power(dmu[wbg],2))/len(dmu[wbg]))
    #print (filter[i],'91T','%0.3f'%rms91t, len(dmu[wt]))
    #print (filter[i],'91bg','%0.3f'%rms91bg, len(dmu[wt]))

    z = np.arange(0.0001,.15,.001)
    evel= (((2.17*vel)/(z*c))**2) +(sig**2)

    evel = np.sqrt(evel)

    #from operator import itemgetter
    #pick_0 = itemgetter(0)
    #pick_1 = itemgetter(1)
    #x_decor = sorted(enumerate(z), key=pick_1)
    #x_idxs = map(pick_0, x_decor)
    #multi_picker = itemgetter(*x_idxs)
    #multi_picker(evel)


    #pl.plot(multi_picker(z),multi_picker(evel),'r-',lw=3,zorder=3)

    #pl.plot(multi_picker(z),multi_picker(-evel),'r-',lw=3,zorder=3)
    #pl.ylim(-1.5,1.5),pl.xlim(0.0,0.15)
    #pl.axhline(0,color='k')
    #pl.legend(numpoints =1,fontsize=16)
    #pl.ylabel(r'$\Delta \mu \ (mag)$' ,fontsize=18)
    #pl.xlabel(r'$z_{cmb}$' ,fontsize=18)
    #pl.axhline(sig[i],color='g')
    #pl.axhline(-sig[i],color='g')
    #pl.savefig('plots/hd_burns.pdf')
    bin_width = 0.1
    bins = np.arange(min(dmu), max(dmu) + bin_width, bin_width)
    pl.hist(dmu,histtype='step',lw=2,color='k',bins=bins)
    pl.axvline(0,lw=3)
pl.tight_layout()
pl.savefig('../../plots/hist_hr_sbfj21_uddin24.pdf')
#pl.show()






