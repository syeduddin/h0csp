import sys
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import triangle
import random,os
from astropy.cosmology import FlatLambdaCDM
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
#pl.figure(figsize=(20,10))
band=sys.argv[1]

file=['../../results/'+band+'_trgb_resultcsp1.txt','../../results/'+band+'_trgb_resultcsp2.txt']

lab = ['CSPI','CSPII']
lab1 = ['CSPI','CSPII']

for i in range(len(file)):
    print (file[i])
    result = ascii.read(file[i])
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


    
    
    tab = ascii.read('../../data/working/'+band+'_trgb.csv')
    #w= np.where((tab['subtype']!='Ia-91T') & (tab['subtype']!='Ia-91bg') & (tab['caltype']=='none')& (tab['sample']!='bla'))
   
   # Excluding peculiar events
    w = np.where((tab['sn']!='CSP14abk') &  (tab['sn']!='PTF13dyt') &  (tab['sn']!='PTF13dym') & (tab['sn']!='PTF14yw') & (tab['sn']!='PS1-13eao') & (tab['subtype']!='Ia-SC') & (tab['subtype']!='Ia-02cx') & (tab['sn']!='LSQ14fmg')& (tab['sn']!='SN2004dt')& (tab['sn']!='SN2005gj')& (tab['sn']!='SN2005hk')& (tab['sn']!='SN2006bt')& (tab['sn']!='SN2006ot')& (tab['sn']!='SN2007so')& (tab['sn']!='SN2008ae')& (tab['sn']!='SN2008bd')& (tab['sn']!='SN2008ha')& (tab['sn']!='SN2008J')& (tab['sn']!='SN2009dc')& (tab['sn']!='SN2009J')& (tab['sn']!='SN2010ae'))

    #  
    st = tab['st'][w]
    est = tab['est'][w]
    zhel = tab['zhel'][w]
    zcmb = tab['zcmb'][w]
    mmax = tab['Mmax'][w]
    emmax = tab['eMmax'][w]
    bv = tab['BV'][w]
    ebv = tab['eBV'][w]
    ebvmw = tab['EBVmw'][w]
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
    #sample = tab['sample'][w]
    cal = tab['caltype'][w]
   
    subtype = tab['subtype'][w] 
    Ho_dist = tab['dist'][w]<0

    st1 = p1*(st - 1)
    st2 = p2*((st - 1))**2
    red = beta*(bv)
    
    mu_obs = mmax - p0 - st1 - st2 - red - alpha*(m_csp-np.median(m_csp)) 
    absmag = p0 + st1 + st2 + red + alpha*(m_csp-np.median(m_csp)) 

    mu_model = np.where(Ho_dist,distmod(h0,zhel,zcmb), dist)
    #mu_model = 5.0*np.log10(c*zcmb/h0) + 25.0
    
    yval = mmax - red - p0-st1-st2- alpha*(m_csp-np.median(m_csp))-mu_model
    pval = st1 + st2
    
    fac= (p1+(2*p2*st))

    err1 = ((fac*est)**2) +(emmax**2) +((beta*ebv)**2)+(2*beta*c_ms)+(2*fac*c_mbv)+(sig**2)+(((2.17*vel)/(zcmb*c))**2)+(alpha*(em/np.log(10)*m_csp))**2 

    err2 = ((fac*est)**2) +(emmax**2) +((beta*ebv)**2)+(2*fac*c_ms)+(beta*c_mbv)+(edist**2)



    err = np.where(Ho_dist,err1,err2)
    err = np.sqrt(err)
    dmu=mu_obs-mu_model

    wl = np.where(dmu >=0)
    wh = np.where(dmu <0)
    print (len(dmu[wl]),len(dmu[wh]))
    

    #pl.subplot(1,2,i+1)
   
    #pl.grid()
    w1 = np.where((dist<0) & (zcmb>0.01) & (st>.5) & (bv<.5))
    #pl.hist(dmu,histtype='step',lw=2,label=lab[i],bins=20)
    pl.hist(absmag[w1],histtype='step',lw=2,label=lab1[i],bins=20)

   
    #pl.ylim(-1.5,1.5),pl.xlim(0.0,0.15)
    pl.legend(numpoints =1,fontsize=12,loc='upper right')
    #pl.ylabel(r'$p0+p1*(s_{BV}-1)+p2*(s_{BV}-1)^2$' ,fontsize=12)
    #pl.label(r'$m-\beta(B-V)-\alpha M-p0-\mu$',fontsize=18)
    pl.xlabel(band+'_M (mag)' ,fontsize=14)
    #pl.axhline(sig[i],color='g')
    #pl.axhline(-sig[i],color='g')
    #pl.savefig('plots/hd_burns.pdf')


#pl.axvline(0,color='k',lw=3)
pl.tight_layout()
pl.savefig('../../plots/histabsmag_'+band+'_trgb.pdf')
#pl.show()






