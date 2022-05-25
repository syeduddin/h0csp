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
import linmix
from astropy.io import ascii
from astropy.table import Table

# Getting results

sl=[]
esl=[]
off=[]
eoff=[]
pl.figure(figsize=(20,10))

filter = ['u','B','g','V','r','i','Y','J','H']

#filter='B'
for j in range(len(filter)):
    
    
    
    
    tab = ascii.read('../../results/Ceph_res_nohm'+filter[j]+'.csv')
   
    w =np.where((tab['sample']=='CSPI')& (tab['cal']=='none'))
    #w = np.where(tab['cal']=='none')
    mass = tab['m'][w]
    print len(tab['sn'])
    
    
    ml =  tab['m'][w]-tab['ml'][w]
    mu = tab['mu'][w]- tab['m'][w]
    em = (ml+mu)/2.
    res=tab['res'][w]
    eres=tab['eres'][w]

    for n,i in enumerate(mass):
        if i==11.5: mass[n]=random.uniform(7.1,7.9)
            
    
    wl=np.where(mass<np.median(mass))
    wh=np.where(mass>np.median(mass))
    
    err_int=(np.std(res)) - (np.mean(eres))
    
    wt= np.sqrt(1/((eres**2)+(err_int**2))) # weights
    #wt= 1/((eres**2.)+(err_int**.2)) # weights
    
    # Low
    mean_x1_low= np.sum(res[wl]*wt[wl])/np.sum(wt[wl])
    error_x1_low= np.sqrt((1/np.sum(wt[wl])))
    
    #high
    
    mean_x1_high= np.sum(res[wh]*wt[wh])/np.sum(wt[wh])
    error_x1_high =np.sqrt((1/np.sum(wt[wh])))
    
    off.append('%6.2f'%(mean_x1_high-mean_x1_low))
    eoff.append('%6.2f'%(np.sqrt((error_x1_low**2)+(error_x1_high**2))))
    

    #LINMIX
    lm = linmix.LinMix(mass, res, em, eres, K=2)
    lm.run_mcmc(silent=True)
    sl.append('%6.3f'%(np.mean(lm.chain['beta'])))
    esl.append('%6.3f'%np.std(lm.chain['beta']))


    #print '& $',filter[j],'$', '&%.3f'%np.mean(lm.chain['beta']),'(%.3f)'%(np.std(lm.chain['beta'])),'&', '%.3f'%(mean_x1_high-mean_x1_low), '(%.3f)'%(np.sqrt((error_x1_low**2)+(error_x1_high**2))),'&','%.2f'%np.median(mass),'\\\\' 
    #sys.exit()
    xl = np.array([5, 15])    
    
    pl.subplot(3,3,j+1)
   
    pl.grid()
    pl.axvline(np.median(mass),c='k',ls='-',lw=3)
    
    pl.errorbar(mass,res,xerr=[ml,mu],yerr=eres,fmt='o',color='gray',mfc='white')
    pl.xlim(7,12),pl.ylim(-1.5,1.5)
    pl.xlabel(r'$Log \ host \ M_{\rm stellar} \ (M_{\odot})$',fontsize=18)
    pl.ylabel(r'$\Delta \mu ('+filter[j]+') \ (mag)$' ,fontsize=18)
    pl.tick_params(axis='both' )
    for i in range(0, len(lm.chain), 25):
        xs =  np.array([5, 15])
        ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
        pl.plot(xs, ys, color='b', alpha=0.02)
    pl.plot(xs, np.mean(lm.chain['alpha']) + xs *np.mean(lm.chain['beta']), color='b',lw=3)
    pl.plot([7,np.median(mass)],[mean_x1_low,mean_x1_low],ls='-',lw=3,c='r')
    pl.plot([np.median(mass),12],[mean_x1_high,mean_x1_high],ls='-',lw=3,c='r')
    #pl.legend()
pl.tight_layout()

pl.savefig('../../plots/AllMassCorrnoHM_csp1.pdf',bbox_inches='tight', dpi=100)
    
   

print ('slopes')
print (', '.join(sl))
print (', '.join(esl))

print ('Offsets')
print (', '.join(off))
print (', '.join(eoff))



