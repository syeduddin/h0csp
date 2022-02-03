import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import numpy as np
import glob
import pandas as pd
import sys
from numpy import polyfit
import random
from scipy.stats import ks_2samp

from astropy.io import ascii
from astropy.table import join

csp=ascii.read('../../data/hosts/proj_dist.csv')



indir = '../../results/'
#path = ['u','B','g','V','r','i','Y','J','H']

path = 'B'

#pl.figure(figsize=(20,10))
#plot all residuals
for j in range(len(path)):
    data = ascii.read(indir+'/Ceph_res_'+path[j]+'.csv')
   
  
    df3 = join(csp,data,keys='sn')
    
    # slopes and means
    z = np.array(df3['zcmb'])
    proj = np.array(df3['proj'])
    res = np.array(df3['res'])
    eres = np.array(df3['eres'])
    bv=np.array(df3['B-V'])
    mass = np.array(df3['m'])
    st = np.array(df3['st'])
    w= np.where(proj<100.)
    
    proj =proj[w]
    res=res[w]
    eres=eres[w]
    z =z[w]
    bv = bv[w]
    mass=mass[w]
    st =st[w]
    #print np.max(ebv)
    
    wl  = np.where(proj<10)
    wh  = np.where(proj>10)
    diff = (np.std(res[wl])-np.std(res[wh]))  

    #print '%6.3f'%diff
    #print np.median(z[wl]), np.median(z[wh])
    
    #print path[j],'&' '%6.3f'%(np.mean(res[wl])),'(','%6.3f'%(np.std(res[wl])),')','&' ,'%6.3f'%(np.mean(res[wh])),'(','%6.3f'%(np.std(res[wh])),')'
    #print '$',path[j],'$','&', '%.3f'%(np.std(res[wl])),'&' ,'%.3f'%(np.std(res[wh])), '&', '%.3f'%(np.std(res[wl])-np.std(res[wh]))

    ks= ks_2samp(res[wl],res[wh])
    m = float(len(res[wl]))
    n = float(len(res[wh]))
    mn= np.sqrt((n+m)/(n*m))

    c_alpha=[1.073,1.138,1.224,1.358,1.48,1.628,1.731,1.949]
    
    print ('$',path[j],'$','&' '%.3f'%ks[0],'&''%.3f'%ks[1],'&''%.3f'%(1.358*mn))
    print (len(z[wh]))

# Monte Carlo
    #diff=[]

    #for i in range(0,10):
     #   nl=random.sample(res[wl], len(res[wh])) # Raondom low-z
      #  nh = res[wh]
       # n1 = np.random.normal(0,np.std(res),len(res))    
        #n2 = np.random.normal(0,np.std(res),len(res[wh]))
    
        #ks= ks_2samp(n1,n2)[1]
        #if ks>0.05:
         #   diff.append(np.std(nl)-np.std(nh))
    #print ('Mone Carlo:')
    #print (len(diff)/10.)
    #print (np.mean(diff), np.std(diff))
    
    #print ('%.3f'%(np.std(res)), '%.3f'%np.sqrt(np.mean((res**2))))
    #pl.subplot(3,3,j+1)
   
    pl.grid()

    pl.errorbar(proj,res,yerr=eres,ls='None',color='k',alpha=.3)
    #pl.plot(proj,bv,'ro')

    cm = pl.cm.get_cmap('RdYlBu')
    pl.scatter(proj,res,c=st,cmap=cm,s=100)
    #pl.colorbar().set_label(r'$Log \ host \ M_{\rm stellar} \ (M_{\odot})$', fontsize=14)
    pl.colorbar().set_label(r'$s_{BV}$', fontsize=14)


    pl.xlabel(r'$Projected \ Distance \ (kpc)$',fontsize=14)
    pl.ylabel(r'$\Delta \mu \ ('+path[j][-1:]+') \ (mag)$' ,fontsize=14)
    pl.tick_params(axis='both', labelsize=14)
    pl.ylim(-1.5,1.5),pl.xlim(0,60)
    #pl.axvline(10,c='k',ls='--',lw=2)
    
pl.tight_layout()
pl.savefig('../../plots/mu_sepSt.pdf',bbox_inches='tight', dpi=100)

#pl.show()






