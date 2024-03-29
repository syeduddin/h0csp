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
path = ['u','B','g','V','r','i','Y','J','H']

#path = 'B'

pl.figure(figsize=(20,10))
#plot all residuals
for j in range(len(path)):
    data = ascii.read(indir+'/Ceph_res_'+path[j]+'_update2.csv')
   
  
    df3 = join(csp,data,keys='sn')
    
    # slopes and means
    sn = np.array(df3['sn'])
    z = np.array(df3['zcmb'])
    proj = np.array(df3['proj'])
    res = np.array(df3['res'])
    eres = np.array(df3['eres'])
    bv=np.array(df3['B-V'])
    mass = np.array(df3['m'])
    st = np.array(df3['st'])
    w= np.where((proj<100.) & (sn!='SN2011jn'))
    
    proj =proj[w]
    res=res[w]
    eres=eres[w]
    z =z[w]
    bv = bv[w]
    mass=mass[w]
    st =st[w]
    
    
    
    med = np.median(proj)
    wl  = np.where(proj<10)
    wh  = np.where(proj>10)
    #print (len(proj[wl]), len(proj[wh]))
    
    #wl  = np.where(proj<med)
    #wh  = np.where(proj>med)

    diff = (np.std(res[wl])-np.std(res[wh]))  

    #print '%6.3f'%diff
    #print np.median(z[wl]), np.median(z[wh])
    
    #print path[j],'&' '%6.3f'%(np.mean(res[wl])),'(','%6.3f'%(np.std(res[wl])),')','&' ,'%6.3f'%(np.mean(res[wh])),'(','%6.3f'%(np.std(res[wh])),')'
    print ('$',path[j],'$','&', '%.3f'%(np.std(res[wl])),'&' ,'%.3f'%(np.std(res[wh])), '&', '%.3f'%(np.std(res[wl])-np.std(res[wh])))

   

   
    pl.subplot(3,3,j+1)
   
    pl.grid()
    #pl.hist(bv[wl],bins=50,range=[-0.25,1.25],histtype='stepfilled',label='<10 kpc',color='r',alpha=.5)
    #pl.hist(bv[wh],bins=50,range=[-0.25,1.25],histtype='stepfilled',label='>10 kpc',color='b',alpha=.5)
    #pl.xlabel(r'$B-V \ (mag)$' ,fontsize=12),
    #pl.legend(loc='upper right')
    pl.errorbar(proj,res,yerr=eres,ls='None',color='k',alpha=.3)
    #pl.plot(proj,bv,'ro')

    cm = pl.cm.get_cmap('plasma')
    pl.scatter(proj,res,c=z,cmap=cm,s=100)
    pl.colorbar()
    #pl.colorbar().set_label(r'$Log \ host \ M_{\rm stellar} \ (M_{\odot})$', fontsize=14)
    #pl.colorbar().set_label(r'${(B-V) \ (mag)}$', fontsize=14)
    #pl.colorbar().set_label(r'$\Delta \mu \ ('+path[j][-1:]+') \ (mag)$', fontsize=14)
    #pl.colorbar().set_label(r'$s_{BV}$', fontsize=14)


    pl.xlabel(r'$Projected \ Distance \ (kpc)$',fontsize=14)
    pl.ylabel(r'$\Delta \mu \ ('+path[j][-1:]+') \ (mag)$' ,fontsize=14)
    #pl.ylabel(r'$Number$',fontsize=12)

    pl.tick_params(axis='both', labelsize=12)
    pl.ylim(-1.5,1.5),pl.xlim(0,60)
    #pl.axvline(10,c='k',ls='--',lw=2)
    
pl.tight_layout()
pl.savefig('../../plots/mu_sep_update2.pdf',bbox_inches='tight', dpi=100)
#pl.savefig('../../plots/bv_dist.pdf',bbox_inches='tight', dpi=100)

#pl.show()






