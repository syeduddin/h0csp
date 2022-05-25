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
path = ['g','r']

data = ascii.read(indir+'/AllRes.csv')
df3 = join(csp,data,keys='sn')

#pl.plot(df3['res_B'],df3['res_H'],'ko',mfc='white')
#pl.xlabel('B'), pl.ylabel('H')
#pl.savefig('../../plots/BH_res.pdf')
w = np.where(df3['proj']>50)
print (df3['res_H'][w])
#sys.exit()

#path = 'B'

pl.figure(figsize=(20,10))
#plot all residuals
for j in range(len(path)):
    
   
    # slopes and means
    sn = np.array(df3['sn'])
    z = np.array(df3['zcmb'+'_'+path[j]])
    proj = np.array(df3['proj'])
    res = np.array(df3['res'+'_'+path[j]])
    eres = np.array(df3['eres'+'_'+path[j]])
    #bv=np.array(df3['B-V'])
    #mass = np.array(df3['m'])
    #st = np.array(df3['st'])
    w= np.where(proj<100.)
    
    proj =proj[w]
    res=res[w]
    eres=eres[w]
    z =z[w]
    #bv = bv[w]
    #mass=mass[w]
    #st =st[w]
   
    
    med = np.median(proj)
    wl  = np.where(proj<10)
    wh  = np.where(proj>10)
    #print (len(proj[wl]), len(proj[wh]))
    #sys.exit()
    #wl  = np.where(proj<med)
    #wh  = np.where(proj>med)

    diff = (np.std(res[wl])-np.std(res[wh]))  

    #print '%6.3f'%diff
    #print np.median(z[wl]), np.median(z[wh])
    
    #print path[j],'&' '%6.3f'%(np.mean(res[wl])),'(','%6.3f'%(np.std(res[wl])),')','&' ,'%6.3f'%(np.mean(res[wh])),'(','%6.3f'%(np.std(res[wh])),')'
    
    print ('$',path[j],'$','&', '%.3f'%(np.std(res[wl])),'&' ,'%.3f'%(np.std(res[wh])), '&', '%.3f'%(np.std(res[wl])-np.std(res[wh])))

    pl.subplot(2,1,j+1)
    pl.grid()

    pl.errorbar(proj,res,yerr=eres,color='k',ls='none',ms=10)
    #pl.plot(proj,bv,'ro')

    cm = pl.cm.get_cmap('plasma')
    pl.scatter(proj,res,c=z,cmap=cm,s=100)
    pl.colorbar()
    #pl.colorbar().set_label(r'$Log \ host \ M_{\rm stellar} \ (M_{\odot})$', fontsize=14)
    #pl.colorbar().set_label(r'${(B-V) \ (mag)}$', fontsize=14)
    #pl.colorbar().set_label(r'$\Delta \mu \ ('+path[j][-1:]+') \ (mag)$', fontsize=14)
    #pl.colorbar().set_label(r'$sz$', fontsize=14)


    pl.xlabel(r'$Projected \ Distance \ (kpc)$',fontsize=14)
    pl.ylabel(r'$\Delta \mu \ ('+path[j][-1:]+') \ (mag)$' ,fontsize=14)
    #pl.ylabel(r'$(B-V) \ (mag)$' ,fontsize=14)

    pl.tick_params(axis='both', labelsize=14)
    #pl.ylim(-2,2),pl.xlim(0,60)
    #pl.axvline(10,c='k',ls='--',lw=2)
    
pl.tight_layout()
pl.savefig('../../plots/mu_sepsame.pdf',bbox_inches='tight', dpi=100)


