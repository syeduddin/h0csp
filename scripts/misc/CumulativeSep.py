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

#path = 'J'

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
    #print (np.min(proj))
    
    
    med = np.median(proj)
    wl  = np.where(proj<10)
    wh  = np.where(proj>10)
    #print (len(proj[wl]), len(proj[wh]))
    pl.subplot(3,3,j+1)
    pl.grid()
#https://www.youtube.com/watch?v=ap4mfGvgDsM
    xl = np.sort(res[wl])
    yl = np.arange(1, len(xl)+1)/len(xl)
    _= pl.plot(xl,yl,ls='-',lw=2,label='<10 kpc')
    pl.margins(0.02)

    xh = np.sort(res[wh])
    yh = np.arange(1, len(xh)+1)/len(xh)
    _= pl.plot(xh,yh,ls='-',lw=2,label='>10 kpc')
    pl.margins(0.02)
    pl.legend()
    pl.xlabel(r'$\Delta \mu \ ('+path[j][-1:]+') \ (mag)$' ,fontsize=14)
    pl.ylabel(r'$CDF$',fontsize=14)
    #print (np.max(xl-xh))
    #values, base = np.histogram(res[wl], bins=100)
    #cumulative = np.cumsum(values)
    #pl.plot(base[:-1], cumulative, c='blue',)
    #values, base = np.histogram(res[wh], bins=100)
    #cumulative = np.cumsum(values)
    #pl.plot(base[:-1], cumulative, c='red')

    
    #pl.xlabel(r'$Projected \ Distance \ (kpc)$',fontsize=14)
    #pl.ylabel(r'$\Delta \mu \ ('+path[j][-1:]+') \ (mag)$' ,fontsize=14)
    #pl.ylabel(r'$(B-V) \ (mag)$' ,fontsize=14)

    pl.tick_params(axis='both', labelsize=14)
    #pl.ylim(-1.5,1.5),pl.xlim(0,60)
    #pl.axvline(10,c='k',ls='--',lw=2)
    
pl.tight_layout()
pl.savefig('../../plots/cdfResAll_update2.pdf',bbox_inches='tight', dpi=100)

#pl.show()






