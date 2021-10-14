import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import numpy as np
import glob
import pandas as pd
import sys
from scipy.stats import pearsonr
from numpy import polyfit
import random
from scipy.odr import *
from astropy.io import ascii
from astropy.table import join

csp=ascii.read('../data/proj_snhost.csv')


#name=[]
#for i in range(len(sn)):
 #   name.append('SN'+'20'+sn[i][2:])

#name=np.array(name)
#proj =np.array(proj)




#df1 = pd.DataFrame({'name':name.byteswap().newbyteorder(),'proj':proj.byteswap().newbyteorder()})


indir = '/Users/suddin/Dropbox/CSP/residuals_CSPI+CSPII/'
#path = glob.glob(indir+'*') # each SN
path = ['u','B','g','V','r','i','Y','J','H']

pl.figure(figsize=(20,10))
#plot all residuals
for j in range(len(path)):
    data = ascii.read(indir+path[j]+'/resids_cv.dat')
   
  
    df3 = join(csp,data,keys='name')
    
    # slopes and means
    z = np.array(df3['zcmb'])
    proj = np.array(df3['proj'])
    res = np.array(df3['res'])
    eres = np.array(df3['eDMcv'])
    ebv=np.array(df3['EBV'])

    std =3* np.std(res)
    w= np.where(proj<100.)
    
    proj =proj[w]
    res=res[w]
    eres=eres[w]
    z =z[w]
    #print np.max(ebv)
    
    wl  = np.where(proj<10)
    wh  = np.where(proj>10)
#    diff = (np.std(res)-np.std(res[wh]))/np.std(res)
    diff = (np.std(res[wl])-np.std(res[wh]))  
    #print '%6.3f'%diff
    #print np.median(z[wl]), np.median(z[wh])
    
    #print path[j],'&' '%6.3f'%(np.mean(res[wl])),'(','%6.3f'%(np.std(res[wl])),')','&' ,'%6.3f'%(np.mean(res[wh])),'(','%6.3f'%(np.std(res[wh])),')'
    print ('$',path[j],'$','&', '%.2f'%(np.std(res[wl])),'&' ,'%.2f'%(np.std(res[wh])))

    #print ('%.3f'%(np.std(res)), '%.3f'%np.sqrt(np.mean((res**2))))
    pl.subplot(3,3,j+1)
   
    pl.grid()

    pl.errorbar(proj,res,yerr=eres,fmt='o',color='k',ms=10,mfc='white')
    
    #pl.errorbar(df3['proj'],df3['res'],yerr=[df3['eDMcv'],df3['eDMcv']],fmt='o',color='m',ms=10)

    pl.xlabel(r'$Projected \ Distance \ (kpc)$',fontsize=18)
    pl.ylabel(r'$\Delta \mu ('+path[j][-1:]+') \ (mag)$' ,fontsize=18)
    pl.tick_params(axis='both', labelsize=14)
    #pl.axvline(10,c='k',ls='--',lw=2)
    
pl.tight_layout()
pl.savefig('../plots/mu_sep.pdf',bbox_inches='tight', dpi=100)

#pl.show()






