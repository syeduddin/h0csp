import astropy.io.fits as pyfits
#import matplotlib.pylab as pl
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

#path = 'u'

#pl.figure(figsize=(20,10))
#plot all residuals
for j in range(len(path)):
    data = ascii.read(indir+'/Ceph_res_'+path[j]+'.csv')
   
  
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
    w= np.where(proj<100.)
    
    proj =proj[w]
    res=res[w]
    eres=eres[w]
    z =z[w]
    bv = bv[w]
    mass=mass[w]
    st =st[w]
    #ww = np.where((res<-.8) & (proj>15)) # The outlier
    #print (sn[ww])
    
    
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
    #print '$',path[j],'$','&', '%.3f'%(np.std(res[wl])),'&' ,'%.3f'%(np.std(res[wh])), '&', '%.3f'%(np.std(res[wl])-np.std(res[wh]))

    ks= ks_2samp(res[wl],res[wh])
    m = float(len(res[wl]))
    n = float(len(res[wh]))
    mn= np.sqrt((n+m)/(n*m))
   
    d_stat = ks[0] 
    
    def c_alpha(alpha):
        #a1 = 1+(m/n); a2 = 2*m; a3 = -np.log(alpha/2.)
        c_alpha= np.sqrt((-np.log(alpha/2.)*.5))
        return(c_alpha*mn)
    
    #print (d_stat)
    #print (path[j])
    for a in (np.arange(0.01,1.0,0.01)):
        
        c_a = c_alpha(a)
        #print(a,c_a)
        if d_stat > c_a:
            print ('$',path[j],'$','&','&','%.3f'%d_stat,'&','%.3f'%c_a,'%d'%((1-a)*100),'$\%$')
            break 

        #c_alpha=[1.073,1.138,1.224,1.358,1.48,1.628,1.731,1.949]
    
    #print ('$',path[j],'$','&' '%.3f'%ks[0],'&''%.3f'%ks[1],'&''%.3f'%(1.358*mn))
    #print (len(z[wh]))


   
    








