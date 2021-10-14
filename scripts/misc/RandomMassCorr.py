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
import linmix


tbg=['SN2005M','SN2005eq','SN2008fw','SN2007S','SN2007ai','SN2005bl','SN2005ke','SN2006bd','SN2006gt','SN2006mr','SN2007N','SN2007ax','SN2007ba','SN2008bd','SN2008bi','SN2008bt','SN2009F']


f1  = pd.read_csv('CSPHostMass.csv',delimiter=',')
d1=pd.DataFrame(f1)
m = d1['m']
eml = d1['m']-d1['ml']
emu = d1['mu']-d1['m']
snname=d1['sn']
snname1=d1['sn1']

#wp = np.where(np.in1d(snname, tbg))[0]
#snname= snname[~np.in1d(range(len(snname)),wp)]
#snname1= snname1[~np.in1d(range(len(snname1)),wp)]
#m=m[~np.in1d(range(len(m)),wp)]
#eml=eml[~np.in1d(range(len(eml)),wp)]
#emu=emu[~np.in1d(range(len(emu)),wp)]
em = (eml+emu)/2.
#print len(m)
#sys.exit()
    
    
indir = '/Users/suddin/Dropbox/CSP/residuals_CSPI+CSPII/'
#path = glob.glob(indir+'*') # each SN
path = ['u','B','g','V','r','i','Y','J','H']
#pl.figure(figsize=(20,10))
sl=[]
esl=[]
off=[]
eoff=[]
for j in range(len(path)):
    res=[]
    eres=[]
    mass=[]
    emass=[]
    emassl=[]
    emassu=[]
   
    data = pd.read_csv(indir+path[j]+'/resids_cv.dat',delim_whitespace=True)
    d2=pd.DataFrame(data)
    sn =d2['name'] 
    for i in range(len(sn)):
        
        if sn[i] in snname.values:w = np.where(sn[i]==snname)
        
        if sn[i] in snname1.values:w = np.where(sn[i]==snname1)
        res.append(d2['res'].values[i])
        eres.append(d2['eDMcv'].values[i])
        mass.append(*d1['m'].values[w])
        emass.append(*em.values[w])
        emassl.append(*eml.values[w])
        emassu.append(*emu.values[w])
        
    #Weigheted means
   
    #ws = np.where(np.isfinite(res))
    
    mass =np.array(mass)
    emassl =np.array(emassl)
    emassu =np.array(emassu)
    emass =np.array(emass)
    eres =np.array(eres)
    res=np.array(res)

    
    #ws = np.where((res>-1.) & (res<1.))
    ws = np.where(mass>0.)
    resmean= np.mean(res)
    cut = 3*np.std(res)
    cutl = resmean-cut
    cutu = resmean+cut

    #ws = np.where((res>cutl) & (res<cutu))
    mass=mass[ws]
    emass=emass[ws]
    emassl=emassl[ws]
    emassu=emassu[ws]
    res=res[ws]
    eres=eres[ws]
    
    wl=np.where(mass<10.5)
    wh=np.where(mass>10.5)
    
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
    

    
    
    

print 'Offsets'
print (', '.join(off))
print (', '.join(eoff))
