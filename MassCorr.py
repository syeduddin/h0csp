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


f1  = pd.read_csv('data/CSPHostMass.csv',delimiter=',')
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
pl.figure(figsize=(20,10))
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
    lm = linmix.LinMix(mass, res, emass, eres, K=2)
    lm.run_mcmc(silent=True)
    sl.append('%6.3f'%(np.mean(lm.chain['beta'])))
    esl.append('%6.3f'%np.std(lm.chain['beta']))


    print '& $'+path[j][-1:]+'$', '&%.3f'%np.mean(lm.chain['beta']),'(','%.3f'%(np.std(lm.chain['beta'])),')&', '%.3f'%(mean_x1_high-mean_x1_low),'(', '%.3f'%(np.sqrt((error_x1_low**2)+(error_x1_high**2))),')'+'\ \\'

    
    xl = np.array([5, 15])    
    
    pl.subplot(3,3,j+1)
   
    pl.grid()
    pl.axvline(np.median(m),c='k',ls='-',lw=3)
    
    pl.errorbar(mass,res,xerr=[emassl,emassu],yerr=eres,fmt='o',color='gray',mfc='white')
    pl.xlim(7,12),pl.ylim(-1.0,1.0)
    pl.xlabel(r'$Log \ host \ M_{\rm stellar} \ (M_{\odot})$',fontsize=18)
    pl.ylabel(r'$\Delta \mu ('+path[j][-1:]+') \ (mag)$' ,fontsize=18)
    pl.tick_params(axis='both' )
    for i in range(0, len(lm.chain), 25):
        xs =  np.array([5, 15])
        ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
        pl.plot(xs, ys, color='b', alpha=0.02)
    pl.plot(xs, np.mean(lm.chain['alpha']) + xs *np.mean(lm.chain['beta']), color='b',lw=3)
    pl.plot([7,np.median(mass)],[mean_x1_low,mean_x1_low],ls='-',lw=3,c='r')
    pl.plot([np.median(mass),12],[mean_x1_high,mean_x1_high],ls='-',lw=3,c='r')
    pl.legend()
#pl.tight_layout()

#pl.savefig('plots/AllMassCorr.pdf',bbox_inches='tight', dpi=100)

   
#pl.show()        

print 'slopes'
print (', '.join(sl))
print (', '.join(esl))

print 'Offsets'
print (', '.join(off))
print (', '.join(eoff))
