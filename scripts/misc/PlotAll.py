import matplotlib.pylab as pl
from astropy.io import ascii
import numpy as np
from astropy.table import join
import sys
tab = ascii.read('../../data/lc/B_max.csv')
tab1 = ascii.read('../../data/working/B_ceph.csv')
tab2 = ascii.read('../../data/working/B_trgb.csv')
tab3 = ascii.read('../../data/working/B_sbf.csv')
w1=np.where(tab1['dist']>1.0)
w2=np.where(tab2['dist']>1.0)
w3=np.where(tab3['dist']>1.0)
w=np.where(tab1['dist']==-1.0)

s = np.where(tab1['dist']==-1)

s1 = np.where((tab1['sample']=='CSPI') &(tab1['dist']==-1.0) )
s2 = np.where((tab1['sample']=='CSPII') & (tab1['dist']==-1.0))


pl.figure(1)
pl.errorbar(tab['st'],tab['BV'],xerr=tab['est'],yerr=tab['eBV'],mfc='white',fmt='k.', alpha=.5,ms=10,label='$CSP$')
pl.errorbar(tab1['st'][w1],tab1['BV'][w1],xerr=tab1['est'][w1],yerr=tab1['eBV'][w1],color='b', fmt='o', ms=8,label='$Cepheid$')
pl.errorbar(tab2['st'][w2],tab2['BV'][w2],xerr=tab2['est'][w2],yerr=tab2['eBV'][w2],color='r',fmt='*',ms=10,label='$TRGB$')
pl.errorbar(tab3['st'][w3],tab3['BV'][w3],xerr=tab3['est'][w3],yerr=tab3['eBV'][w3],color='#377eb8',fmt='d',ms=6,label='$SBF$')
pl.legend(numpoints=1,loc='upper left')
pl.grid()
pl.xlabel(r'$Color-stretch \ parameter \ (s_{BV})$',fontsize=14)
pl.ylabel(r'$Color \ (B-V)$',fontsize=14)
pl.savefig('../../plots/st_bv.pdf')


pl.figure(2)
pl.subplot(1,2,1)
pl.hist(tab1['m'][w],bins=50,range=[7,12],histtype='stepfilled',label='CSP',color='k', alpha=.3)
pl.hist(tab1['m'][s1],bins=50,range=[7,12],histtype='step',label='CSPI',color='r',ls='-',lw=2)
pl.hist(tab1['m'][s2],bins=50,range=[7,12],histtype='step',label='CSPII',color='b',ls='-',lw=2)
pl.xlabel(r'$Log \ host \ M_{stellar} \ (M_{\odot})$',fontsize=18), pl.ylabel('$Number$',fontsize=18)
pl.legend(loc='upper left'),pl.grid()

pl.subplot(1,2,2)
pl.hist(tab1['m'][w1],bins=10,range=[7,12],histtype='step',label='Cepheid',color='b',lw=2)
pl.hist(tab2['m'][w2],bins=10,range=[7,12],histtype='step',label='TRGB',color='r',lw=2)
pl.hist(tab3['m'][w3],bins=10,range=[7,12],histtype='step',label='SBF',color='#377eb8',lw=2)

pl.xlabel(r'$Log \ host \ M_{stellar} \ (M_{\odot})$',fontsize=18), pl.ylabel('$Number$',fontsize=18)
#pl.ylim(0,30),pl.xlim(7,12)
pl.legend(loc='upper left'),pl.grid()
pl.savefig('../../plots/massdist.pdf')


pl.figure(3)



pl.hist(tab1['zcmb'][s],range=[0,.14],bins=20,histtype='stepfilled',label='CSP',color='k', alpha=.3)
pl.hist(tab1['zcmb'][s1],range=[0,.14],bins=20,histtype='step',label='CSPI',color='r',ls='-',lw=2)
pl.hist(tab1['zcmb'][s2],range=[0,0.14],bins=20,histtype='step',label='CSPII',color='b',ls='-',lw=2)
pl.legend(loc='upper right'),pl.grid()
pl.xlabel(r'$z_{CMB}$',fontsize=14), pl.ylabel('$Number$',fontsize=14)

pl.savefig('../../plots/zdist.pdf')