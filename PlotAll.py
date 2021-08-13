import matplotlib.pylab as pl
from astropy.io import ascii
import numpy as np
from astropy.table import join

tab = ascii.read('data/CSPI+II_max/B_max.dat')
tab1 = ascii.read('data/B_ceph.dat')
tab2 = ascii.read('data/B_trgb.dat')
tab3 = ascii.read('data/B_sbf.dat')
w1=np.where(tab1['dist']>1.0)
w2=np.where(tab2['dist']>1.0)
w3=np.where(tab3['dist']>1.0)
w=np.where(tab1['dist']==-1.0)



pl.figure(1)
pl.errorbar(tab['st'],tab['BV'],xerr=tab['est'],yerr=tab['eBV'],fmt='ko',alpha=.5,ms=10,label='$CSP$')
pl.errorbar(tab1['st'][w1],tab1['BV'][w1],xerr=tab1['est'][w1],yerr=tab1['eBV'][w1],fmt='ro',ms=10,label='$Cepheid$')
pl.errorbar(tab2['st'][w2],tab2['BV'][w2],xerr=tab2['est'][w2],yerr=tab2['eBV'][w2],fmt='go',ms=10,label='$TRGB$')
pl.errorbar(tab3['st'][w3],tab3['BV'][w3],xerr=tab3['est'][w3],yerr=tab3['eBV'][w3],fmt='bo',ms=10,label='$SBF$')
pl.legend(numpoints=1,loc='upper left')
pl.grid()
pl.xlabel(r'$Color-stretch \ parameter \ (s_{BV})$',fontsize=20)
pl.ylabel(r'$Color \ (B-V)$',fontsize=20)
pl.savefig('plots/st_bv.pdf')
#pl.show()

pl.figure(2)

pl.hist(tab1['m'][w],bins=50,histtype='stepfilled',alpha=.5,label='SN Ia',color='k')
pl.hist(tab1['m'][w1],bins=20,histtype='stepfilled',alpha=.5,label='Cepheid',color='r')
pl.hist(tab2['m'][w2],bins=20,histtype='stepfilled',alpha=.5,label='TRGB',color='g')
pl.hist(tab3['m'][w3],bins=20,histtype='stepfilled',alpha=.5,label='SBF',color='b')

pl.xlabel(r'$Log \ host \ M_{stellar} \ (M_{\odot})$',fontsize=20), pl.ylabel('$Number$',fontsize=20)
pl.ylim(0,30),pl.xlim(7,12)
pl.legend(loc='upper left'),pl.grid()
pl.savefig('plots/massdist.pdf')


t1 = join(tab3,tab2,keys='sn')
p1 = np.where((t1['dist_1']>1.0) & (t1['dist_2']>1.0))

print t1['sn'][p1]#, t1['dist_1'][p1], t1['dist_2'][p1]

pl.close('all')
pl.figure(3)

pl.errorbar(t1['dist_1'][p1],t1['dist_2'][p1],xerr=t1['edist_1'][p1],yerr=t1['edist_2'][p1],fmt='ro',ms=10,label='$SBF-TRGB$')
pl.grid()
x = np.arange(30.5,32,.1)
pl.plot(x,x,'k')
pl.xlim(30.8,31.8),pl.ylim(30.8,31.8)
pl.xlabel(r'$\mu_{SBF} \ (mag)$',fontsize=20)
pl.ylabel(r'$\mu_{TRGB} \ (mag)$',fontsize=20)
pl.savefig('plots/caldiff.pdf')

pl.show()


