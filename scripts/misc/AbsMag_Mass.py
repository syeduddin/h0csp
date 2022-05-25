from astropy.io import ascii
import matplotlib.pylab as pl
import numpy as np
from scipy.stats import norm
from astropy.table import join, hstack, vstack,Table

t1 = ascii.read('../../data/working/B_all.csv')
t0 = ascii.read('../../data/working/calmorph.txt',format='csv')

t = join(t1,t0,keys='host')


w1 = np.where(t['caltype']=='c')
Mc = t['Mmax'][w1]-t['dist'][w1] +(t['st'][w1]-1)*1.17 -t['BV'][w1]*3.04
eMc = np.sqrt((t['eMmax'][w1])**2+(t['edist'][w1])**2)
mc = t['m'][w1]
emcl = t['m'][w1]-t['ml'][w1]
emcu = t['mu'][w1]-t['m'][w1]
print ('Ceph ',len(mc))
w2 = np.where(t['caltype']=='t')
Mt = t['Mmax'][w2]-t['dist'][w2] +(t['st'][w2]-1)*1.41 -t['BV'][w2]*3.07
eMt = np.sqrt((t['eMmax'][w2])**2+(t['edist'][w2])**2)
mt = t['m'][w2]
emtl = t['m'][w2]-t['ml'][w2]
emtu = t['mu'][w2]-t['m'][w2]
print ('TRGB ',len(mt))
w3 = np.where(t['caltype']=='s')
Ms = t['Mmax'][w3]-t['dist'][w3] +(t['st'][w3]-1)*1.07 -t['BV'][w3]*2.88
eMs = np.sqrt((t['eMmax'][w3])**2+(t['edist'][w3])**2)
ms = t['m'][w3]
emsl = t['m'][w3]-t['ml'][w3]
emsu = t['mu'][w3]-t['m'][w3]
print ('SBF: ',len(ms))


#t['host'][w3].more()
#print (t['host'][w3])




wl = np.where(t['morph'][w1]=='l')
we = np.where(t['morph'][w1]=='e')

pl.errorbar(mc,Mc,yerr=eMc,xerr=[emcl,emcu],fmt='s',color='b',label='Cepheid',ms=8)
#pl.errorbar(mc[wl],Mc[wl],yerr=eMc[wl],xerr=[emcl[wl],emcu[wl]],color='b',fmt='s',label='Late',ms=8)
#pl.errorbar(mc[we],Mc[we],yerr=eMc[we],xerr=[emcl[we],emcu[we]],color='r',fmt='s',label='Early',ms=8)

##w2
wl = np.where(t['morph'][w2]=='l')
we = np.where(t['morph'][w2]=='e')

pl.errorbar(mt,Mt,yerr=eMt,xerr=[emtl,emtu],fmt='d',color='g',label='TRGB',ms=10)
#pl.errorbar(mt[wl],Mt[wl],yerr=eMt[wl],xerr=[emtl[wl],emtu[wl]],color='b',fmt='d',ms=8)
#pl.errorbar(mt[we],Mt[we],yerr=eMt[we],xerr=[emtl[we],emtu[we]],color='r',fmt='d',ms=8)

## w3
wl = np.where(t['morph'][w3]=='l')
we = np.where(t['morph'][w3]=='e')

pl.errorbar(ms,Ms,yerr=eMs,xerr=[emsl,emsu],fmt='o',color='r',label='SBF',ms=8)
#pl.errorbar(ms[wl],Ms[wl],yerr=eMs[wl],xerr=[emsl[wl],emsu[wl]],color='b',fmt='o',ms=8)
#pl.errorbar(ms[we],Ms[we],yerr=eMs[we],xerr=[emsl[we],emsu[we]],color='r',fmt='o',ms=8)



pl.xlabel('$Host \ Mass$'), pl.ylabel('$M_B$')
pl.grid(),pl.legend(numpoints=1,loc='lower left')
pl.savefig('../../plots/absmag_mass.pdf')






