from astropy.io import ascii
import numpy as np
from astropy.table import join,vstack,unique
import matplotlib.pylab as pl

csp=ascii.read('data/CSPI+II_max/B_max.dat')

sbf = ascii.read('data/khetan/SBF.dat')



sbf.rename_column('SN_name','name')

t2 = join(csp,sbf, keys='name')
w=np.where(t2['mu_sbf']>0.)

print np.mean(t2['st']-t2['stretch'])
print np.mean(t2['BV']-t2['color'])
print t2
sys.exit()
pl.figure(figsize=(20,10))

pl.subplot(1,3,1)
pl.errorbar(t2['st'],t2['stretch'],xerr=t2['est'],yerr=t2['err_stretch'],fmt='o',color='.65',ms=10)
pl.errorbar(t2['st'][w],t2['stretch'][w],xerr=t2['est'][w],yerr=t2['err_stretch'][w],fmt='ro',ms=10)
x=np.arange(0,1.5,.1)
pl.plot(x,x)
pl.xlim(.5,1.2),pl.ylim(.5,1.4)
pl.xlabel('$s_{BV} \ (This \ Work)$'),pl.ylabel('$s_{BV} \ (Khetan \ 2021)$' )
pl.grid()


pl.subplot(1,3,2)
pl.errorbar(t2['Mmax'],t2['Bmax'],xerr=t2['eMmax'],yerr=t2['eBmax'],fmt='o',color='.65',ms=10)
pl.errorbar(t2['Mmax'][w],t2['Bmax'][w],xerr=t2['eMmax'][w],yerr=t2['eBmax'][w],fmt='ro',ms=10)
x=np.arange(10,20,1)
pl.plot(x,x)
pl.xlim(11,19),pl.ylim(11,19)
pl.xlabel('$B_{max} \ (This \ Work)$'),pl.ylabel('$B_{max} \ (Khetan \ 2021)$' )
pl.grid()

pl.subplot(1,3,3)
pl.errorbar(t2['BV'],t2['color'],xerr=t2['eBV'],yerr=t2['e_color'],fmt='o',color='.65',ms=10)
pl.errorbar(t2['BV'][w],t2['color'][w],xerr=t2['eBV'][w],yerr=t2['e_color'][w],fmt='ro',ms=10)
x=np.arange(-.2,5,.1)
pl.plot(x,x)
pl.xlim(-0.1,.4),pl.ylim(-.1,.4)
pl.xlabel('$B-V \ (This \ Work)$'),pl.ylabel('$B-V \ (Khetan \ 2021)$' )
pl.grid()
pl.savefig('data/comp_lc.pdf')
pl.show()
