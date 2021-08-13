from astropy.io import ascii
from astropy.table import join,vstack,unique
import matplotlib.pylab as pl


t1 = ascii.read('MarkList.txt')
t2= ascii.read('CSPHostMass.csv')

t3 = join(t1,t2, keys='sn')


t1.rename_column('sn','sn1')

t4 = join(t1,t2, keys='sn1')

final = vstack([t3,t4])

new = unique(final, keys='sn')

print new

new.write('91TMasses.dat', format='ascii.fixed_width', delimiter=' ',overwrite=True)



table = ascii.read('CSPIIHostMass.csv')

pl.hist(table['m'],color='.65')
pl.xlabel(r'$Log \ Host \ M_{\odot}$',fontsize=20), pl.ylabel('$Count$',fontsize=20)

pl.savefig('CSPIIHostMass.pdf')
pl.show()
