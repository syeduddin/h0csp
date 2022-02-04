from astropy.io import ascii
from astropy.table import join
import matplotlib.pylab as pl

t1 = ascii.read('../../results/Ceph_res_B.csv')

t2 = ascii.read('/Users/suddin/Dropbox/CSP/residuals_CSPI+CSPII/B/resids_cv.dat')

t3 = ascii.read('../../data/lc/Spreadsheet.csv')

t3.rename_column('Name','name')

t = join(t2,t3, keys='name')

pl.plot(t['res_1'],t['res_2'],'ro')
pl.show()
