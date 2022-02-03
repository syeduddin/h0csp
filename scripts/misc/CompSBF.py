from astropy.io import ascii
from astropy.table import join,vstack, Table, unique
import numpy as np
import sys
import matplotlib.pylab as pl
from scipy.stats import linregress
from collections import Counter
filter = 'B'

t1 = ascii.read('../../data/lc/'+filter+'_max_SBF_old.dat')
t2 = ascii.read('../../data/lc/'+filter+'_max_SBF.dat')


         




t =join(t1,t2,keys='name')

pl.errorbar(t['BV_1'],t['BV_2'],xerr=t['eBV_1'],yerr=t['eBV_2'],color='r',fmt='o')
#pl.errorbar(t['Mmax_1'],t['Mmax_2'],xerr=t['eMmax_1'],yerr=t['eMmax_2'],color='r',fmt='o')

pl.xlabel('old'),pl.ylabel('new')
pl.savefig('../../plots/compSBFBmax.pdf')

print (np.mean(t['BV_1']-t['BV_2']))
m, c, r_value, p_value, std_err = linregress(t['Mmax_1'],t['Mmax_2'])
print ('Mmax:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)

m, c, r_value, p_value, std_err = linregress(t['BV_1'],t['BV_2'])
print ('B-V:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)



m, c, r_value, p_value, std_err = linregress(t['st_1'],t['st_2'])
print ('st:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)

