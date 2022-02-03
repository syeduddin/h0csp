from astropy.io import ascii
import matplotlib.pylab as pl
import numpy as np
from scipy.stats import norm
from astropy.table import join
import os,sys

from scipy.stats import linregress

c0 = 300000./1.029e+14


tab = ascii.read('../../data/working/B_ceph.csv')

w = np.where(tab['caltype']=='c')  
dist = 10**((tab['dist'][w]/5)+1)/1e6

print (tab['dist'][w])
edist =0.461*dist*tab['edist'][w]
z  = tab['zcmb'][w]

#edist = edist
m, c, r_value, p_value, std_err = linregress(dist,z*c0)

print ((np.mean(z*c0/dist))*1.029e+14)


pl.errorbar(dist,z*c0,xerr=edist,fmt='o',color='r')
pl.grid
pl.savefig('../../directH0SBF.pdf')
#os.system('open ../../directH0SBF.pdf')

