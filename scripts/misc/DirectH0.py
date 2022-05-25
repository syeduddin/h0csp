from astropy.io import ascii
import matplotlib.pylab as pl
import numpy as np
from scipy.stats import norm
import os,sys
from astropy.table import hstack,join
from scipy.stats import linregress




t1 = ascii.read('../../data/calibrators/need_sbfj21.csv',format='csv')
t2 = ascii.read('../../data/working/B_sbf.csv')

t1.rename_column('gal','name')


#t2.rename_column('Host','name')
#t2['name'].fill_value ='a b'

#print(t2['name'])
#host = join(t1,t2.filled(),keys='name')
#print (host)

#match= ((set(t1['sn']).difference(t2['sn'])))

match = join(t1,t2,keys='sn')

print (match)

tab = ascii.read('../../data/calibrators/sbf_vel_dist.csv')
name =tab['name']
dist = tab['dist']
edist = tab['edist']
vel = tab['vel']
evel = tab['evel']

match = join(t1,tab,keys='name')

#print (match['sn'])

match= ((set(t1['name']).difference(name)))

#print (match)
sys.exit()

m, c, r_value, p_value, std_err = linregress(dist,vel)

print (m,std_err)


pl.errorbar(dist,vel,xerr=edist,yerr=evel,fmt='o',color='r')
pl.grid
pl.savefig('../../directH0SBF.pdf')
os.system('open ../../directH0SBF.pdf')

