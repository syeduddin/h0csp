from astropy.io import ascii
from astropy.table import join,vstack, Table, unique, hstack
import numpy as np
import sys
import matplotlib.pylab as pl
from scipy.stats import linregress
from collections import Counter

filter = sys.argv[1]

t0 =  ascii.read('../../data/calibrators/sbf_jensen21.csv') # only SBF dist
t1 = ascii.read('../../data/lc/'+filter+'_max_SBF.dat') # khetan LC fit
t2 = ascii.read('../../data/lc/'+filter+'_max_SBFJ21.csv') # Jensen LC fit
t3 = ascii.read('../../data/working/'+filter+'_ceph.csv') # distant sample 

t0['ml'] = t0['m']- t0['em']
t0['mu']= t0['m']+ t0['em']


t3.remove_column('cosmo')
t3.remove_column('phys')

t1.rename_column('name','sn')

t2.rename_column('name','sn')
t2=unique(t2,keys='sn')

tt =vstack([t1,t2]) # LC

t0.rename_column('gal','host')
#t0['m']=0.0
#t0['ml'] = 0.0
#t0['mu'] = 0.0
t0['sample'] ='N/A'
t0['subtype'] ='N/A'
t0['quality'] =0
#t0['cosmo'] =9
#t0['phys'] ='N/A'
t0['caltype'] ='s'

t =join(tt,t0,keys='sn') # 17 from jensen 21

t.remove_column('Ho')
t.remove_column('eHo')
#t.remove_column('method')

t=unique(t,keys='sn')

#print (hstack([t['sn'],t['m']]))

for n,i in enumerate(t3['caltype']):
        if i=='c': t3['dist']=-1

w = np.where(t3['caltype']!='c')
          
t4 = vstack([t,t3[w]]) # 17 SBF

t4 = unique(t4,keys='sn',keep='first')
print (t)

#t.remove_column('method')


t.write('../../data/calibrators/calibrators_'+filter+'_sbfj21.csv', format='ascii.csv', delimiter=',',overwrite=True)
t4.write('../../data/working/'+filter+'_sbfj21.csv', format='ascii.csv', delimiter=',',overwrite=True)







 
