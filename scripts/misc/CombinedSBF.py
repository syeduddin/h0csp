from astropy.io import ascii
from astropy.table import join,vstack, Table, unique,hstack
import numpy as np
import sys
from collections import Counter


filter = sys.argv[1]

k =  ascii.read('../../data/calibrators/calibrators_sbf_final.csv') # only dist
j =  ascii.read('../../data/calibrators/sbf_jensen21.csv') # only dist
t1 = ascii.read('../../data/lc/'+filter+'_max_SBF.dat') # to get 7 
t2 = ascii.read('../../data/lc/'+filter+'_max_SBFJ21.csv') # for SBF  lc data
t3 = ascii.read('../../data/working/'+filter+'_ceph.csv') # distant sample 




j.rename_column('gal','host')

del k['zhel','zcmb','st', 'est', 'Mmax', 'eMmax', 't0', 'EBVmw', 'eEBVmw', 'covMs','BV','eBV','covBV_M','caltype','sample', 'subtype', 'quality', 'cosmo', 'phys']

print (k.colnames)
print (j.colnames)

j['ml'] = j['m']- j['em']
j['mu']= j['m']+ j['em']

j.remove_column('em')

t0 = vstack([j,k]) # dist+mass

t0=unique(t0,keys='sn',keep='first') # distances

kem = ((k['m']-k['ml'])+ (k['mu']-k['m']))/2.

#print (hstack([k['sn'],k['m'],kem]))

t3.remove_column('cosmo')
t3.remove_column('phys')

t1.rename_column('name','sn')
t2.rename_column('name','sn')

t2=unique(t2,keys='sn')

tt =vstack([t1,t2]) # LC 

t0['sample'] ='N/A'
t0['subtype'] ='N/A'
t0['quality'] =0
#t0['cosmo'] =9
#t0['phys'] ='N/A'
t0['caltype'] ='s'

t =join(tt,t0,keys='sn') # 10 from jensen 21


for n,i in enumerate(t3['caltype']):
        if i=='c': t3['dist']=-1

w = np.where(t3['caltype']!='c')
          
t4 = vstack([t,t3[w]]) # 17 SBF

t4 = unique(t4,keys='sn',keep='first')

print (t)


t4.remove_column('kmag')
t4.remove_column('ekmag')
t4.remove_column('dsbf')
t4.remove_column('ed')



t.write('../../data/calibrators/calibrators_'+filter+'_sbfcombined.csv', format='ascii.csv', delimiter=',',overwrite=True)
t4.write('../../data/working/'+filter+'_sbfcombined.csv', format='ascii.csv', delimiter=',',overwrite=True)



 
