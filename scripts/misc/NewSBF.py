from astropy.io import ascii
from astropy.table import join,vstack, Table
import numpy as np
t1 = ascii.read('../data/calibrators/calibrators_sbf_new.csv')
t2 = ascii.read('../data/CSPI+II_max/B_max_SBF.csv')

t2.remove_column('dist')
t2.remove_column('edist')

t3 = Table()

t3['sn'] = t1['sn']
t3['dist'] = t1['dist']
t3['edist'] = t1['edist']
t3['host'] = t1['host']
t3['m'] = t1['m']
t3['ml'] = t1['m']-t1['ml']
t3['mu'] = t1['mu']+t1['m']

t =join(t2,t3,keys='sn')
print (t2)
print (t)

table =Table()

table['sn']=t['sn']
table['zhel']=t['zhel']
table['zcmb']=t['zcmb']
table['st']=t['st']
table['est']=t['est']
table['Mmax']=t['Mmax']
table['eMmax']=t['eMmax']
table['t0']=t['t0']
table['EBVmw']=t['EBVmw']
table['eEBVmw']=t['eEBVmw']
table['covMs']=t['covMs']
table['sample']='N/A'
table['subtype']='N/A'
table['quality']=0
table['cosmo']='N/A'
table['phys']='N/A'
table['BV']=t['BV']
table['eBV']=t['eBV']
table['covBV_M']=t['covBV_M']
table['ml']=t['ml']
table['m']=t['m']
table['mu']=t['mu']
table['host']=t1['host']
table['dist']=t['dist']
table['edist']=t['edist']
table['caltype']='s'

print (table)
table.write('../data/calibrators/calibrators_sbf_final.csv', format='ascii.csv', delimiter=',',overwrite=True)
