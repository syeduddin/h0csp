from astropy.io import ascii
from astropy.table import join,vstack, Table, unique
import numpy as np
import sys
import matplotlib.pylab as pl
from scipy.stats import linregress
from collections import Counter
filter = sys.argv[1]

t1 = ascii.read('../../data/calibrators/calibrators_sbf_final.csv')
t2 = ascii.read('../../data/lc/'+filter+'_max_SBF.dat')
t4 = ascii.read('../../data/working/'+filter+'_ceph.csv')

for n,i in enumerate(t4['caltype']):
        if i=='c': t4['dist']=-1


         
t2.rename_column('name','sn')



t3 = Table()

t3['sn'] = t1['sn']
t3['dist'] = t1['dist']
t3['edist'] = t1['edist']
t3['host'] = t1['host']
t3['m'] = t1['m']
t3['ml'] = t1['m']-t1['ml']
t3['mu'] = t1['mu']+t1['m']

t =join(t1,t2,keys='sn')

#pl.errorbar(t['BV_1'],t['BV_2'],xerr=t['eBV_1'],yerr=t['eBV_2'],color='r',fmt='o')
pl.errorbar(t['Mmax_1'],t['Mmax_2'],xerr=t['eMmax_1'],yerr=t['eMmax_2'],color='r',fmt='o')

pl.xlabel('old'),pl.ylabel('new')
pl.savefig('../../plots/compSBFBmax.pdf')

#print (np.mean(t['BV_1']-t['BV_2']))
m, c, r_value, p_value, std_err = linregress(t['Mmax_1'],t['Mmax_2'])
print ('Mmax:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)

m, c, r_value, p_value, std_err = linregress(t['BV_1'],t['BV_2'])
print ('B-V:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)



m, c, r_value, p_value, std_err = linregress(t['st_1'],t['st_2'])
print ('st:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)

table =Table()

table['sn']=t['sn']
table['zhel']=t['zhel_2']
table['zcmb']=t['zcmb_2']
table['st']=t['st_2']
table['est']=t['est_2']
table['Mmax']=t['Mmax_2']
table['eMmax']=t['eMmax_2']
table['t0']=t['t0_2']
table['EBVmw']=t['EBVmw_2']
table['eEBVmw']=t['eEBVmw_2']
table['covMs']=t['covMs_2']
table['sample']='N/A'
table['subtype']='N/A'
table['quality']=0
table['cosmo']=9
table['phys']=9
table['BV']=t['BV_2']
table['eBV']=t['eBV_2']
table['covBV_M']=t['covBV_M_2']
table['ml']=t['ml']
table['m']=t['m']
table['mu']=t['mu']
table['host']=t['host']
table['dist']=t['dist']
table['edist']=t['edist']
table['caltype']='s'


final = vstack([t4,table])
final = unique(final,keys='sn',keep='last')


table.write('../../data/calibrators/calibrators_sbf_'+filter+'.csv', format='ascii.csv', delimiter=',',overwrite=True)

#final.write('../../data/working/'+filter+'_sbf.csv', format='ascii.csv', delimiter=',',overwrite=True)
