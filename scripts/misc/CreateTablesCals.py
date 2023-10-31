from astropy.io import ascii
from astropy.table import join,vstack,Table
import sys
import numpy as np


filter = sys.argv[1]


tab = ascii.read('../../data/working/'+filter+'_sbfj21_update3.csv')




st = tab['st']
est = tab['est']
zhel = tab['zhel']
zcmb = tab['zcmb']
mmax = tab['Mmax']
emmax = tab['eMmax']
bv = tab['BV']
ebv = tab['eBV']
m_csp = tab['m']
dist = tab['dist']
edist = tab['edist']
c_ms = tab['covMs']
c_mbv = tab['covBV_M']
sn = tab['sn']
host = tab['host']



w = np.where(dist>1.0)

data = Table()
data['$Name$']=sn[w]
data['$Host$'] = host[w]
data['$z_{cmb}$']=zcmb[w].round(4)
data['$z_{zhel}$']=zhel[w].round(5)

data['$B_{max}$']=mmax[w].round(5)
data['$eB_{max}$']=emmax[w].round(5)
data['$s_{BV}$']=st[w]
data['$es_{BV}$']=est[w].round(5)
data['$B-V$']=bv[w].round(5)
data['$eB-V$']=ebv[w].round(5)
data['$cov_{ms}$'] =c_ms[w].round(5)
data['$cov_{mbv}$'] =c_mbv[w].round(5)
data['$M_{host}$']=m_csp[w]
#data['$Sample$']=tab['sample_1'][w]
#data['$Subclass$']=tab['subtype_1'][w]
data['$Distant$']= dist[w]
data['$eDistant$']= edist[w]

print (data)
print (len(data))

ascii.write(data,'../../data/appendix/'+filter+'_sbf21.csv',format='csv',overwrite=True)
