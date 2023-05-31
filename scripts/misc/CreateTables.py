from astropy.io import ascii
from astropy.table import join,vstack,Table

tab = ascii.read('../../data/working/B_sbf_update2.csv')
import numpy as np

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

w = np.where(dist>1.0)

data = Table()
data['$Name$']=sn[w]
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
#data['$Sample$']=tab['sample'][w]
#data['$Subclass$']=tab['subtype'][w]
data['$Distant$']= dist[w]
data['$eDistant$']= edist[w]


print (data)
ascii.write(data,'../../results/B_sbfk21cal_update2.txt',format='latex',overwrite=True,)
