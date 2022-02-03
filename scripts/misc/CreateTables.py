from astropy.io import ascii
from astropy.table import join,vstack,Table

tab = ascii.read('../../data/working/B_ceph.csv')
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
data['$z_{cmb}$']=zcmb[w]
data['$z_{zhel}$']=zhel[w]
data['$B_{max}$']=mmax[w]
data['$eB_{max}$']=emmax[w]
data['$s_{BV}$']=st[w]
data['$es_{BV}$']=est[w]
data['$B-V$']=bv[w]
data['$eB-V$']=ebv[w]
data['$cov_{ms}$'] =c_ms[w]
data['$cov_{mbv}$'] =c_mbv[w]
data['$M_{host}$']=m_csp[w]
#data['$Sample$']=tab['sample'][w]
#data['$Subclass$']=tab['subtype'][w]
data['$Distant$']= dist[w]
data['$eDistant$']= edist[w]



ascii.write(data,'../../results/B_ceph.txt',format='latex',overwrite=True)
