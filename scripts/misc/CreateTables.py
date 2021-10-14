from astropy.io import ascii
from astropy.table import join,vstack,Table
tab = ascii.read('../data/working/B_ceph.csv')

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


data = Table()

data['$Name$']=sn
data['$z_{cmb}$']=zcmb
data['$z_{zhel}$']=zhel
data['$B_{max}$']=mmax
data['$s_{BV}$']=st
data['$B-V$']=bv
data['$M_{host}$']=m_csp
data['$Sample$']=tab['sample']
data['$Subclass$']=tab['subtype']


ascii.write(data[0:10],'../results/table.txt',format='latex',overwrite=True)
