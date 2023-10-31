from astropy.io import ascii
from astropy.table import join,vstack,Table
import sys
import numpy 


filter = sys.argv[1]

host = ascii.read('../../data/hosts/V_ceph_update2+host.csv')

lc = ascii.read('../../data/working/'+filter+'_ceph_update3.csv')

tab = join(host,lc,keys='sn')



st = tab['st_1']
est = tab['est_1']
zhel = tab['zhel_1']
zcmb = tab['zcmb_1']
mmax = tab['Mmax_1']
emmax = tab['eMmax_1']
bv = tab['BV_1']
ebv = tab['eBV_1']
m_csp = tab['m_1']
dist = tab['dist_2']
edist = tab['edist_1']
c_ms = tab['covMs_1']
c_mbv = tab['covBV_M_1']
c_sbv = tab['covBVs']

sn = tab['sn']
host = tab['host_1']



w = numpy.where(zcmb>0)

data = Table()
data['Name']=sn[w]
data['Host'] = host[w]
data['z_{cmb}']=zcmb[w].round(4)
data['z_{zhel}']=zhel[w].round(5)

data['B_{max}']=mmax[w].round(5)
data['eB_{max}']=emmax[w].round(5)
data['s_{BV}']=st[w]
data['es_{BV}']=est[w].round(5)
data['B-V']=bv[w].round(5)
data['eB-V']=ebv[w].round(5)
data['cov_{ms}'] =c_ms[w].round(5)
data['cov_{mbv}'] =c_mbv[w].round(5)
data['cov_{sbv}'] =c_sbv[w].round(5)

data['M_{host}']=m_csp[w].round(5)
data['Sample']=tab['sample_1'][w]
data['Subclass']=tab['subtype_1'][w]
data['Distant']= dist[w]
data['eDistant']= edist[w]

print (data)
ascii.write(data,'../../data/appendix/'+filter+'_update3.csv',format='csv',overwrite=True)
