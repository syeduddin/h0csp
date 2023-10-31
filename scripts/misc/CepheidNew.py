from astropy.io import ascii
from astropy.table import join,vstack,unique
import sys

filter = sys.argv[1]
old = ascii.read('../../data/working/'+filter+'_ceph_update3.csv')
ceph22 = ascii.read('../../data/calibrators/cepheid_riess22.dat',format='basic')

old.rename_column('sn','SN')      # for consistency
old.remove_column('host')
old.remove_column('caltype')

SN=[]
name = ceph22['SN']
for i in range(len(name)):
        SN.append('SN'+name[i])
    

ceph22['SN']=SN
ceph22.remove_column('N')
#ceph22.remove_column('Host')
ceph22.remove_column('${m}_{B,i}^{0}$')
ceph22.remove_column('sigma1')
ceph22.remove_column('mu_Ceph')
ceph22.remove_column('sigma2')
ceph22.remove_column('${M}_{B,i}^{0}$')
ceph22.remove_column('sigma3')
ceph22.remove_column('R')
#ceph22.rename_column('Host','host')


t = join(old,ceph22,keys='SN')
t.rename_column('dist_2','dist')
t.rename_column('edist_2','edist')
t.remove_column('dist_1')
t.remove_column('edist_1')
t['caltype']='c'

table = vstack([old,t])

table =unique(table,keys='SN',keep='last')

table.rename_column('SN','sn')      # for consistency
table.rename_column('Host','host')      # for consistency


print (len(table))
print (len(t))
table.write('../../data/working/'+filter+'_ceph_update3.csv', format='ascii.csv', delimiter=',',overwrite=True)
