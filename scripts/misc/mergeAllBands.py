from astropy.table import join,hstack,vstack
from astropy.io import ascii
import sys
import numpy as np

u = ascii.read('../data/working/u_ceph.csv')
B = ascii.read('../data/working/B_ceph.csv')
g = ascii.read('../data/working/g_ceph.csv')
V = ascii.read('../data/working/V_ceph.csv')
r = ascii.read('../data/working/r_ceph.csv')
i = ascii.read('../data/working/i_ceph.csv')
Y = ascii.read('../data/working/Y_ceph.csv')
J = ascii.read('../data/working/J_ceph.csv')
H = ascii.read('../data/working/H_ceph.csv')

data =join(H,Y,keys='sn')

w=data['caltype_2']=='c'
print (data['st_2'][w],data['st_1'][w])
sys.exit()

table = hstack([u,B,g,V,r,i,Y,J,H])
#w = np.where((table['caltype_7']=='none') | (table['caltype_9']=='c'))

#table= (table[w])



table.write('../data/working/Allbands_ceph.csv', format='ascii.csv', delimiter=',',overwrite=True)
#newt.write('test.csv', format='ascii.csv', delimiter=' ')

tab = ascii.read('../data/working/Allbands_ceph.csv')

for i in range(1,10,1):
    w =tab['caltype_'+str(i)]=='c'
    print (len(tab['caltype_'+str(i)][w]), len(tab['caltype_'+str(i)]))
    
