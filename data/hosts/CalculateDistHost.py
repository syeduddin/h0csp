import numpy as np
from astropy.io import ascii
from astropy.table import join, vstack,unique, Table
import cosmolopy.distance as cd
cosmo = {'omega_M_0' : .3, 'omega_lambda_0' : .7, 'h' : 0.70}
cosmo1 = cd.set_omega_k_0(cosmo)


c1 = ascii.read('../../data/hosts/CSP1_SNIa_coord.csv')
c2 = ascii.read('../../data/hosts/CSP2_SN_host_coord.csv')
b = ascii.read('../../data/lc/Spreadsheet.csv')

c1.remove_column('zc')
c1.remove_column('gal')
c1.rename_column('ra','snra')
c1.rename_column('dec','sndec')
c = vstack([c1,c2])

b.replace_name('Name','sn')
tab1 = join(c,b,keys='sn')
b.remove_column('sn')
b.replace_name('Name(P19)','sn')

tab2 = join(c,b,keys='sn')
tab = vstack([tab1,tab2])

#w = np.where((tab['Sub-type']!='Ia-02cx') & ((tab['Sub-type']!='Ia-SC')))
#tab=tab[w]
print (tab)
sra=tab['snra']
sdec=tab['sndec']
hra=tab['hra']
hdec=tab['hdec']
z =tab['zcmb']

impact = np.sqrt((((sra-hra)*3600)*np.cos(sdec*np.pi/180))**2 + ((sdec-hdec)*3600)**2)

proj=(cd.angular_diameter_distance(np.abs(z),**cosmo))*(impact/206.264806) # projected distance in kpc

table = Table()
table['sn']= tab['sn']
table['proj'] = proj

table.write('../../data/hosts/proj_dist.csv',format='ascii.csv', delimiter=',',overwrite=True)
