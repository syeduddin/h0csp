from astropy.table import join,vstack,unique
from astropy.io import ascii

tab1 = ascii.read('data/khetan/lowz.csv')
tab2 = ascii.read('data/khetan/cal.csv')

# tab3 has no 'sn1' so just copy 'sn'
#tab3['sn1'] = tab3['sn']
# tab2 has no host, dist, edist or caltype, no make invalid entries
tab2['redshift']=-1.0
tab1['host']='N/A'
tab1['mu_ceph'] = -1.0
tab1['e_mu'] = -1.0


# now that tab2 and tab3 have consistent columns, add them together
#masses = vstack([tab2,tab3])

# Now we join tab1 to masses using 'sn' as the join key
#t1 = join(tab1,masses, keys='sn')
#t1.remove_column('sn1')
# Next, join tab1 to masses using 'sn1' as the join key
#tab1.rename_column('sn','sn1')
#t2 = join(tab1,masses, keys='sn1')
#t2.remove_column('sn')
#t2.rename_column('sn1','sn')

# Now we put both the sn-indexed and sn1-indexed tables together
bigt = vstack([tab1,tab2])

#Lastly, get rid of the repeated entries
newt = unique(bigt, keys='SN_name')
newt.write('data/khetan/Shoes.dat', format='ascii.fixed_width', delimiter=' ')
#newt.write('data/B_sbf.csv', format='ascii.csv', delimiter=',')
