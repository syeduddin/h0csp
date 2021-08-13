from astropy.table import join,vstack,unique
from astropy.io import ascii
import sys
filter = sys.argv[1]
tab1 = ascii.read('../CSPI+II_max/'+filter+'_max.dat')
tab1.rename_column('name','sn')      # for consistency
tab1.remove_column('dist')           # avoid collision
tab2 = ascii.read('../CSPHostMass.csv')
tab3 = ascii.read('cal_ceph.csv')

# tab3 has no 'sn1' so just copy 'sn'
tab3['sn1'] = tab3['SN_name']
# tab2 has no host, dist, edist or caltype, no make invalid entries
tab2['host'] = 'N/A'
tab2['dist'] = -1
tab2['edist'] = -1
tab2['caltype'] = 'none'

tab3.rename_column('log_Mass','m')
tab3.remove_column('')

# now that tab2 and tab3 have consistent columns, add them together
masses = vstack([tab2,tab3])
print masses
sys.exit()
# Now we join tab1 to masses using 'sn' as the join key
t1 = join(tab1,masses, keys='sn')
t1.remove_column('sn1')
# Next, join tab1 to masses using 'sn1' as the join key
tab1.rename_column('sn','sn1')
t2 = join(tab1,masses, keys='sn1')
t2.remove_column('sn')
t2.rename_column('sn1','sn')

# Now we put both the sn-indexed and sn1-indexed tables together
bigt = vstack([t1,t2])

#Lastly, get rid of the repeated entries
newt = unique(bigt, keys='sn')


#newt.write('data/'+filter+'_ceph.dat', format='ascii.fixed_width', delimiter=' ',overwrite=True)
#newt.write('test.csv', format='ascii.csv', delimiter=' ')
