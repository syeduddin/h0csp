from astropy.io import ascii
from astropy.table import join,vstack, Table, unique
import numpy as np
import sys
import matplotlib.pylab as pl
from scipy.stats import linregress
from collections import Counter
filter = 'B'



j = ascii.read('../../data/calibrators/calibrators_B_sbfj21.csv') # only dist
k = ascii.read('../../data/calibrators/calibrators_sbf_final.csv',format='csv') # only dist



print (np.mean(j['dist']))
print (np.mean(k['dist']))

range = [np.min(k['dist']),np.max(k['dist'])]
pl.hist(j['dist'],histtype='stepfilled',bins=10,range=range,label='Jensen+21',color='r',alpha=.5)
pl.hist(k['dist'],histtype='stepfilled',bins=10,range=range,label='Khetan+21',color='b',alpha=.5)
pl.legend()
pl.xlabel('$Distacne \ moduli \ (\mu)$ [mag]',fontsize=14),pl.ylabel('$Number$',fontsize=14)
pl.grid()

pl.savefig('../../plots/compSBFdist.pdf')
t = join(j,k,keys='sn')

err_int1 = np.std(t['dist_1']) - np.mean(t['edist_1'])
err_int2 = np.std(t['dist_2']) - np.mean(t['edist_2'])

wt1= np.sqrt(1/((t['edist_1']**2)+(err_int1**2))) # weights
wt2= np.sqrt(1/((t['edist_2']**2)+(err_int2**2))) # weights


mean1= np.sum(t['dist_1']*wt1)/np.sum(wt1)
err1 = np.sqrt((1/np.sum(wt1)))

mean2= np.sum(t['dist_2']*wt2)/np.sum(wt2)
err2 = np.sqrt((1/np.sum(wt2)))


print (mean1-mean2)
print (np.mean(t['dist_1']-t['dist_2']))

#print (np.mean(t['BV_1']-t['BV_2']))
print (np.mean(t['edist_1']-t['edist_2']))

print(np.sqrt((err1)**2+(err2)**2))

sys.exit()








pl.errorbar (t['dist_1'],t['dist_2'], xerr=t['edist_1'],yerr=t['edist_2'],fmt='ko',ms=8)
pl.xlabel('Jensen+21'),pl.ylabel('Kehtan+21')
pl.grid()

pl.savefig('../../plots/compSBFmu.pdf')
sys.exit()

pl.xlabel('old'),pl.ylabel('new')
pl.savefig('../../plots/compSBFBmax.pdf')


j =  ascii.read('../../data/calibrators/sbf_jensen21.csv') # only dist
c= ascii.read('../../data/calibrators/calibrators_B_sbfj21.csv') # only dist

print (set(j['gal']).difference(set(c['host'])))
sys.exit()

t1 = ascii.read('../../data/lc/'+filter+'_max_SBF_old.dat')
t2 = ascii.read('../../data/lc/'+filter+'_max_SBF.dat')


         




t =join(t1,t2,keys='name')

pl.errorbar(t['BV_1'],t['BV_2'],xerr=t['eBV_1'],yerr=t['eBV_2'],color='r',fmt='o')
#pl.errorbar(t['Mmax_1'],t['Mmax_2'],xerr=t['eMmax_1'],yerr=t['eMmax_2'],color='r',fmt='o')

pl.xlabel('old'),pl.ylabel('new')
pl.savefig('../../plots/compSBFBmax.pdf')

print (np.mean(t['BV_1']-t['BV_2']))
m, c, r_value, p_value, std_err = linregress(t['Mmax_1'],t['Mmax_2'])
print ('Mmax:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)

m, c, r_value, p_value, std_err = linregress(t['BV_1'],t['BV_2'])
print ('B-V:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)



m, c, r_value, p_value, std_err = linregress(t['st_1'],t['st_2'])
print ('st:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)

