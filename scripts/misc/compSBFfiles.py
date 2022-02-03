from astropy.io import ascii
from astropy.table import join
import matplotlib.pylab as pl
import os,sys
import numpy as np
from scipy.stats import linregress
import matplotlib.gridspec as gridspec
import linmix

dir ='../../data/working/'

old = ascii.read(dir+'B_sbf_old.csv',format='csv')
new = ascii.read(dir+'B_sbf.csv')



table = join(old,new,keys='sn')


print (np.mean(table['BV_1']-table['BV_2']))
print (np.mean(table['st_1']-table['st_2']))
print (np.mean(table['Mmax_1']-table['Mmax_2']))

m, c, r_value, p_value, std_err = linregress(table['Mmax_1'],table['Mmax_2'])
print ('Mmax:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)

m, c, r_value, p_value, std_err = linregress(table['BV_1'],table['BV_2'])
print ('B-V:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)



m, c, r_value, p_value, std_err = linregress(table['st_1'],table['st_2'])
print ('st:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)

sys.exit()

result = linregress(tab['BV'],tab['color'])

print(result.slope, result.stderr)

#print(result.intercept, result.intercept_stderr)


#lm = linmix.LinMix(table['BV'], table['color'], table['eBV'], table['e_color'], K=2)
#lm.run_mcmc(silent=True)
#print('%6.3f'%(np.mean(lm.chain['beta'])))
#print('%6.3f'%np.std(lm.chain['beta']))


#sys.exit()

pl.figure(figsize=(20,10))

pl.subplot(1,3,1,aspect='equal')
pl.subplots_adjust(wspace=0.0,hspace=0.0)
gs=gridspec.GridSpec(2,1,height_ratios=[3,2])
ax1=pl.subplot(gs[0])
pl.setp(ax1.get_xticklabels(),visible=False)

m, c, r_value, p_value, std_err = linregress(table['BV'],table['color'])
print ('B-V:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)
x = np.arange(-.2,.5,.1)
y = m*x+c
pl.plot(x,y,'r',lw=2)
pl.errorbar(table['BV'],table['color'],xerr=table['eBV'],yerr=table['e_color'], color='r',fmt='o',ms=16)
pl.errorbar(tab['BV'],tab['color'],xerr=tab['eBV'],yerr=tab['e_color'], color='b',fmt='o',ms=16)

pl.title(r'$B-V$', fontsize=20),pl.ylabel(r'$Khetan \ et \ al. \ (2021)$', fontsize=20)
pl.grid()
pl.xticks(np.arange(-.2,.5,.1))
pl.xlim(-.2,.5,.1)

m, c, r_value, p_value, std_err = linregress(tab['BV'],tab['color'])
print ('B-V:','%.2f'%m,'%.2f'%c,'%.2f'%std_err)
x = np.arange(-.2,.5,.1)
y = m*x+c
pl.plot(x,y,'b',lw=2)


ax2=pl.subplot(gs[1])
pl.setp(ax2.get_xticklabels(),visible=True)
pl.errorbar(table['BV'],table['BV']-table['color'],xerr=table['eBV'],yerr=np.sqrt((table['eBV']**2)+(table['e_color']**2)), color='r',fmt='o',ms=16)

pl.errorbar(tab['BV'],tab['BV']-tab['color'],xerr=tab['eBV'],yerr=np.sqrt((tab['eBV']**2)+(tab['e_color']**2)), color='b',fmt='o',ms=16)


pl.xlabel(r'$This \ Work$', fontsize=20),pl.ylabel(r'$Residual$', fontsize=20),
pl.axhline(0,color='k',lw=2)
pl.grid()
pl.xticks(np.arange(-.2,.5,.1))
pl.yticks(np.arange(-.2,.2,.1))

pl.xlim(-.2,.5,.1)

pl.savefig('compBV.pdf')
os.system('open compBV.pdf')


