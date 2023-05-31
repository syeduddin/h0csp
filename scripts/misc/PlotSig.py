import matplotlib.pylab as pl
import numpy as np





x=[0,1,2,3,4,5,6,7,8,9,10]
x = np.array(x)

lam = [' ','u','B','g','V','r','i','Y','J','H',' ']

sig =[0,0.219, 0.172, 0.14, 0.176, 0.162, 0.147, 0.117, 0.172, 0.138,0]
esig=[0,0.0205, 0.013, 0.0165, 0.013, 0.0115, 0.0115, 0.01049, 0.0144, 0.013,0]

sig91 = [0,0.195, 0.152, 0.123, 0.155, 0.14, 0.119, 0.098, 0.139, 0.115,0]
esig91 = [0,0.022, 0.013, 0.018, 0.013, 0.0115, 0.010, 0.010, 0.0125, 0.015,0]

sigz = [0,0.241, 0.189, 0.168, 0.192, 0.176, 0.165, 0.134, 0.186, 0.163,0]
esigz = [0,0.0165, 0.018, 0.0135, 0.011, 0.0104, 0.0125, 0.01, 0.0144, 0.0125,0]

sigz1 = [0,0.218, 0.169, 0.142, 0.175, 0.159, 0.144, 0.107, 0.162, 0.127,0]
esigz1 = [0,0.018, 0.0115, 0.0135, 0.011, 0.01, 0.009, 0.01, 0.0115, 0.0125,0]



pl.figure(1)

pl.xticks(x,lam)

pl.errorbar(x-0.02,sig,yerr=esig,fmt='o',color='b',ms=12,markeredgewidth=2,label='No cuts',lw=2)
pl.errorbar(x,sigz,yerr=esigz,fmt='s',color='m',ms=12,markeredgewidth=2,label='z>0.01',lw=2)
#pl.errorbar(x+0.04,sigz1,yerr=esigz1,fmt='s',color='g',ms=10,markeredgewidth=2,label='z>0.01 (fixed Vpec)',lw=2)


pl.errorbar(x+0.05,sig91,yerr=esig91,fmt='d',color='r',ms=12,markeredgewidth=2,label='Excluding 91T & 91bg',lw=2)





pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$\sigma_{int} \ (mag)$',fontsize=14)


#pl.tick_params(axis='both', labelsize=14)
#pl.legend(['$Bin \ split \ = \ 10.48 \ M_{\odot}$','$Bin \ split \ = \ 10 \ M_{\odot}$','$Simulated \ dust \ offsets$'],loc='lower right',numpoints=1)
pl.legend(loc='upper right',numpoints=1)

pl.grid()
pl.xlim(.8,9.2); pl.ylim(0.08,.25)

pl.tight_layout()

pl.savefig('../../plots/sigma_int.pdf')
#pl.show()
