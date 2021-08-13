import matplotlib.pylab as pl
import numpy as np
x=[0,1,2,3,4,5,6,7,8,9,10]
x = np.array(x)

lam = [' ','u','B','g','V','r','i','Y','J','H',' ']





off0=[0, -0.06,  -0.01,  -0.03,  -0.01,  -0.01,  -0.02,   0.00,  -0.01, 0.03,0]

off1=[0,-0.16,  -0.04,  -0.07,  -0.04,  -0.04,  -0.04,  -0.04,  -0.05,  -0.02,0]

off01=[0,-0.136, -0.089, -0.085, -0.076, -0.089, -0.080, -0.122, -0.096, -0.090,0]

eoff0=[0, 0.07,   0.05,   0.06,   0.05,   0.05,   0.05,   0.06,   0.06,   0.07,0]

eoff1=[0, 0.07,   0.05,   0.06,   0.05,   0.05,   0.05,   0.06,   0.06,   0.06,0]

eoff01=[0,0.047,  0.032,  0.027,  0.031,  0.031,  0.033,  0.041,  0.051,  0.048,0]


cb= [0,-.21,-.22,-.22,-.22,-.22,-.17,-.09,-.06,-.03,0]
ecb = [0,0.04,0.04,0.04,0.04,0.04, 0.03,0.02,0.01,0.01, 0]


off=[0,-0.15,  -0.03,  -0.05,  -0.02,  -0.03,  -0.03,  -0.01,  -0.02,  -0.00,0]
eoff=[0,  0.07,   0.05,   0.06,   0.05,   0.05,   0.05,   0.05,   0.06,   0.06,0]
  

cb = np.array(cb)
ecb = np.array(ecb)
pl.figure(1)

pl.xticks(x,lam)

pl.errorbar(x-0.02,off,yerr=eoff,fmt='o',color='b',ms=12,markeredgewidth=2,label='Median',lw=2)
pl.errorbar(x+0.05,off0,yerr=eoff0,fmt='o',color='r',ms=12,markeredgewidth=2,label='10.0',lw=2)
pl.errorbar(x-0.05,off1,yerr=eoff1,fmt='o',color='g',ms=12,markeredgewidth=2,label='10.5',lw=2)



pl.plot(x,cb,'k-',label='Simulated dust offsets',lw=2)
pl.fill_between(x, cb-ecb, cb+ecb, facecolor='k', alpha=0.3)

pl.xlabel('$Bands$',fontsize=20), pl.ylabel('$\Delta_{HR} \ (mag)$',fontsize=20)


pl.tick_params(axis='both', labelsize=14)
#pl.legend(['$Bin \ split \ = \ 10.48 \ M_{\odot}$','$Bin \ split \ = \ 10 \ M_{\odot}$','$Simulated \ dust \ offsets$'],loc='lower right',numpoints=1)
pl.legend(loc='lower right',numpoints=1)

pl.grid()
pl.xlim(.9,9.1); pl.ylim(-.3,.1)

pl.savefig('/Users/suddin/Dropbox/CSP/plots/offsim.pdf')

pl.show()
