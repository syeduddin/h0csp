import matplotlib.pylab as pl
import numpy as np





x=[0,1,2,3,4,5,6,7,8,9,10]
x = np.array(x)

lam = [' ','$u$','$B$','$g$','$V$','$r$','$i$','$Y$','$J$','$H$',' ']


#Sim
cb= [0,-.21,-.22,-.22,-.22,-.22,-.17,-.09,-.06,-.03,0]
ecb = [0,0.04,0.04,0.04,0.04,0.04, 0.03,0.02,0.01,0.01, 0]

#CSPI
off1=[0,-0.147, -0.089, -0.076, -0.074, -0.087, -0.085, -0.137, -0.090, -0.093,0]
eoff1=[0,0.044,  0.031,  0.028,  0.030,  0.029,  0.032,  0.038,  0.046,  0.043,0]

off1=[0,-0.03,  -0.02,  -0.01,  -0.02,  -0.03,  -0.03,  -0.08,  -0.08,  -0.04,0]
eoff1=[0, 0.09,   0.08,   0.08,   0.08,   0.08,   0.08,   0.08,   0.09,   0.09,0]

#CSPII
#off2 = [0, -0.13,  -0.03,  -0.07,  -0.04,  -0.03,  -0.04,  -0.04,  -0.02,  -0.03,0]old
#eoff2=[0,0.10,   0.06,   0.08,   0.06,   0.06,   0.06,   0.07,   0.08,   0.09,0]

off2=[0, -0.03,  -0.10,  -0.05,  -0.10,  -0.09,  -0.09,  -0.05,  -0.03,  -0.03,0]
eoff2=[0, 0.10,   0.07,   0.10,   0.07,   0.07,   0.07,   0.07,   0.08,   0.09,0]

#full
#off=[0,-0.14,  -0.03,  -0.05,  -0.03,  -0.02,  -0.03,   0.01,  -0.02,   0.00,0]old
#eoff=[0, 0.07,   0.05,   0.06,   0.05,   0.05,   0.05,   0.05,   0.06,   0.06,0]


off=[0,-0.08,  -0.06,  -0.05,  -0.06,  -0.05,  -0.05,  -0.01,  -0.03,   0.00,0]
eoff=[0,  0.07,   0.05,   0.06,   0.05,   0.05,   0.05,   0.05,   0.06,   0.06,0]

cb = np.array(cb)
ecb = np.array(ecb)

pl.figure(1)

pl.xticks(x,lam)

pl.errorbar(x-0.1,off1,yerr=eoff1,fmt='s',color='#377eb8',ms=12,markeredgewidth=2,label='$CSP-I$',lw=2)
#pl.errorbar(x+0.2,off1,yerr=eoff1,fmt='*',color='#ff7f00',ms=14,label='$CSP-I \  (Uddin \ 2020)$',lw=2)
pl.errorbar(x+0.1,off2,yerr=eoff2,fmt='d',color='r',ms=12,markeredgewidth=2,label='$CSP-II$',lw=2)

pl.errorbar(x-0.2,off,yerr=eoff,fmt='o',color='b',ms=12,markeredgewidth=2,label='$CSP-I \ & \ II$',lw=2)


pl.plot(x,cb,'k-',label='$Simulated \ dust \ offsets$',lw=2)
pl.fill_between(x, cb-ecb, cb+ecb, facecolor='k', alpha=0.3)

pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$\Delta_{HR} \ (mag)$',fontsize=14)


#pl.tick_params(axis='both', labelsize=14)
#pl.legend(['$Bin \ split \ = \ 10.48 \ M_{\odot}$','$Bin \ split \ = \ 10 \ M_{\odot}$','$Simulated \ dust \ offsets$'],loc='lower right',numpoints=1)
pl.legend(loc='lower right',numpoints=1,fontsize=10)

pl.grid()
pl.xlim(.5,9.5); pl.ylim(-.3,.1)

pl.savefig('../../plots/offsim.pdf')
pl.tight_layout()
#pl.show()
