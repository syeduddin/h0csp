import matplotlib.pylab as pl
import numpy as np





x=[0,1,2,3,4,5,6,7,8,9,10]
x = np.array(x)

lam = [' ','$u$','$B$','$g$','$V$','$r$','$i$','$Y$','$J$','$H$',' ']


#Sim
cb= [0,-.21,-.22,-.22,-.22,-.22,-.17,-.09,-.06,-.03,0]
ecb = [0,0.04,0.04,0.04,0.04,0.04, 0.03,0.02,0.01,0.01, 0]

#CSPI
off1=[0,-0.09,  -0.03,  -0.02,  -0.03,  -0.05,  -0.05,  -0.10,  -0.11,  -0.06,0]
eoff1=[0,0.09,   0.08,   0.08,   0.08,   0.08,   0.08,   0.08,   0.09,   0.09,0]

#slopes
#-0.035,  0.006,  0.017,  0.005, -0.009, -0.014, -0.055, -0.073, -0.005
# 0.041,  0.028,  0.028,  0.028,  0.026,  0.025,  0.025,  0.036,  0.030

#CSPII
off2=[0, -0.11,  -0.12,  -0.07,  -0.11,  -0.11,  -0.11,  -0.08,  -0.07,  -0.06,0]
eoff2=[0, 0.11,   0.07,   0.10,   0.07,   0.07,   0.07,   0.07,   0.08,   0.09,0]
#slopes
#-0.062, -0.043, -0.035, -0.044, -0.044, -0.050, -0.034, -0.045, -0.036
# 0.028,  0.015,  0.023,  0.015,  0.014,  0.012,  0.012,  0.017,  0.019


#full

off=[0,-0.15,  -0.07,  -0.06,  -0.07,  -0.07,  -0.07,  -0.04,  -0.07,  -0.03,0]
eoff=[0, 0.07,   0.05,   0.06,   0.05,   0.05,   0.05,   0.05,   0.06,   0.06,0]
#slopes
#-0.066, -0.038, -0.024, -0.039, -0.038, -0.043, -0.019, -0.039, -0.014
#0.021,  0.012,  0.016,  0.012,  0.011,  0.010,  0.010,  0.015,  0.015

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
pl.xticks(fontsize=14),pl.yticks(fontsize=14)
#pl.tick_params(axis='both', labelsize=14)
#pl.legend(['$Bin \ split \ = \ 10.48 \ M_{\odot}$','$Bin \ split \ = \ 10 \ M_{\odot}$','$Simulated \ dust \ offsets$'],loc='lower right',numpoints=1)
pl.legend(loc='lower right',numpoints=1,fontsize=10)
pl.grid()
pl.xlim(.5,9.5); pl.ylim(-.3,.1)
pl.tight_layout()

pl.savefig('../../plots/offsim.pdf')
#pl.show()

pl.figure(2)
pl.xticks(x,lam)

trgb = [0,67.201, 69.334, 69.33, 69.328, 68.232, 68.012, 72.147, 72.352, 72.406,0]
etrgb=[0,1.0990000000000002, 0.7545, 1.309, 0.7184999999999999, 0.6855, 0.6595, 1.165, 0.8285, 0.804,0]

ceph = [0,75.109, 73.676, 78.004, 73.604, 72.664, 72.558, 74.313, 73.903, 74.902,0]
eceph=[0,1.2625000000000002, 0.7464999999999999, 1.2730000000000001, 0.699, 0.7095, 0.6605000000000001, 1.245, 0.9694999999999999, 0.8795,0]

sbf=[0,69.091, 72.623, 72.919, 72.535, 72.063, 71.956, 71.799, 70.152, 69.126,0]
esbf=[0,1.443, 0.911, 1.399, 0.87, 0.802, 0.781, 1.532, 1.2995, 1.042,0]

pl.errorbar(x-.1,ceph,yerr=eceph,fmt='o',color='b',ms=12,markeredgewidth=2,label='$Cepheid$',lw=2)

pl.errorbar(x,trgb,yerr=etrgb,fmt='d',color='r',ms=12,markeredgewidth=2,label='$TRGB$',lw=2)
pl.errorbar(x+.1,sbf,yerr=esbf,fmt='s',color='#377eb8',ms=12,markeredgewidth=2,label='$TRGB$',lw=2)
pl.legend(loc='upper right',numpoints=1,fontsize=10)
pl.grid()
pl.xlim(.5,9.5),pl.ylim(65,80)
pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)

pl.savefig('../../plots/h0bands.pdf')
