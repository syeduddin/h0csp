import matplotlib.pylab as pl
import numpy as np

x=[0,1,2,3,4,5,6,7,8,9,10]
x = np.array(x)

lam = [' ','$u$','$B$','$g$','$V$','$r$','$i$','$Y$','$J$','$H$',' ']


#Sim
cb= [0,-.21,-.22,-.22,-.22,-.22,-.17,-.09,-.06,-.03,0]
ecb = [0,0.04,0.04,0.04,0.04,0.04, 0.03,0.02,0.01,0.01, 0]

#CSPI
off1=[0,-0.10,  -0.04,  -0.02,  -0.04,  -0.05,  -0.05,  -0.10,  -0.10,  -0.08,0]
eoff1=[0,0.09,   0.08,   0.08,   0.08,   0.08,   0.08,   0.08,   0.09,   0.09,0]

#slopes
#-0.035,  0.006,  0.017,  0.005, -0.009, -0.014, -0.055, -0.073, -0.005
# 0.041,  0.028,  0.028,  0.028,  0.026,  0.025,  0.025,  0.036,  0.030

#CSPII
off2=[0, -0.14,  -0.11,  -0.09,  -0.11,  -0.11,  -0.11,  -0.10,  -0.08,  -0.05,0]
eoff2=[0, 0.11,   0.07,   0.10,   0.07,   0.07,   0.07,   0.07,   0.08,   0.09,0]
#slopes
#-0.062, -0.043, -0.035, -0.044, -0.044, -0.050, -0.034, -0.045, -0.036
# 0.028,  0.015,  0.023,  0.015,  0.014,  0.012,  0.012,  0.017,  0.019


#full

off=[0,-0.13,  -0.08,  -0.08,  -0.08,  -0.07,  -0.07,  -0.07,  -0.04,  -0.01,0]
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

trgb = [0,68.834, 69.879, 69.57, 69.598, 68.713, 69.26, 70.776, 72.591, 71.066,0]
etrgb=[0,1.0775, 0.75550, 1.1995, 0.714, 0.6105, 0.608, 1.0975, 0.8115, 0.7335,0]

ceph = [0,75.576, 73.385, 75.484, 73.442, 72.407, 73.414, 73.434, 75.73, 74.578,0]
eceph=[0,1.168, 0.7315, 1.0415, 0.713, 0.667, 0.668, 0.945, 0.827, 0.7955,0]

sbf=[0,69.779, 72.732, 73.496, 72.586, 72.206, 72.701, 70.136, 70.769, 68.748,0]
esbf=[0,1.479, 0.928, 1.3904, 0.87, 0.8385, 0.7865, 1.233, 1.0955, 0.98,0]

pl.errorbar(x-.15,ceph,yerr=eceph,fmt='o',color='b',ms=12,markeredgewidth=2,label='$Cepheid$',lw=2)

pl.errorbar(x,trgb,yerr=etrgb,fmt='d',color='g',ms=14,markeredgewidth=2,label='$TRGB$',lw=2)
#pl.errorbar(x+.15,sbf,yerr=esbf,fmt='s',color='#377eb8',ms=12,markeredgewidth=2,label='$SBF \ (combined)$',lw=2)
pl.errorbar(x+.15,sbf,yerr=esbf,fmt='s',color='r',ms=10,markeredgewidth=2,label='$SBF \ (combined)$',lw=2)
pl.legend(loc='upper right',numpoints=1,fontsize=10)
pl.grid()
pl.xlim(.5,9.5),pl.ylim(65,80)
pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)

pl.savefig('../../plots/h0bands.pdf')



pl.figure(3)
pl.xticks(x,lam)

trgb1 = [0,69.981, 71.239, 70.679, 71.082, 70.443, 69.93, 70.669, 70.024, 69.956,0]
etrgb1=[0,1.269, 0.8745, 1.4575, 0.8775, 0.8375, 0.829, 1.225, 0.835, 0.8855,0]

trgb2 = [0,65.88, 68.79, 69.053, 68.538, 68.295, 67.95, 72.576, 71.3, 72.974,0]
etrgb2=[0,1.3695, 0.8214, 1.4845, 0.7969, 0.7665, 0.739, 1.234, 0.9584, 0.9544,0]

wendy19 = [0,0, 69.90, 0, 0, 0, 69.00, 0, 69.56, 69.21,0]
ewendy19 = [0,0, 1.41,0, 0, 0, 1.30,0, 1.36, 1.35, 0]

trgbz=[0,67.762, 69.785, 69.512, 69.514, 69.139, 69.002, 71.912, 70.819, 71.067,0]
etrgbz = [0,1.081, 0.681, 1.21, 0.67599, 0.71, 0.605, 1.124, 0.761, 0.725,0]

trgbst = [0,67.411, 69.545, 69.552, 69.293, 68.837, 68.604, 71.76, 69.973, 70.837,0]
etrgbst = [0,1.091, 0.7075, 1.311, 0.6875, 0.6705, 0.655, 1.141, 0.735, 0.793,0]

trgbbv = [0,67.89, 70.005, 69.011, 69.705, 69.351, 68.958, 71.785, 70.113, 70.715,0]
etrgbbv = [0, 1.047, 0.6915, 1.224, 0.692, 0.7015, 0.6345, 1.1055, 0.754, 0.7835,0]

trgbt0 = [0,67.514, 69.482, 69.629, 69.407, 69.034, 68.828, 71.151, 70.21, 71.573,0]
etrgbt0 = [0,1.101, 0.6835, 1.3085, 0.673, 0.716, 0.6535, 1.117, 0.8185, 0.8665,0]

trgball = [0,67.499, 69.761, 68.917, 69.516, 69.264, 68.701, 71.252, 69.288, 70.463,0]

etrgball = [0,1.037, 0.701, 1.2645, 0.698, 0.7164, 0.695, 1.1255, 0.8305, 0.823,0]


pl.errorbar(x,trgb,yerr=etrgb,fmt='d',color='g',ms=14,markeredgewidth=2,label='$CSPI&II$',lw=2)

pl.errorbar(x-.15,trgb1,yerr=etrgb1,fmt='o',color='b',ms=12,markeredgewidth=2,label='$CSPI$',lw=2)
pl.errorbar(x+.15,trgb2,yerr=etrgb2,fmt='s',color='r',ms=10,markeredgewidth=2,label='$CSPII$',lw=2)
pl.errorbar(x+.20,wendy19,yerr=ewendy19,fmt='o',color='k',ms=10,markeredgewidth=2,label='$Freedman+19$',lw=2)



pl.legend(loc='lower right',numpoints=1,fontsize=10)
pl.grid()
pl.xlim(.5,9.5),pl.ylim(65,74)
pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)

pl.savefig('../../plots/h0bandsTRGB.pdf')


pl.figure(4)
pl.xticks(x,lam)

ceph1 = [0,77.286, 74.887, 76.479, 75.044, 73.946, 74.633, 73.185, 73.585, 74.126,0]
eceph1 = [0, 1.416, 0.961, 1.2475, 0.911, 0.8214, 0.8565, 0.994, 0.898, 0.9435 ,0]

ceph2 = [0,72.896, 72.246, 74.132, 72.201, 71.843, 71.947, 75.59, 74.352, 76.572,0]
eceph2 = [0,1.488, 0.853, 1.24, 0.823, 0.801, 0.714, 1.0745, 0.9904, 0.981 ,0]

burns18 = [0,73.98,72.74,72.64,74.32,71.85,72.98,72.25,72.47,73.84,0]
ebursn18 = [0, 3.10,1.6,1.57,2.45,1.48,1.54,2.35,1.74,1.78,0]


pl.errorbar(x,ceph,yerr=eceph,fmt='d',color='g',ms=14,markeredgewidth=2,label='$CSPI&II$',lw=2)
pl.errorbar(x-.15,ceph1,yerr=eceph1,fmt='o',color='b',ms=12,markeredgewidth=2,label='$CSPI$',lw=2)
pl.errorbar(x+.15,ceph2,yerr=eceph2,fmt='s',color='r',ms=10,markeredgewidth=2,label='$CSPII$',lw=2)
pl.errorbar(x+.20,burns18,yerr=ebursn18,fmt='o',color='k',ms=10,markeredgewidth=2,label='$Burns+18$',lw=2)

pl.legend(loc='upper right',numpoints=1,fontsize=10)
pl.grid()
pl.xlim(.5,9.5),pl.ylim(70,80)
pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)

pl.savefig('../../plots/h0bandsCeph.pdf')


pl.figure(5)
pl.xticks(x,lam)

sbf1 = [0,69.645, 71.806, 72.802, 71.487, 71.562, 71.887, 68.922, 68.315, 67.391,0]
esbf1 = [0, 1.68, 1.1875, 1.527, 1.297, 0.9875, 0.956, 1.414, 1.1915, 1.112 ,0]

sbf2 = [0,67.759, 71.984, 73.405, 71.628, 71.369, 70.967, 73.137, 68.898, 69.05,0]
esbf2 = [0,1.831, 0.9924, 1.8405, 0.9929, 0.905, 0.84, 1.5514, 1.233, 1.2865 ,0]


pl.errorbar(x,sbf,yerr=esbf,fmt='d',color='g',ms=14,markeredgewidth=2,label='$CSPI&II$',lw=2)
pl.errorbar(x-.15,sbf1,yerr=esbf1,fmt='o',color='b',ms=12,markeredgewidth=2,label='$CSPI$',lw=2)
pl.errorbar(x+.15,sbf2,yerr=esbf2,fmt='s',color='r',ms=10,markeredgewidth=2,label='$CSPII$',lw=2)


pl.legend(loc='upper right',numpoints=1,fontsize=10)
pl.grid()
pl.xlim(.5,9.5),pl.ylim(65,75)
pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)

pl.savefig('../../plots/h0bandsSBF.pdf')




pl.figure(6)
pl.xticks(x,lam)

cephc=[0,76.583, 72.555, 75.172, 72.446, 71.218, 71.981, 78.715, 75.762, 76.171,0]
ecephc=[0,1.404, 0.908, 1.4815, 0.8945, 0.881, 0.834, 1.5185, 1.0445, 1.0015,0]


trgbc=[0,71.69, 70.324, 69.007, 70.213, 69.452, 69.906, 74.0, 73.384, 74.08,0]
etrgbc=[0,1.4, 0.892, 1.785, 0.822, 0.8145, 0.743, 1.556, 1.0185, 0.9305,0]

pl.errorbar(x,cephc,yerr=ecephc,fmt='o',color='b',ms=14,markeredgewidth=2,label='Cepheid',lw=2)

pl.errorbar(x+.1,trgbc,yerr=etrgbc,fmt='d',color='g',ms=14,markeredgewidth=2,label='TRGB',lw=2)

pl.legend(loc='upper right',numpoints=1,fontsize=10)
pl.grid()
pl.xlim(.5,9.5),pl.ylim(65,80)
pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)

pl.savefig('../../plots/h0common.pdf')


pl.figure(7)
pl.xticks(x,lam)

allt = [0,67.89, 70.489, 70.75, 70.521, 69.259, 69.443, 72.699, 71.388, 71.558,0]
eallt=[0,1.2195, 0.745, 1.2505, 0.7475, 0.74249, 0.646, 1.098, 0.8434, 0.81,0]

noss =[0,67.328, 70.407, 71.903, 70.204, 69.077, 68.942, 72.229, 70.82, 70.908,0]
enoss=[0,1.077, 0.7284, 1.282, 0.704, 0.682, 0.6145, 1.1745, 0.768, 0.7905,0]

pl.errorbar(x,allt,yerr=eallt,fmt='o',color='b',ms=14,markeredgewidth=2,label='All Branch Types',lw=2)

pl.errorbar(x+.1,noss,yerr=enoss,fmt='d',color='g',ms=14,markeredgewidth=2,label='Excluding SS Type',lw=2)
pl.xlim(.5,9.5),pl.ylim(65,75)
pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)
pl.legend(loc='upper left',numpoints=1,fontsize=10)
pl.grid()
pl.savefig('../../plots/h0branchTRGB.pdf')



pl.figure(8)
pl.xticks(x,lam)

allc = [0,74.1, 72.861, 74.801, 72.767, 72.073, 72.409, 74.409, 73.622, 74.768,0]
eallc=[0,1.232, 0.7415, 1.141, 0.752, 0.688, 0.687, 0.999499, 0.821, 0.8205,0]

nossc =[0,74.151, 73.133, 77.518, 73.186, 72.247, 72.163, 75.833, 73.855, 75.068,0]
enossc=[0,1.2445, 0.786, 1.335, 0.735, 0.732499, 0.6645, 1.093, 0.8545, 0.866,0]

pl.errorbar(x,allc,yerr=eallc,fmt='o',color='b',ms=14,markeredgewidth=2,label='All Branch Types',lw=2)

pl.errorbar(x+.1,nossc,yerr=enossc,fmt='d',color='g',ms=14,markeredgewidth=2,label='Excluding SS Type',lw=2)
pl.xlim(.5,9.5),pl.ylim(70,80)
pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)
pl.legend(loc='upper right',numpoints=1,fontsize=10)
pl.grid()
pl.savefig('../../plots/h0branchCepheids.pdf')


## from updata2 files

pl.figure(9)
pl.xticks(x,lam)

trgbu = [0,68.834, 69.879, 69.57, 69.598, 68.713, 69.26, 70.776, 72.591, 71.066,0]
etrgbu=[0,1.0775, 0.7555, 1.1995, 0.714, 0.6105, 0.608, 1.097, 0.8115, 0.7335,0]

cephu = [0,75.576, 73.381, 75.484, 73.442, 72.407, 73.414, 73.434, 75.73, 74.578,0]
ecephu=[0,1.168, 0.7164, 1.0415, 0.713, 0.667, 0.668, 0.945, 0.827, 0.7955,0]

k21u = [0,66.235, 69.548, 68.874, 69.731, 69.373, 69.975, 64.337, 67.63, 65.101,0]
ek21u=[0,1.559, 0.960, 1.73349, 0.907, 0.915, 0.895, 1.808, 1.1535, 1.108,0]


j21u = [0,72.993, 77.022, 76.84, 77.274, 76.209, 76.641, 73.895, 75.629, 76.021,0]
ej21u=[0,1.6745, 1.1935, 1.7575, 1.067, 1.064, 1.004, 1.6, 1.496, 1.4175,0]


pl.errorbar(x,cephu,yerr=ecephu,fmt='o',color='b',ms=14,markeredgewidth=2,label='Cepheids',lw=2)

pl.errorbar(x-.15,trgbu,yerr=etrgbu,fmt='d',color='g',ms=14,markeredgewidth=2,label='TRGB',lw=2)



pl.errorbar(x+.15,k21u,yerr=ek21u,fmt='s',color='#377eb8',ms=12,markeredgewidth=2,label='$SBF \ (K21)$',lw=2)
pl.errorbar(x+.20,j21u,yerr=ej21u,fmt='s',color='r',ms=10,markeredgewidth=2,label='$SBF \ (J21)$',lw=2)


pl.xlim(.5,9.5),pl.ylim(62,80)


pl.xlabel(r'$Bands$',fontsize=14), pl.ylabel(r'$H_0 \ (km \ s^{-1}\ Mpc^{-1})$',fontsize=14)
pl.legend(loc='lower center',numpoints=1,fontsize=10)
pl.grid()
pl.savefig('../../plots/h0update2.pdf')
