import matplotlib.pyplot as pl
from scipy import stats
import numpy as np

x=np.linspace(50,90,1000)
b=stats.norm.pdf(x, 71.76, 1.19) # B-band
h=stats.norm.pdf(x, 73.22, 1.28) # H-band

reiss = stats.norm.pdf(x, 73.04, 1.04)
wendy = stats.norm.pdf(x, 69.80, 1.6)
peter = stats.norm.pdf(x, 74.60, 2.7)
nandita = stats.norm.pdf(x, 70.50, 4.12)
cmb = stats.norm.pdf(x, 67.40, 0.5)


pl.figure(figsize=(20,10))
pl.plot(x,b,lw=5,color='b')
pl.fill_between(x,b,color='b',alpha=.3)
pl.plot(x,h,lw=5,color='#e41a1c')
pl.fill_between(x,h,color='#e41a1c',alpha=.3)

pl.plot(x,reiss,lw=3,color='#4daf4a')
pl.fill_between(x,reiss,color='#4daf4a',alpha=.3)

pl.plot(x,wendy,lw=3,color='#f781bf')
pl.fill_between(x,wendy,color='#f781bf',alpha=.3)

pl.plot(x,peter,lw=3,color='y')
pl.fill_between(x,peter,color='y',alpha=.3)


pl.plot(x,nandita,lw=3,color='m')
pl.fill_between(x,nandita,color='m',alpha=.3)

pl.plot(x,cmb,lw=3,color='k')
pl.fill_between(x,cmb,color='k',alpha=.3)





pl.ylim(0,0.8),pl.xlim(55,85)
pl.xlabel(r'$H_0 \ (km \ s^{-1} \ Mpc^{-1})$',fontsize=20), pl.ylabel(r' $Normalized \ Probability$',fontsize=20)
pl.grid(True,alpha=.5)
pl.legend(['$This \ Work \ (B)$','$This \ Work \ (H)$','$Riess \ et \ al. \ (2022)$', '$Freedman \ et \ al. \ (2021)$','$Garnavich \ et \ al. \ (2022)$', '$Khetan \ et \ al. \ (2021)$', '$Plank Collaboraion \ et \ al. \ (2018)$'],fontsize=18,loc ='upper right')

l1 = pl.legend(['$This \ Work \ (B)$','$This \ Work \ (H)$','$Riess \ et \ al. \ (2022)$', '$Freedman \ et \ al. \ (2021)$','$Garnavich \ et \ al. \ (2022)$', '$Khetan \ et \ al. \ (2021)$', '$Plank \ Collaboration \ et \ al. \ (2018)$'],fontsize=18,loc ='upper right')
pl.gca().add_artist(l1)



pl.legend(['$71.76 \pm 0.58 \ (stat)  \pm 1.19 \ (sys)$','$73.22 \pm 0.68 \ (stat) \pm 1.28 \ (sys)$', '$73.04 \pm 1.04 \ (total)$','$69.80 \pm 0.60 \ (stat) \pm 1.60 \ (sys)$','$74.60 \pm 0.90 \ (stat) \pm 2.70 \ (sys)$', '$70.50 \pm 2.37 \ (stat) \pm 3.38 \ (sys)$','$67.40 \pm 0.50 \ (total)$' ],fontsize=18,loc ='upper left')



pl.xticks(fontsize=18),pl.yticks(fontsize=18)

pl.savefig('../../plots/pdfH0.pdf')

