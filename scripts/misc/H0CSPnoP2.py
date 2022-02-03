import sys
import numpy as np
import emcee
import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import random,os
from multiprocessing import Pool
from multiprocessing import cpu_count
import time
from astropy.io import ascii
import corner
#from astropy.table import vstack
#filter = raw_input("Please enter a filter name:\n")

file = sys.argv[1]

ncpu = cpu_count()
print("{0} CPUs".format(ncpu))

os.environ["OMP_NUM_THREADS"] = "1"
c = 300000.
q=-0.53


tab = ascii.read('../../data/working/'+file)

#st = tab['st']
#est = tab['est']
#zhel = tab['zhel']
#zcmb = tab['zcmb']
#mmax = tab['Mmax']
#emmax = tab['eMmax']
#bv = tab['BV']
#ebv = tab['eBV']
#m_csp = tab['m']
#eml = (tab['m']-tab['ml'])
#emu = (tab['mu']-tab['m'])
#em = (emu+eml)/2.
#dist = tab['dist']
#edist = tab['edist']
#c_ms = tab['covMs']
#c_mbv = tab['covBV_M']


st = tab['stretch']
est = tab['err_stretch']
zhel = tab['redshift']
zcmb = tab['redshift']
mmax = tab['Bmax']
emmax = tab['eBmax']
bv = tab['color']
ebv = tab['e_color']
m_csp = tab['log_Mass']
em = tab['e_logMass']
dist = tab['mu_SBF']
edist = tab['e_mu']
c_ms = tab['color']*0.
c_mbv =tab['color']*0.

#w0 = ((tab['sn']!='CSP14abk') &  (tab['sn']!='PTF13dyt') &  (tab['sn']!='PTF13dym') &  (tab['sn']!='PS1-13eao')& (tab['zcmb']>0.01) & (tab['dist']<0)& (tab['t0']<5) & (tab['st']> 0.5) & (tab['BV']<0.5))

#w0 = ((tab['sn']!='CSP14abk') &  (tab['sn']!='PTF13dyt') &  (tab['sn']!='PTF13dym') &  (tab['sn']!='PS1-13eao') & (tab['dist']<0) & (tab['zcmb']>0.0) & (tab['zcmb']<0.75))
w0 = ((dist<0))# & (tab['redshift']>0.02) & (tab['redshift']<0.075))
f1 =open('../../results/'+file[:-4]+'_results_nop2.txt','w')

w1 = (dist>0)

print (file, len(st[w0]),len(st[w1]))

#mmax[w0] = (mmax[w0]*.99)+.22
#st[w0] = (st[w0]*1.02)-0.03
#bv[w0] = (bv[w0]*1.0)+0.04
#mmax[w1] = (mmax[w1]*1.)+.06
#st[w1] = (st[w1]*1.27)-0.21
#bv[w1] = (bv[w1]*0.98)+0.01


#initial guess
plim=-19.3, -19.2
p1lim =-1.2,-1.1
p2lim=-.055,-0.05
rvlim =2.7,2.71
alphalim=-0.11,-0.1
siglim=0.1,.12
vellim =300.,310
h0lim= 71.0,71.1

# Eqn 9 of Bruns 2018
def distmod(h,z1,z2):
    t1 = (1+z1)/(1+z2)
    t2 = (c*z2)/h
    t3 = 1 + (((1-q)/2)*z2)
    return (5*np.log10(t1*t2*t3)) +25

# Liklihood function
def like(par):
    p,p1,rv,alpha,sig1,sig2,h0 = par

    mu_obs=[]
    mu_model=[]
    mu_stat=[]
    
    if  -25.0<p<14.0  and -10.0<p1<10.0 and 0.<rv<10.0 and -1.<alpha<1. and 0.<sig1<1. and 0.<sig2<1. and  0< h0 < 1000.0: # priors
        
        
        mu_sn = mmax[w0] - p - p1*(st[w0] - 1.)  - rv*bv[w0] - alpha*(m_csp[w0]-np.median(m_csp[w0]))

        mu_cal = mmax[w1] - p -p1*(st[w1] - 1.) - rv*bv[w1] - alpha*(m_csp[w1]-np.median(m_csp[w1]))

        for mu_sn in mu_sn:
            mu_obs.append(mu_sn)
        for mu_cal in mu_cal:
            mu_obs.append(mu_cal)
        
        mu_model_sn = distmod(h0,zhel[w0],zcmb[w0])
        mu_model_cal = dist[w1]
        for mu_model_sn in mu_model_sn:
            mu_model.append(mu_model_sn)
        for mu_model_cal in mu_model_cal:
            mu_model.append(mu_model_cal)
        
        
        err = (p1*est[w0])**2 +emmax[w0]**2 +(rv*ebv[w0])**2+sig1**2-2*rv*emmax[w0]**2+(alpha*em[w0])**2
        err1 = (p1*est[w1])**2 +emmax[w1]**2 +(rv*ebv[w1])**2+sig2**2-2*rv*emmax[w1]**2 +edist[w1]**2+(alpha*em[w1])**2
            
        for err in err:
            mu_stat.append(err)
        for err1 in err1:
            mu_stat.append(err1)
        mu_stat=np.sqrt(mu_stat)

        mu_obs = np.array(mu_obs)
        mu_model = np.array(mu_model)
        mu_stat = np.array(mu_stat)
        dmu=mu_obs-mu_model
        
        chi =np.sum((dmu)**2/mu_stat**2)
        return -0.5*chi - (0.5*np.sum(np.log(2*np.pi*(mu_stat)**2))) 
    else:
        return -np.inf
# EMCEE
ndim, nwalkers = 7, 80
ssize=1000
burnin = 200


p00 = np.random.rand(nwalkers) * (plim[1] - plim[0]) + plim[0]
p10 = np.random.rand(nwalkers) * (p1lim[1] - p1lim[0]) + p1lim[0]
p20 = np.random.rand(nwalkers) * (p2lim[1] - p2lim[0]) + p2lim[0]
rv0 = np.random.rand(nwalkers) * (rvlim[1] - rvlim[0]) + rvlim[0]
alpha0 = np.random.rand(nwalkers) * (alphalim[1] - alphalim[0]) + alphalim[0]
sig01 = np.random.rand(nwalkers) * (siglim[1] - siglim[0]) + siglim[0]
sig02 = np.random.rand(nwalkers) * (siglim[1] - siglim[0]) + siglim[0]
h00 = np.random.rand(nwalkers) * (h0lim[1] - h0lim[0]) + h0lim[0]


p0 = np.array([p00,p10,rv0,alpha0,sig01,sig02,h00]).T



sampler = emcee.EnsembleSampler(nwalkers, ndim, like)
print ("running mcmc..")
start = time.time()
sampler.run_mcmc(p0,ssize,progress=True)
samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
end = time.time()
serial_time = end - start


# Chains
fig, axes = pl.subplots(7, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["$M_B$","P1", "$R_V$",r"$\alpha$", r"$\sigma_{int1}$",r"$\sigma_{int2}$", r"$H_0$"]
for j in range(ndim):
    ax = axes[j]
    ax.plot(samples[:, :, j], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[j])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")

fig.savefig("../../plots/steps_sbf_"+str(nwalkers)+"_"+str(ssize)+".pdf")

samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))


 # Printing results
p0_mcmc,p1_mcmc,rv_mcmc,alpha_mcmc,sig1_mcmc,sig2_mcmc,H0_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))


print("""MCMC result:
    P0 = {0[0]} +{0[1]} -{0[2]} 
    P1 = {1[0]} +{1[1]} -{1[2]} 
    RV = {2[0]} +{2[1]} -{2[2]} 
    Alpha = {3[0]} +{3[1]} -{3[2]}
    Sigma1 = {4[0]} +{4[1]} -{4[2]}
    Sigma2 = {5[0]} +{5[1]} -{5[2]}
    H0 = {6[0]} +{6[1]} -{6[2]}

   
""".format(p0_mcmc, p1_mcmc,rv_mcmc,alpha_mcmc,sig1_mcmc,sig2_mcmc, H0_mcmc))





# Triangle plot
figure = corner.corner(samples,labels=["$P0$","$P1$", "$R_V$",r"$\alpha$", r"$\sigma_{int1}$",r"$\sigma_{int2}$", r"$H_0$"],quantiles=[0.16, 0.5, 0.84],truths=[p0_mcmc[0],p1_mcmc[0],rv_mcmc[0],alpha_mcmc[0],sig1_mcmc[0],sig2_mcmc[0],H0_mcmc[0]],show_titles=True)

figure.savefig("../../plots/mcmcH0_sbf_nop2"+str(nwalkers)+"_"+str(ssize)+".pdf")



#pl.hist(samples[:, 4], 500, color="k", histtype="step")
#pl.savefig("H0.pdf")


print ("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
print("Serial took {0:.1f} minutes".format(serial_time/60.))


    









