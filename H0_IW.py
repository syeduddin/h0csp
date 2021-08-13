import sys
import numpy as np
import emcee
import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import triangle
import random,os
from multiprocessing import Pool
from multiprocessing import cpu_count
import time
from astropy.io import ascii
from scipy.stats import multivariate_normal
from scipy.stats import invgamma
#filter = raw_input("Please enter a filter name:\n")

filter = sys.argv[1]
print  'Working with filter:'+filter

ncpu = cpu_count()
print("{0} CPUs".format(ncpu))

os.environ["OMP_NUM_THREADS"] = "1"
c = 300000.
q=-0.53


tab = ascii.read('data/'+filter+'_all.dat')
st = tab['st']
est = tab['est']
zhel = tab['zhel']
zcmb = tab['zcmb']
mmax = tab['Mmax']
emmax = tab['eMmax']
bv = tab['BV']
ebv = tab['eBV']
m_csp = tab['m']
dist = tab['dist']
edist = tab['edist']
c_ms = tab['covMs']
c_mbv = tab['covBV_M']
methods = tab['caltype']

w1 = np.where(dist<0.) # SNe Ia
w2 = np.where(dist>0.) # Calibrators 
#initial guess
plim=-19.3, -19.2
p1lim =-1.2,-1.1
p2lim=-.055,-0.05
rvlim =2.7,2.71
alphalim=-0.11,-0.1
siglim=0.1,.12
vellim =300.,310
h0lim= 71.0,71.1
sd_callim = .1,.11
alim = .1,.11
blim = .1,.11

# Eqn 9 of Bruns 2018
def distmod(h,z1,z2):
    t1 = (1+z1)/(1+z2)
    t2 = (c*z2)/h
    t3 = 1 + (((1-q)/2)*z2)
    return (5*np.log10(t1*t2*t3)) +25

# Covariance matrix for calibrtors
def covmat(p1,p2,rv,sd_cal):
    # Error in calibrators
    fac= (p1+(2*p2*st[w2]))
    error = ((fac*est[w2])**2) +(emmax[w2]**2) +((rv*ebv[w2])**2)+(2*fac*c_ms[w2])+(rv*c_mbv[w2])+(edist[w2]**2)
    
    num_vars = len(dist[w2])
    cov = np.zeros((num_vars, num_vars), dtype=float)          
    np.fill_diagonal(cov,error)

    method_same = (methods[w2][:,np.newaxis] == methods[w2][np.newaxis,:])
    cov = cov + method_same*sd_cal**2 
    return cov

#print covmat(-1,-1,2,1)
#sys.exit()


# Liklihood function
def like(par):
    p,p1,p2,rv,alpha,sig,vel,h0, sd_cal,a,b = par
    if  -25.0<p<14.0  and -10.0<p1<10.0 and -10.0<p2<10.0 and 0.<rv<10.0 and -1.<alpha<1. and 0.<sig<1. and 0.<vel<1000. and  0< h0 < 1000.0 and 0.<sd_cal<10. and 0<a<10 and 0<b<10:  # priors

        # SNe Ia
        mu_obs = mmax[w1] - p - p1*(st[w1] - 1.) -  p2*(st[w1] - 1.)**2 - rv*bv[w1] - alpha*(m_csp[w1]-np.median(m_csp[w1])) 
        mu_model =  distmod(h0,zhel[w1],zcmb[w1])
        fac= (p1+(2*p2*st[w1]))
        mu_stat = (fac*est[w1])**2 +emmax[w1]**2 +(rv*ebv[w1])**2+2*fac*c_ms[w1]+rv*c_mbv[w1]+sig**2+(0.00000723*vel/zcmb[w1])**2
        mu_stat = np.sqrt(mu_stat)
        dmu=mu_obs-mu_model
        chi =np.sum((dmu)**2/mu_stat**2)
        
        lnlike1 = -0.5*chi - (0.5*np.sum(np.log(2*np.pi*(mu_stat)**2)))

        ## Calibrator
        mu_obscal = mmax[w2] - p - p1*(st[w2] - 1.) -  p2*(st[w2] - 1.)**2 - rv*bv[w2] - alpha*(m_csp[w2]-np.median(m_csp[w2])) 

        # New covmat+ inverse-gamma distribution
       
        lnlike2= np.sum(np.log(multivariate_normal.pdf(dist[w2], mean=mu_obscal, cov=covmat(p1,p2,rv,sd_cal))))\
+invgamma.logpdf(sd_cal**2, a=a, scale=b)
        
        return lnlike1+ lnlike2
    else:
        return -np.inf

# EMCEE
ndim, nwalkers = 11, 110
ssize=1000
burnin = 500


p00 = np.random.rand(nwalkers) * (plim[1] - plim[0]) + plim[0]
p10 = np.random.rand(nwalkers) * (p1lim[1] - p1lim[0]) + p1lim[0]
p20 = np.random.rand(nwalkers) * (p2lim[1] - p2lim[0]) + p2lim[0]
rv0 = np.random.rand(nwalkers) * (rvlim[1] - rvlim[0]) + rvlim[0]
alpha0 = np.random.rand(nwalkers) * (alphalim[1] - alphalim[0]) + alphalim[0]
sig0 = np.random.rand(nwalkers) * (siglim[1] - siglim[0]) + siglim[0]
vel0 = np.random.rand(nwalkers) * (vellim[1] - vellim[0]) + vellim[0]
h00 = np.random.rand(nwalkers) * (h0lim[1] - h0lim[0]) + h0lim[0]
sd_cal0 = np.random.rand(nwalkers) * (sd_callim[1] - sd_callim[0]) + sd_callim[0]
a0 = np.random.rand(nwalkers) * (alim[1] - alim[0]) + alim[0]
b0 = np.random.rand(nwalkers) * (blim[1] - blim[0]) + blim[0]

p0 = zip(*[p00,p10,p20,rv0,alpha0,sig0,vel0,h00,sd_cal0,a0,b0])


sampler = emcee.EnsembleSampler(nwalkers, ndim, like,pool=Pool())
print "running mcmc.."
start = time.time()
sampler.run_mcmc(p0,ssize,progress=True)
samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
end = time.time()
serial_time = end - start


# Chains
fig, axes = pl.subplots(11, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["$M_B$","P1", "P2", r"$\beta$",r"$\alpha$", r"$\sigma_{int}$",r"$V_{pec}$", r"$H_0$",r"$\sigma_{cal}$","$a$", "$b$"]

for j in range(ndim):
    ax = axes[j]
    ax.plot(samples[:, :, j], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[j])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")

fig.savefig("plots/steps_sd"+filter+"_"+str(nwalkers)+"_"+str(ssize)+".pdf")

samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))


 # Printing results
p0_mcmc,p1_mcmc,p2_mcmc,rv_mcmc,alpha_mcmc,sig_mcmc,vel_mcmc, H0_mcmc, sd_cal_mcmc,a_mcmc,b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))


print("""MCMC result:
    P0 = {0[0]} +{0[1]} -{0[2]} 
    P1 = {1[0]} +{1[1]} -{1[2]} 
    P2 = {2[0]} +{2[1]} -{2[2]} 
    Beta = {3[0]} +{3[1]} -{3[2]}
    Alpha = {4[0]} +{4[1]} -{4[2]}
    Sigma = {5[0]} +{5[1]} -{5[2]}
    Vpec = {6[0]} +{6[1]} -{6[2]}
    H0 = {7[0]} +{7[1]} -{7[2]}
    sd_cal = {8[0]} +{8[1]} -{8[2]}
a = {9[0]} +{9[1]} -{9[2]}
b = {10[0]} +{10[1]} -{10[2]}


   
""".format(p0_mcmc, p1_mcmc, p2_mcmc,rv_mcmc,alpha_mcmc,sig_mcmc,vel_mcmc, H0_mcmc, sd_cal_mcmc,a_mcmc,b_mcmc))





# Triangle plot
figure = triangle.corner(samples,labels=["$P0$","$P1$", "$P2$", r"$\beta$",r"$\alpha$", r"$\sigma_{int}$","$V_{pec}$", r"$H_0$",r"$\sigma_{cal}$","$a$", "$b$" ],quantiles=[0.16, 0.5, 0.84],truths=[p0_mcmc[0],p1_mcmc[0],p2_mcmc[0],rv_mcmc[0],alpha_mcmc[0],sig_mcmc[0],vel_mcmc[0],H0_mcmc[0], sd_cal_mcmc[0],a_mcmc[0],b_mcmc[0]],show_titles=True)

#figure.savefig("plots/mcmcH0_"+filter+"_"+str(nwalkers)+"_"+str(ssize)+".pdf")
figure.savefig("plots/mcmcH0_sd"+filter+"_"+str(nwalkers)+"_"+str(ssize)+".pdf")



#pl.hist(samples[:, 4], 500, color="k", histtype="step")
#pl.savefig("H0.pdf")


print "Mean acceptance fraction:", np.mean(sampler.acceptance_fraction)
print("Serial took {0:.1f} minutes".format(serial_time/60.))

os.system('say "your program has finished."')

    













#num_vars = len(dist[w2])
 #   cov = np.zeros((num_vars, num_vars), dtype=float)
  #  for i in range (0,num_vars):
   #     for j in range (0,num_vars):
    #        if i==j:
     #           cov[i][j]= error[i] + sd_cal**2
      #      else:
       #         cov[i][j]=sd_cal**2
            
