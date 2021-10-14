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
#filter = raw_input("Please enter a filter name:\n")

file = sys.argv[1]

ncpu = cpu_count()
print("{0} CPUs".format(ncpu))

os.environ["OMP_NUM_THREADS"] = "1"
c = 300000.
q=-0.53


tab = ascii.read('../data/working/'+file)



w = np.where((tab['sn']!='CSP14abk') &  (tab['sn']!='PTF13dyt') &  (tab['sn']!='PTF13dym') &  (tab['sn']!='PS1-13eao'))



st = tab['st'][w]
est = tab['est'][w]
zhel = tab['zhel'][w]
zcmb = tab['zcmb'][w]
mmax = tab['Mmax'][w]
emmax = tab['eMmax'][w]
bv = tab['BV'][w]
ebv = tab['eBV'][w]
m_csp = tab['m'][w]
em = (tab['mu'][w]+tab['ml'][w])/2.
dist = tab['dist'][w]
edist = tab['edist'][w]
c_ms = tab['covMs'][w]
c_mbv = tab['covBV_M'][w]
sn = tab['sn'][w]

Ho_dists = (dist < 0) 
print (file, len(st), len(st[Ho_dists]))

ss= np.where(dist>1)
print (len(st[ss]))

result = ascii.read('../results/'+file[:-4]+'_result_cal.txt')
print (result)
p=result['p0'][0]
p1=result['p1'][0]
p2=result['p2'][0]
alpha=result['alpha'][0]
rv=result['beta'][0]
sig=result['sig_int'][0]
vel=result['vel'][0]
    

#initial guess
h0lim= 71.0,71.1

# Eqn 9 of Bruns 2018
def distmod(h,z1,z2):
    t1 = (1+z1)/(1+z2)
    t2 = (c*z2)/h
    t3 = 1 + (((1-q)/2)*z2)
    return (5*np.log10(t1*t2*t3)) +25

# Liklihood function
def like(par):
    h0 = par
    if    0< h0 < 1000.0: # priors
        

        
        mu_obs = mmax - p - p1*(st - 1.) -  p2*(st - 1.)**2 - rv*bv - alpha*(m_csp-np.median(m_csp))
        #mu_obs = mmax - p - p1*(st - 1.) -  p2*(st - 1.)**2 - rv*bv - (alpha*m_csp) # slope 


        mu_model = np.where(Ho_dists, distmod(h0,zhel,zcmb), dist)
        fac= (p1+(2*p2*st))
        velterm = (2.17*vel)**2/(c*zcmb)**2
        err = (fac*est)**2 +emmax**2 +(rv*ebv)**2+2*fac*c_ms+rv*c_mbv+sig**2+(0.00000723*vel/zcmb)**2 +(alpha*em)**2
        err1 = ((fac*est)**2) +(emmax**2) +((rv*ebv)**2)+(2*fac*c_ms)+(rv*c_mbv)+(edist**2)+(alpha*em)**2
    
        mu_stat = np.where(Ho_dists,err,err1)
        
      
        mu_stat=np.sqrt(mu_stat)
        dmu=mu_obs-mu_model
        
        chi =np.sum((dmu)**2/mu_stat**2)
        return -0.5*chi - (0.5*np.sum(np.log(2*np.pi*(mu_stat)**2))) 
    else:
        return -np.inf
# EMCEE
ndim, nwalkers = 1, 10
ssize=1000
burnin = 500



h00 = np.random.rand(nwalkers) * (h0lim[1] - h0lim[0]) + h0lim[0]

#p0 = zip(*[p00,p10,p20,rv0,alpha0,sig0,vel0,h00])
p0 = np.array([h00]).T




sampler = emcee.EnsembleSampler(nwalkers, ndim, like)
print ("running mcmc on "+file)
start = time.time()
sampler.run_mcmc(p0,ssize,progress=True)
samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
end = time.time()
serial_time = end - start


# Chains
#fig, axes = pl.figure(figsize=(10, 7))
samples = sampler.get_chain()
labels = [ r"$H_0$"]
#ax = axes
#pl.plot(samples[:,:,0], "k", alpha=0.3)
#pl.set_xlim(0, len(samples))
#pl.set_ylabel(labels)
#pl.yaxis.set_label_coords(-0.1, 0.5)

#axes[-1].set_xlabel("step number")

#pl.savefig("../plots/steps_H0_"+str(nwalkers)+"_"+str(ssize)+".pdf")

sample = sampler.chain[:, burnin:, :].reshape((-1, ndim))


 # Printing results
#H0_mcmc = (lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             #(*np.percentile(sample, [16, 50, 84],
                              #                  axis=0)))

mcmc = np.percentile(sample, [16, 50, 84])
q = np.mean(np.diff(mcmc))
print ('H0=','%0.2f'%mcmc[1],'+\-','%0.2f'%q)
#print (H0_mcmc)
#print("""MCMC result:
    
 #   H0 = {0[0]} +{0[1]} -{0[2]}

   
#""".format(H0_mcmc))





print ("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))




pl.hist(sample[:,0], 100, color="k", histtype="step")
pl.savefig("../plots/H0_"+file[:-4]+".pdf")


print("Serial took {0:.1f} minutes".format(serial_time/60.))

#os.system('say "your program has finished."')

    












