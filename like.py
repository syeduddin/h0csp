# Testing this out on the new master table:
from astropy.io import ascii
import numpy as np

tab = ascii.read('Bmax+Mass.dat')


Ho_dists = tab['dist'] < 0
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
c1 = edist*0  # not sure what these are
c2 = edist*0


def distmod(Ho, zhel,zcmb):
   # Just a quick one for testing
   return 5*np.log10(3e5*zcmb/Ho) + 25


# Liklihood function
def like(par):
    p,p1,p2,rv,alpha,sig,h0 = par
    if  -25.0<p<-15.0  and -10.0<p1<10.0 and -10.0<p2<10.0 and -10.0<rv<10.0 \
      and -10<alpha<10 and 0<sig<10 and 0< h0 < 100.0: # priors
        
        #
        #cosmo=FlatLambdaCDM(H0=h0, Om0=.30)

        st1 = p1*(st - 1)
        st2 = p2*((st - 1)**2)
        red = rv*(bv)
        mu_obs = mmax - p - st1 - st2 - red - alpha*(m_csp > 10.27),
        mu_model = np.where(Ho_dists, distmod(h0,zhel,zcmb), dist)
        mu_stat = np.where(Ho_dists,
                  ((est**2) +(emmax**2) +(ebv**2)+(c1**2)+(c2**2)+(sig**2)),
                  ((est**2) +(emmax**2) +(ebv**2)+(c1**2)+(c2**2)+(sig**2)))
      
        mu_stat=np.sqrt(mu_stat)
        dmu=mu_obs-mu_model
        chi =np.sum((dmu)**2/mu_stat**2)
        return -0.5*chi - 0.5*np.sum(np.log(mu_stat))
    else:
        return -np.inf

if __name__ == "__main__":
    # TEst it out:

    p = (-19, -0.6, 0.0, 2.0, 0.05, 0.08, 72.0)
    print(like(p))
