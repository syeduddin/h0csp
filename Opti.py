import sys
import numpy as np
from cosmolopy import *
import astropy.io.fits as pyfits
import matplotlib.pylab as pl
import triangle
import random,os
from astropy.cosmology import FlatLambdaCDM
import pandas as pd
from scipy import optimize

filter = raw_input("Please enter a filter name:\n")

print  'Working with filter:'+filter

# Data
## CSP LC
file ='CSPI+II_max/'+filter+'_max.dat'
print file
data = pd.read_csv(file,delim_whitespace=True)
df=pd.DataFrame(data)
mmax = df['Mmax']
emmax =df['eMmax']
st = df['st']
est = df['est']
zcmb= df['zcmb']
zhel = df['zhel']
distype = df['dist']
bv = df['BV']
ebv =df['eBV']
sn = df['name']
c1 = df['covMs']
c2 = df['covBV_M']
## CSP host mass
f1  = pd.read_csv('CSPHostMass.csv',delimiter=',')
d1=pd.DataFrame(f1)
m_csp = d1['m']
## Calibrator mass and distance
f2  = pd.read_csv('calibrators.csv',delimiter=',')
d2=pd.DataFrame(f2)
m_cal = d2['m']


m_sp1 =np.median(m_csp.values)
m_sp2 =np.median(m_cal.values)




# Likelihood for distance 
def like(par):
    p,p1,p2,rv,alpha,h0 = par
    
    cosmo=FlatLambdaCDM(H0=h0, Om0=.30)

    mu_obs=[]
    mu_model=[]
    mu_stat=[]
    for i in range(len(sn)):
        if df['dist'][i]=='Ho':
            if sn[i] not  in d1['sn'].values:w = np.where(sn[i]==d1['sn'])
            if sn[i] not in d1['sn1'].values:w = np.where(sn[i]==d1['sn1'])
            st1= p1*(st[i]-1)
            st2 =p2*((st[i]-1)**2)
            red = rv*(bv[i])
            mass = alpha*(m_csp.values[w]-m_sp1)
            mu_obs.append(mmax[i]-p-st1-st2-red-mass)
            distance = cosmo.distmod(zhel[i])
            mu_model.append(distance.value)
               
            mu_stat.append((est[i]**2) +(emmax[i]**2) +(ebv[i]**2)+(c1[i]**2)+(c2[i]**2))
            
        
        if df['dist'][i]<>'Ho':
            if sn[i] in d2['sn'].values:
                w1 = np.where(sn[i]==d2['sn'])
                st1= p1*(st[i]-1)
                st2 =p2*((st[i]-1)**2)
                red = rv*(bv[i])
                #mass = alpha*(m_csp.values[w]-m_sp2)
                #mu_obs.append(mmax[i]-p-st1-st2-red-mass)
                #mu_model.append(d2['dist'].values[w1])
                #mu_stat.append((est[i]**2) +(emmax[i]**2) +(ebv[i]**2)+(c1[i]**2)+(c2[i]**2))
        print len(mu_obs)
        mu_stat=np.sqrt(mu_stat)
        mu_stat=np.array(mu_stat)
        mu_obs=np.array(mu_obs)
        mu_model=np.array(mu_model)
        dmu=mu_obs-mu_model
        chi =np.sum((dmu)**2/mu_stat**2)
        return chi
    else:
        return -np.inf

result =optimize.fmin(like,[-19.0,-1.0,1.0,3.0,-.1,70.0])

