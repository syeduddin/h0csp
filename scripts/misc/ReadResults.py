from astropy.io import ascii
import math, sys
import numpy as np

filter = ['u','B','g','V','r','i','Y','J','H']
#filter =['B','H']

h=[]
eh=[]
cal = sys.argv[1]
for i in range(len(filter)):
    
    
    result = ascii.read('../../results/'+filter[i]+'_'+cal+'_update2_result.txt')
    
    #print (result)
    p0=result['p0'][0]
    ep0 = (result['p0'][1]+result['p0'][2])/2.
    ep0 = str(ep0).split('.')
    p1=result['p1'][0]
    ep1 = (result['p1'][1]+result['p1'][2])/2.
    ep1 = str(ep1).split('.')
    p2=result['p2'][0]
    ep2 = (result['p2'][1]+result['p2'][2])/2.
    ep2 = str(ep2).split('.')
    alpha=result['alpha'][0]
    ealpha= (result['alpha'][1]+result['alpha'][2])/2.
    ealpha = str(ealpha).split('.')
    beta=result['beta'][0]
    ebeta = (result['beta'][1]+result['beta'][2])/2.
    ebeta = str(ebeta).split('.')
    sig=result['sig_int'][0]
    esig = (result['sig_int'][1]+result['sig_int'][2])/2.
    esig = str(esig).split('.')
    vel=result['vel'][0]
    evel = (result['vel'][1]+result['vel'][2])/2.
    h0=result['H0'][0]
    eh0 = (result['H0'][1]+result['H0'][2])/2.

    h.append(h0)
    eh.append(eh0)
    #print ('& $',filter[i],'$', '&', '%0.2f'%h0,'(%0.2f)'%eh0,'&','%0.2f'%sig, '(%s)'%esig[1][0:2],'&', '%d'%vel,'(%d)'%evel, '&','%0.3f'%p0,'(%s)'%ep0[1][0:3],'&','%0.2f'%p1,'(%s)'%ep1[1][0:2],'&','%0.2f'%p2,'(%s)'%ep2[1][0:2],'&','%0.2f'%alpha,'(%s)'%ealpha[1][0:2],'&','%0.2f'%beta,'(%s)'%ebeta[1][0:2],'\\\\')

    #print (filter[i], '%0.2f'%h0, '(%0.2f)'%eh0)
print ((h))
print (eh)
#print (np.std([74.75,73.07,75.27,73.04,72.30,72.69,74.86,74.25,75.25]))




# All cuts
cephB =[72.39,73.13,73.43,73.07,73.21] # B
#print (np.std(cephB))
cephH =[75.27,75.06,74.84,75.69,74.28] # B
#print (np.std(cephH))

trgbB = [69.73,69.55,69.92,69.52,69.67]
trgbH = [71.11,71.00,70.84,71.42,70.42]
#print (np.std(trgbB))
#print (np.std(trgbH))

sbfB = [72.18,72.14,71.66,71.97,71.70 ]
sbfH = [66.87,67.09,67.21,66.99,66.64 ]
#print (np.std(sbfB))
#print (np.std(sbfH))

allB = [71.24,71.43,71.49,71.03,71.20]
allH = [70.61,70.98,71.85,70.12,70.27]

#print (np.std(allB))
#print (np.std(allH))

methodB = [73.07,69.76,72.62]
#print (np.std(methodB))
methodH = [75.25,71.83,69.13]
#print (np.std(methodH))

