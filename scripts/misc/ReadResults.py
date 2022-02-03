from astropy.io import ascii
import math
filter = ['u','B','g','V','r','i','Y','J','H']
#filter ='B'

scatter=[]
escatter=[]

for i in range(len(filter)):
    
    
    #result = ascii.read('../../results/'+filter[i]+'_sbf_results_cut.txt')
    result = ascii.read('../../results/SBF_result_'+filter[i]+'.txt')
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

    scatter.append(sig)
    escatter.append(esig)
    print ('& $',filter[i],'$', '&', '%0.2f'%h0,'(%0.2f)'%eh0,'&','%0.2f'%sig, '(%s)'%esig[1][0:2],'&', '%d'%vel,'(%d)'%evel, '&','%0.2f'%p0,'(%s)'%ep0[1][0:2],'&','%0.2f'%p1,'(%s)'%ep1[1][0:2],'&','%0.2f'%p2,'(%s)'%ep2[1][0:2],'&','%0.2f'%alpha,'(%s)'%ealpha[1][0:2],'&','%0.2f'%beta,'(%s)'%ebeta[1][0:2])

#print (escatter)
