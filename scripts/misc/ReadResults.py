from astropy.io import ascii

filter = ['u','B','g','V','r','i','Y','J','H']
for i in range(len(filter)):
    
    
    result = ascii.read('../results/'+filter[i]+'_trgb_result_cut.txt')
    #result = ascii.read('../results/TRGB_result_'+filter[i]+'.txt')
    p0=result['p0'][0]
    p1=result['p1'][0]
    p2=result['p2'][0]
    alpha=result['alpha'][0]
    beta=result['beta'][0]
    sig=result['sig_int'][0]
    vel=result['vel'][0]
    h0=result['H0'][0]
    h0h=result['H0'][1]
    h0l=result['H0'][2]
    
    print ('$',filter[i],'$', '&', '$','%0.2f'%h0,'\pm','%0.2f'%h0h,'$&','%0.3f'%sig)
    #print ('%0.3f'%sig)

