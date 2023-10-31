import os

cal = ['ceph','trgb','sbf','sbfj21','sbfcombined']
#cal = ['trgb']

for cal in cal:
    filter = ['u','B','g','V','r','i','Y','J','H']
    #filter = ['B','H']
    for filter in filter:        
        file = filter+'_'+cal+'_update3.csv'
        os.system("python H0CSP_vpec.py "+file)

os.system('say "your program has finished."')
