import os

cal = ['ceph','trgb','sbfcombined']
#cal = ['sbfcombined']

for cal in cal:
    filter = ['u','B','g','V','r','i','Y','J','H']
    #filter = ['B','H']
    for filter in filter:        
        file = filter+'_'+cal+'_update2.csv'
        os.system("python H0CSPcuts.py "+file)

os.system('say "your program has finished."')
