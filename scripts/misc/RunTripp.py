import os


filter = ['u','B','g','V','r','i','Y','J','H']


for filter in filter:
    file = filter+'_ceph.csv'
    #os.system("python H0CSP_noHM.py "+file)
    os.system("python CalibrationTripp.py "+file)


os.system('say "your program has finished."')
