import os


filter = ['u','B','g','V','r','i','Y','J','H']

for filter in filter:
    file = filter+'_trgb.dat'
    os.system("python H0CSP.py "+file)
