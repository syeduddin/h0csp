import os


filter = ['u','B','g','V','r','i','Y','J','H']

for filter in filter:
    file = filter+'_sbf.csv'
    os.system("python H0CSP.py "+file)


#os.system('say "your program has finished."')
