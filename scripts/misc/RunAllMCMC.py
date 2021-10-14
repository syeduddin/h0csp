import os


filter = ['u','B','g','V','r','i','Y','J','H']

for filter in filter:
    file = filter+'_trgb.csv'
    os.system("python H0CSPcuts.py "+file)


os.system('say "your program has finished."')
