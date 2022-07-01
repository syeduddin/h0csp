from astropy.io import ascii
from astropy.table import join

t1 = ascii.read('trgbdist_for_chris_2021-2.csv')
t2 = ascii.read('calibrators_trgb.csv')

t = join(t1,t2,keys='sn')

print (t['muTRGBF20']-t['dist'])
