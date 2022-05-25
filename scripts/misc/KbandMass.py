import numpy as np
import sys


k = float(sys.argv[1])
ek = float(sys.argv[2])
d = float(sys.argv[3])
ed = float(sys.argv[4])

#mass = (-0.4*(k-mu))-1.04

#emass = np.log(10)*mass*(ek**2+emu**2)
#emass = np.sqrt(emass)
#print ('%.3f'%mass,',%.3f'%emass)

#K−5log10D−25−0.11AV

M_k =k - (5*np.log10(d)) -25

#log10(M∗) = 10.58 − 0.44(MK + 23) .

mass = 10.58 -(0.44*(M_k+23))

emass = ek**2 + ((ed/d*np.log(10))**2)


emass = np.sqrt(emass)
print ('%.3f'%mass,',%.3f'%emass)

