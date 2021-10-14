#!/usr/bin/env python

from astropy.io import ascii
from astropy.table import Table,join
import sys,os
from numpy import *

tab = ascii.read('st_fits_CSPI+II+Ceph+TRGB_gen3.dat')
tab2 = ascii.read('Spreadsheet.csv')

tab2 = tab2['Name','sample','Sub-type','q','Cosmology?','Physics?']
tab2.rename_column('Name','name')
tab2.rename_column('Sub-type','subtype')
tab2.rename_column('q','quality')
tab2['cosmo'] = where(tab2['Cosmology?'] == 'Yes', 1, 0)
tab2['phys'] = where(tab2['Physics?'] == 'Yes', 1, 0)
tab2.remove_column('Cosmology?')
tab2.remove_column('Physics?')
tab = join(tab, tab2, join_type='left', keys='name')

Btab = tab[tab['f'] == 'B']['name','Mmax','eMmax']
Btab.rename_column('Mmax','B')
Btab.rename_column('eMmax','eB')
Vtab = tab[tab['f'] == 'V']['name','Mmax','eMmax']
Vtab.rename_column('Mmax','V')
Vtab.rename_column('eMmax','eV')

tab = join(tab, Btab, keys='name')
tab = join(tab, Vtab, keys='name')
tab['BV'] = tab['B'] - tab['V']
tab['eBV'] = sqrt(power(tab['eB'],2)+power(tab['eV'],2))
tab.remove_columns(['B','V'])
tab['BV'].info.format='%.3f'
tab['eBV'].info.format='%.3f'

for filt in ['u','B','V','g','r','i','Y','J','H']:
   ids = (tab['f'] == filt)
   newtab = tab[ids]
   if filt == 'B':
      newtab['covBV_M'] = power(newtab['eB'],2)
   elif filt == 'V':
      newtab['covBV_M'] = -power(newtab['eV'],2)
   else:
      newtab['covBV_M'] = 0
   newtab['covBV_M'].info.format = '%.6f'
   newtab.remove_columns(['f','eB','eV'])
   newtab.write(filt+"_max.csv", format='ascii.csv', delimiter=',',
         overwrite=True)
