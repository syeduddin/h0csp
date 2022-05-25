import pandas as pd
import numpy as np

indir = '../../results/'

u = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_u.csv'))
B = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_B.csv'))
g = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_g.csv'))
V = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_V.csv'))
r = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_r.csv'))
i = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_i.csv'))
Y = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_Y.csv'))
J = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_J.csv'))
H = pd.DataFrame(pd.read_csv(indir+'/Ceph_res_H.csv'))


d = ['u','B','g','V','r','i','Y','J','H']


data=u.merge(B,on='sn',how='left',suffixes=['_u','_B']).merge(g,on='sn',how='left',suffixes=['','_g']).merge(V,on='sn',how='left',suffixes=['','_V']).merge(r,on='sn',how='left',suffixes=['','_r']).merge(i,on='sn',how='left',suffixes=['','_i']).merge(Y,on='sn',how='left',suffixes=['','_Y']).merge(J,on='sn',how='left',suffixes=['','_J']).merge(H,on='sn',how='left',suffixes=['','_H'])




data.to_csv('../../results/AllRes.csv')
