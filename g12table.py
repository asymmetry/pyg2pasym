#!/usr/bin/env python

import sys
import numpy as np

try:
    from tools import zload
except:
    from pyg2pasym import zload

Ebeam = int(sys.argv[1])
data = zload('data/asym_{}.pkl'.format(Ebeam))
fortg12 = 'data/WQg12_{}.txt'.format(Ebeam)

W=data['W']
Q2=data['Q2GeV']
g1=data['g1']
g2=data['g2']

N=len(W)
g1[np.isnan(g1)]=0
g2[np.isnan(g2)]=0

with open(fortg12, 'w') as f:
    for i in xrange(N):
        f.write('{:7.2f} {:7.5f} {:12.9f} {:12.9f}\n'.format(W[i], Q2[i], g1[i], g2[i]))
