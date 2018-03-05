#!/usr/bin/env python

import numpy as np

try:
    from asymmetry import asymmetry
    from unpolxs import unpolxs
    from tools import zload
except:
    from pyg2prad import asymmetry, unpolxs, zload

for E in [2253.5, 1710.5, 1157.0, 3350.5]:
    xs = unpolxs(E)
    xs.calxs()
    xs.save2pkl()
    del xs

    asym = asymmetry(E)
    asym.calxs()
    asym.calxsdf()
    asym.save2pkl()
    del asym

    data = zload('data/asym_{}.pkl'.format(int(E)))
    fortg12 = 'data/WQg12_{}.txt'.format(int(E))

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

for E in [2253.5, 1710.5, 1157.0, 3350.5]:
    xs = unpolxs(E)
    xs.calxs()
    xs.calxsrad()
    xs.save2pkl()
    del xs

    asym = asymmetry(E)
    asym.calxs()
    asym.calxsrad()
    asym.calxsdf()
    asym.calxsdfrad(0)
    asym.calxsdfrad(1)
    asym.calasym()
    asym.save2pkl()
    del asym
