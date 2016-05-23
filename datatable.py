#!/usr/bin/env python

import numpy as np

try:
    from asymmetry import asymmetry
    from unpolxs import unpolxs
    from tools import zload
except:
    from pyg2pasym import asymmetry, unpolxs, zload

for E in [2253.5, 1710.5, 1157.0, 3350.5]:
    xs = unpolxs(E)
    xs.calxs()
    xs.calxsrad()
    xs.save2pkl()
    del xs
    asym = asymmetry(E)
    #asym.data = zload('data/asym_{}.pkl'.format(int(E)))
    asym.calxs()
    asym.calxsrad()
    asym.calxsdf()
    asym.calxsdfrad(0)
    asym.calxsdfrad(1)
    asym.calasym()
    asym.save2pkl()
    del asym
