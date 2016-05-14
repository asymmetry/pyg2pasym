#!/usr/bin/env python

import os, re, zlib
import cPickle as pkl
import numpy as np

def zdump(value, filename):
    with open(filename, 'wb', -1) as fpz:
	fpz.write(zlib.compress(pkl.dumps(value, -1), 9))

def zload(filename):
    with open(filename, 'rb') as fpz:
	value  = fpz.read()
	try:
            return pkl.loads(zlib.decompress(value))
	except:
            return pkl.loads(value)

def Q2range():
    Q2_0 = np.asarray([0])
    Q2_1 = np.arange(0.0001, 0.01, 0.0001)
    Q2_2 = np.arange(0.01, 0.1, 0.001)
    Q2_3 = np.arange(0.1, 5.01, 0.01)
    Q2 = np.concatenate((Q2_0, Q2_1, Q2_2, Q2_3))
    return Q2, len(Q2)

def Wrange():
    Wmin, Wmax, Wstep = 1080, 2000, 5
    NW = int((Wmax - Wmin) / Wstep) + 1
    W = np.arange(Wmin, Wmax + Wstep, Wstep)
    return W, NW, Wmin, Wmax, Wstep

# generate hash for W(MeV), Q2(MeV2) pairs
def WQhash(W, Q2):
    return str(hash(np.array_str(W) + np.array_str(Q2)))[-5:]

# generate ids for W(MeV), Q2(GeV2) pairs for bilinear interpolation
def id_WQ2(W, Q2GeV):
    # remove out range value
    W[np.isnan(W)] = 1080
    Q2GeV[np.isnan(Q2GeV)] = 0
    W[W<1080] = 1080
    W[W>2000] = 2000
    Q2GeV[Q2GeV<0] = 0
    Q2GeV[Q2GeV>5] = 5
    # floor, ceil, id for W
    W_floor = np.floor(W * 2 / 10.) * 5
    W_ceil = np.ceil(W * 2 / 10.) * 5
    W_floor_id = ((W_floor - 1080) / 5.).astype('int')
    W_ceil_id = ((W_ceil - 1080) / 5.).astype('int')
    # floor, ceil, id for Q2
    Q2GeV_cp, Q2GeV_floor_cp, Q2GeV_ceil_cp = [0]*4, [0]*4, [0]*4
    Q2GeV_floor_id_cp, Q2GeV_ceil_id_cp = [0]*4, [0]*4

    Q2GeV_cp[0] = np.logical_and(Q2GeV < 0.0001, Q2GeV >= 0)
    Q2GeV_cp[1] = np.logical_and(Q2GeV < 0.01, Q2GeV >= 0.0001)
    Q2GeV_cp[2] = np.logical_and(Q2GeV < 0.1, Q2GeV >= 0.01)
    Q2GeV_cp[3] = np.logical_and(Q2GeV <= 5, Q2GeV >= 0.1)

    Q2GeV_floor_cp[0] = np.zeros(len(Q2GeV))
    Q2GeV_ceil_cp[0] = np.ones(len(Q2GeV)) * 0.0001
    Q2GeV_floor_id_cp[0] = np.zeros(len(Q2GeV), 'int')
    Q2GeV_ceil_id_cp[0] = np.ones(len(Q2GeV), 'int')

    Q2GeV_floor_cp[1] = np.floor(Q2GeV * 1.e4) / 1.e4
    Q2GeV_ceil_cp[1] = np.ceil(Q2GeV * 1.e4) / 1.e4
    Q2GeV_floor_id_cp[1] = ((Q2GeV_floor_cp[1] - 0.0001) / 0.0001 + 1).astype('int')
    Q2GeV_ceil_id_cp[1] = ((Q2GeV_ceil_cp[1] - 0.0001) / 0.0001 + 1).astype('int')

    Q2GeV_floor_cp[2] = np.floor(Q2GeV * 1.e3) / 1.e3
    Q2GeV_ceil_cp[2] = np.ceil(Q2GeV * 1.e3) / 1.e3
    Q2GeV_floor_id_cp[2] = ((Q2GeV_floor_cp[2] - 0.01) / 0.001 + 100).astype('int')
    Q2GeV_ceil_id_cp[2] = ((Q2GeV_ceil_cp[2] - 0.01) / 0.001 + 100).astype('int')

    Q2GeV_floor_cp[3] = np.floor(Q2GeV * 1.e2) / 1.e2
    Q2GeV_ceil_cp[3] = np.ceil(Q2GeV * 1.e2) / 1.e2
    Q2GeV_floor_id_cp[3] = ((Q2GeV_floor_cp[3] - 0.1) / 0.01 + 191).astype('int')
    Q2GeV_ceil_id_cp[3] = ((Q2GeV_ceil_cp[3] - 0.1) / 0.01 + 191).astype('int')

    Q2GeV_floor, Q2GeV_ceil, Q2GeV_floor_id, Q2GeV_ceil_id = 0, 0, 0, 0
    for i in range(4):
        Q2GeV_floor += Q2GeV_floor_cp[i] * Q2GeV_cp[i]
        Q2GeV_ceil += Q2GeV_ceil_cp[i] * Q2GeV_cp[i]
        Q2GeV_floor_id += Q2GeV_floor_id_cp[i] * Q2GeV_cp[i]
        Q2GeV_ceil_id += Q2GeV_ceil_id_cp[i] * Q2GeV_cp[i]

    # four id arrays
    NW = 185
    id_WfQf = Q2GeV_floor_id * NW + W_floor_id
    id_WfQc = Q2GeV_ceil_id * NW + W_floor_id
    id_WcQf = Q2GeV_floor_id * NW + W_ceil_id
    id_WcQc = Q2GeV_ceil_id * NW + W_ceil_id

    return [W_floor, W_ceil], [Q2GeV_floor, Q2GeV_ceil], [id_WfQf, id_WfQc, id_WcQf, id_WcQc]
