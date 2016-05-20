#!/usr/bin/env python

import os, re
import numpy as np
import scipy.special,scipy.integrate
import subprocess as sp
from multiprocessing.dummy import Pool as tp
from os.path import dirname, exists, join, realpath

from tools import id_WQ2, Q2range, WQhash, Wrange, zload, zdump
from unpolxs import unpolxs

# 1001.3898v1 1
# PhysRevD.59.094021 2
# PhysRevD.58.112003 3
# arxiv:hep-ex/9903055 4
# Physics Reports.378.99 5
# constant from pdg 2010

class asym(unpolxs):
    def __init__(self, Ebeam=2253):
        # constant
        unpolxs.__init__(self, Ebeam)
        self.datadir = 'data'
        self.save = join(self.datadir, 'asym_{}.pkl'.format(int(self.E))) # pkl filename for calculated data
        self.datatype = '2007tot'
        self.DMAID = zload(join(self.datadir, 'MAID1D{}.pdt'.format(self.datatype)))
        self.polin = join(dirname(realpath(__file__)), 'polin')
        self.polel = join(dirname(realpath(__file__)), 'polel')

    def calasym(self, nodump=False):
        if not self.data:
            self.calxsdfrad()
            self.calxsrad()

        if not 'dxs11_rad' in self.data:
            self.calxsrad()
        if not 'dxsL11_rad' in self.data:
            self.calxsdfrad(0)
        if not 'dxsT11_rad' in self.data:
            self.calxsdfrad(1)

        self.data['asymL11'] = self.data['dxsL11'] / (2 * self.data['dxs11'])
        self.data['asymT11'] = self.data['dxsT11'] / (2 * self.data['dxs11'])

        self.data['asymL11_rad'] = self.data['dxsL11_rad'] / (2 * self.data['dxs11_rad'])
        self.data['asymT11_rad'] = self.data['dxsT11_rad'] / (2 * self.data['dxs11_rad'])

    # do not call after choose_WQ2
    def calxsdf(self):
        # choose chan and sum
        ma = self.mp
        chan = [0, 2] # proton
        chanfilt = [self.DMAID['chan'] == chan[0], self.DMAID['chan'] == chan[1]]

        if not self.data:
            self.data = {}
        for k in ["sigTT'", "sigLT'", "sigL", "sigT"]:
            self.data[k] = self.DMAID[k][chanfilt[0]] + self.DMAID[k][chanfilt[1]] # ub
        for k in ['W', 'Q2', 'w(lab)']:
            self.data[k] = self.DMAID[k][chanfilt[0]] # MeV, GeV, MeV
        self.data['nu'] = self.data['w(lab)']
        self.data['sigTT'] = self.data["sigTT'"]
        self.data['sigLT'] = self.data["sigLT'"]
        del self.data["sigTT'"], self.data["sigLT'"], self.data['w(lab)']
        # Q2 GeV to MeV
        self.data['Q2GeV'] = self.data['Q2']
        self.data['Q2'] = self.data['Q2'] * 1e6 # GeV to MeV
        self.ohash = WQhash(self.data['W'], self.data['Q2']) # hash for original data

        # g1 g2 f1 f2
        K = (self.data['W']**2 - ma**2) / (2 * ma) # (1,38),(5,67) MeV
        gamma = np.sqrt(self.data['Q2']) / self.data['nu'] # (1,24)
        gamma2 = 1 + gamma**2
        Par = 4 * np.pi**2 * self.alpha / ma * self.hc2
        Parnu = 4 * np.pi**2 * self.alpha / self.data['nu'] * self.hc2
        y = self.data['nu'] / self.E
        self.data['x'] = self.data['Q2'] / (2 * ma * self.data['nu'])
        self.data['g1'] = K / (Par * gamma2) * (self.data['sigTT'] + gamma * self.data['sigLT']) # (5,73)
        self.data['g2'] = -K / (Par * gamma2) * (self.data['sigTT'] - self.data['sigLT'] / gamma) # (5,73)
        #self.data['F1'] = self.data['sigT'] * K / Par # (5,73)
        #self.data['F2'] = (self.data['sigT'] + self.data['sigL']) * K * gamma**2 / (Parnu * gamma2) # (5,73)

        # cross section differences
        E2 = self.E - self.data['nu']
        self.data['theta'] = 2 * np.arcsin(np.sqrt(self.data['Q2'] / (4 * E2 * self.E))) # scat angle, rad
        sin2thd2, cos2thd2, tan2thd2 = self.trigo(self.data['theta'])
        mott = self._mott(self.E, self.data['theta'])
        alphaEQ = 4 * self.alpha**2 * E2 * self.hc2 / (self.data['nu'] * self.E * self.data['Q2'] * ma)
        self.data['dxsL11'] = alphaEQ * ((self.E + E2 * np.cos(self.data['theta'])) * self.data['g1'] - 2 * ma * self.data['x'] * self.data['g2']) # (4,1)
        self.data['dxsT11'] = alphaEQ * E2 * np.sin(self.data['theta']) * (self.data['g1'] + 2*self.E * self.data['g2'] / self.data['nu']) # (4,2)

    # radiate polarized dxsL or dxsT, LT = 0 for long, 1 for tran
    def calxsdfrad(self, LT, nodump=False):
        if not self.data:
            self.calxsdf()

        xskey = 'dxsL11' if LT == 0 else 'dxsT11'
        xsin = self.inelastic(xskey, 1, 1, pol=1, nodump=nodump) # inelastic extern
        xsin_in = self.polin_intern(xskey, xsin, nodump=nodump)
        xsel_ex, tmp = self.elastic(xskey, 1, 1, pol=1) # elastic extern
        xsel_in = self.polel_intern(xskey)
        self.data['{}_rad'.format(xskey)] = 2 * (xsel_ex + xsel_in) + xsin_in
        self.data['{}_rad_in'.format(xskey)] = xsin
        self.data['{}_rad_el'.format(xskey)] = 2 * (xsel_in + xsel_ex)
        self.data['{}_rad_in_in'.format(xskey)] = xsin_in
        self.data['{}_rad_el_ex'.format(xskey)] = xsel_ex
        self.data['{}_rad_el_in'.format(xskey)] = xsel_in

    # polarized radiative correction for inelastic internal
    # xsin is the xs after external rc
    def polin_intern(self, xskey, xsin, nodump=False):
        LT = 3 if 'T' in xskey else 1

        fname = 'dxs'
        ohash = WQhash(self.data['W'], self.data['Q2'])
        pklfile = join(self.datadir, 'pol/polin{}{}.pkl'.format(fname, ohash))
        if not nodump and exists(pklfile):
            return zload(pklfile)

        Ethlist = 'polinEiEflist{}.dat'.format(fname)
        outlist = 'polin{}.dat'.format(fname)
        E2 = self.E - self.data['nu']
        N = len(E2)
        theta=self.data['theta']
        theta[np.isnan(theta)] = 0
        xsin[np.isnan(xsin)] = 0
        with open(Ethlist, 'w') as f:
            for i in xrange(N):
                if i % 10000 == 0:
                    print 'save {} {}'.format(Ethlist, i)
                f.write('{:7.2f} {:7.2f} {:7.5f} {:7.5f}\n'.format(self.E, E2[i], theta[i], xsin[i] * 1000 / 2.))
        command = [self.polin, str(self.E), str(LT), fname]
        sp.call(command)

        xs_rad=np.zeros(N)
        # read value
        with open(outlist, 'r') as f:
            i = 0
            for l in f:
                xs_rad[i] = -float(l) * 2. / 1000
                i += 1
        os.remove(outlist)
        os.remove(Ethlist)

        if not nodump:
            if not exists(join(self.datadir, 'pol')):
                os.makedirs(join(self.datadir, 'pol'))
            zdump(xs_rad, pklfile)

        return xs_rad

    # polarized radiative correction for elastic tail internal
    def polel_intern(self, xskey, nodump=False):
        LT = 3 if 'T' in xskey else 1

        fname = 'dxs'
        ohash = WQhash(self.data['W'], self.data['Q2'])
        pklfile = join(self.datadir, 'pol/polel{}{}.pkl'.format(fname, ohash))
        if not nodump and exists(pklfile):
            return zload(pklfile)

        Ethlist = 'polelEiEflist{}.dat'.format(fname)
        outlist = 'polel{}.dat'.format(fname)
        E2 = self.E - self.data['nu']
        N = len(E2)
        theta = self.data['theta']
        theta[np.isnan(theta)] = 0
        with open(Ethlist, 'w') as f:
            for i in xrange(N):
                if i % 10000 == 0:
                    print 'save {} {}'.format(Ethlist, i)
                f.write('{:7.2f} {:7.2f} {:7.5f}\n'.format(self.E, E2[i], theta[i]))
        command = [self.polel, str(LT), fname]
        sp.call(command)

        xs_rad = np.zeros(N)
        # read value
        with open(outlist, 'r') as f:
            i = 0
            for l in f:
                tmp = [x for x in re.split('\s',l) if len(x) > 0]
                xs_rad[i] = -np.float64(tmp[1]) / 1.e3
                i += 1
        os.remove(outlist)
        os.remove(Ethlist)

        if not nodump:
            if not exists(join(self.datadir, 'pol')):
                os.makedirs(join(self.datadir, 'pol'))
            zdump(xs_rad, pklfile)

        return xs_rad
