#!/usr/bin/env python

import os, re
import numpy as np
import scipy.special,scipy.integrate
import subprocess as sp
from multiprocessing.dummy import Pool as tp
from os.path import dirname, exists, join, realpath

from tools import id_WQ2, Q2range, WQhash, Wrange, zload, zdump

# 1001.3898v1 1
# PhysRevD.59.094021 2
# PhysRevD.58.112003 3
# arxiv:hep-ex/9903055 4
# Physics Reports.378.99 5
# constant from pdg 2010

class unpolxs():
    def __init__(self, Ebeam=2253):
        # constant
        self.alpha = 1 / 137.035999679 # fine-structure constant
        self.alphapi = self.alpha / np.pi
        self.mp = 938.272013 # mass of proton, MeV
        self.mn = 939.565346 # mass of neutron, MeV
        self.me = 0.51099891 # mass of electron, MeV
        self.mup = 2.792847356 # magnetic moment for proton
        self.E = Ebeam # energy of electron, MeV
        self.hc = 197.3269631 # MeV fm
        self.hc2 = 0.389379304e9 # conversion constant, MeV^2ub
        self.allQ2, self.NQ2 = Q2range()
        self.allW, self.NW, self.Wmin, self.Wmax, self.Wstep = Wrange()
        self.datadir = 'data'
        self.save = join(self.datadir, 'unpolxs_{}.pkl'.format(int(self.E))) # pkl filename for calculated data
        self.data = False
        self.Ndata = 126170 # datas in raw database
        self.relatedid = np.asarray([], 'int') # used for data radiate speed up, only care about the related data
        self.pbF12 = join(dirname(realpath(__file__)), 'pbF12')
        # default setting
        #self.ta = 0.00694 + 0.02805/2. # rad length in the target before
        #self.tb = 0.00956 + 0.02805/2. # rad length in the target after
        self.ta = 0.02511 + 0.04138/2. * 0
        self.tb = 0.02249 + 0.04138/2. * 0

    # do not call after choose_WQ2
    def calxs(self):
        M = self.mp

        data = {}
        data['W'] = np.zeros(self.NQ2 * self.NW, dtype = 'float32')
        data['Q2'] = np.zeros(self.NQ2 * self.NW, dtype = 'float32')
        data['nu'] = np.zeros(self.NQ2 * self.NW, dtype = 'float32')
        ndata = 0
        for Q2 in self.allQ2:
            for W in self.allW:
                data['Q2'][ndata] = Q2 * 1e6 # GeV to MeV
                data['W'][ndata] = W
                data['nu'][ndata] = (W**2 + Q2 * 1e6 - M**2) / (2 * M)
                ndata += 1
        self.ohash = WQhash(data['W'], data['Q2']) # hash for original data

        E2 = self.E - data['nu']
        data['theta'] = 2 * np.arcsin(np.sqrt(data['Q2'] / (4 * E2 * self.E))) # scat angle, rad

        # PB model for H1, He4, C12, N14, Al27
        self.data = data
        for ZA in [[1, 1], [2, 4], [6, 12], [7, 14], [13, 27]]:
            pbdxs = self.xspb(self.E, E2, data['theta'], ZA[0], ZA[1])
            data['dxs{}{}'.format(ZA[0], ZA[1])] = pbdxs
        self.data = data

        return data

    def calxsrad(self):
        if not self.data:
            self.calxs()

        xskey = 'dxs{}{}'
        for ZA in [[1, 1], [2, 4], [6, 12], [7, 14], [13, 27]]:
            xsin = self.inelastic(xskey.format(ZA[0], ZA[1]), ZA[0], ZA[1])
            xsel_ex, xsel_ep = self.elastic(xskey.format(ZA[0], ZA[1]), ZA[0], ZA[1])
            self.data['dxs{}{}_rad'.format(ZA[0], ZA[1])] = xsin + xsel_ex + xsel_ep
            self.data['dxs{}{}_rad_in'.format(ZA[0], ZA[1])] = xsin
            self.data['dxs{}{}_rad_el'.format(ZA[0], ZA[1])] = xsel_ex + xsel_ep

    def save2pkl(self, filename=False):
        if not filename:
            zdump(self.data, self.save)
        else:
            zdump(self.data, join(self.datadir, filename))

    # choose data, input is W(MeV), Q2(MeV2) array. if no keys will overwrite self.data
    def choose_WQ2(self, W, Q2, keys, nodump=False):
        if not self.data:
            return False
        Q2GeV = Q2 / 1.e6
        Ws, Q2s, ids = id_WQ2(W, Q2GeV)

        dhash = WQhash(W, Q2)
        addkeys = ''
        if keys:
            addkeys = '_' + ''.join(keys)
        pklpath = join(self.datadir, 'tmp/unpolxs{}_{}.pkl'.format(addkeys, dhash))
        if not nodump and exists(pklpath):
            data = zload(pklpath)
            if not keys:
                self.data = data
                return True
            else:
                return data

        if len(self.data['W']) != self.Ndata:
            odata = zload(self.save) # reload defaults
        else:
            odata = self.data

        keysbak = keys
        if not keys:
            keys = odata.keys()

        data = {}
        # speed up radiation
        self.relatedid = np.concatenate((self.relatedid, ids[0], ids[1], ids[2], ids[3]))
        self.relatedid = np.unique(self.relatedid)

        for k in keys:
            # bilinear interpolation
            d00 = odata[k][ids[0]]
            d01 = odata[k][ids[1]]
            d10 = odata[k][ids[2]]
            d11 = odata[k][ids[3]]
            Weq = (Ws[0] == Ws[1])
            d0 = (Ws[1] - W) / (Ws[1] - Ws[0]) * d00 + (W - Ws[0]) / (Ws[1] - Ws[0]) * d10
            d1 = (Ws[1] - W) / (Ws[1] - Ws[0]) * d01 + (W - Ws[0]) / (Ws[1] - Ws[0]) * d11
            d0[Weq] = d00[Weq]
            d1[Weq] = d01[Weq]
            Q2eq = (Q2s[1] == Q2s[0])
            d = (Q2s[1] - Q2GeV) / (Q2s[1] - Q2s[0]) * d0 + (Q2GeV - Q2s[0]) / (Q2s[1] - Q2s[0]) * d1
            d[Q2eq] = d0[Q2eq]
            data[k] = d
        del odata

        if not nodump:
            zdump(data, pklpath)
        if not keysbak:
            self.data = data
            return True
        else:
            return data

    # choose data, input is W(MeV), theta(rad) array
    def choose_Wth(self, W, theta, keys=False, nodump=False):
        sin2thd2 = np.sin(theta / 2.)**2
        E2 = (W**2 - self.mp**2 - 2 * self.mp * self.E) / (-2 * self.mp - 4 * self.E * sin2thd2)
        Q2 = 4 * self.E * E2 * sin2thd2
        return self.choose_WQ2(W, Q2, keys, nodump=nodump)

    # choose data, input is nu(MeV), theta(rad) array
    def choose_nuth(self, nu, theta, keys=False, nodump=False):
        E2 = self.E - nu
        sin2thd2 = np.sin(theta / 2.)**2
        Q2 = 4 * self.E * E2 * sin2thd2
        W = np.sqrt(self.mp**2 + 2 * self.mp * nu - Q2)
        return self.choose_WQ2(W, Q2, keys, nodump=nodump)

    # choose data, input is nu(MeV), Q2(MeV2) array
    def choose_nuQ2(self, nu, Q2, keys=False, nodump=False):
        W = np.sqrt(self.mp**2 + 2 * self.mp * nu - Q2)
        return self.choose_WQ2(W, Q2, keys, nodump=nodump)

    # mass of target
    def mass(self, Z, A):
        MA = self.mp / 1.007276 # mass of nucleon, MeV
        if Z == 1 and A == 1:
            MT = MA * 1.007276 # H-1
        elif Z == 2 and A == 4:
            MT = MA * 4.002602 # He-4
        elif Z == 6 and A == 12:
            MT = MA * 12.011 # C-12
        elif Z == 7 and A == 14:
            MT = MA * 14.007 # N-14
        else:
            MT = MA * A
        return MT

    def trigo(self, theta):
        sin2thd2 = np.sin(theta / 2.)**2
        cos2thd2 = np.cos(theta / 2.)**2
        tan2thd2 = np.tan(theta / 2.)**2
        return sin2thd2, cos2thd2, tan2thd2

    # F1 F2 from PBosted model
    def __pbf12(self, Q2, W, Z, A, fname='', dumpfile=True):
        ohash = WQhash(self.data['W'], self.data['Q2'])
        FDPB = join(self.datadir, 'pb/PB{}{}{}_{}.pkl'.format(Z, A, fname, ohash))

        PBdata = {}
        Q2GeV = Q2 / 1.e6
        if exists(FDPB):
            DPB = zload(FDPB)
            PBdata['PBF1_q'] = DPB['PBF1_q'] # PBosted F1, quasi
            PBdata['PBF2_q'] = DPB['PBF2_q'] # PBosted F2, quasi
            PBdata['PBF1_i'] = DPB['PBF1_i'] # PBosted F1, inelastic
            PBdata['PBF2_i'] = DPB['PBF2_i'] # PBosted F2, inelastic
        else:
            N = len(Q2GeV)
            for k in ['PBF1_i', 'PBF2_i', 'PBF1_q', 'PBF2_q']:
                PBdata[k] = np.zeros(N, 'float64')
            WQ2list = 'WQ2list{}{}{}.dat'.format(Z, A, fname)
            PBF12 = 'PBF12{}{}{}.dat'.format(Z, A, fname)
            W2 = (W / 1000.)**2
            iq=['q', 'i']

            # write W Q2 list
            with open(WQ2list, 'w') as f:
                for i in xrange(N):
                    if i%10000 == 0:
                        print 'save {} {} ......'.format(WQ2list, i)
                    f.write('{} {}\n'.format(W2[i], Q2GeV[i]))

            for j in range(2):
                command = [self.pbF12, str(j), str(Z), str(A), fname]
                sp.call(command)

                # read F1 F2
                with open(PBF12, 'r') as f:
                    i = 0
                    for l in f:
                        f12 = re.split('\s', l)
                        PBdata['PBF1_{}'.format(iq[j])][i] = np.float64(f12[0])
                        PBdata['PBF2_{}'.format(iq[j])][i] = np.float64(f12[1])
                        i += 1
                os.remove(PBF12)

            if dumpfile:
                if not exists(join(self.datadir, 'pb')):
                    os.makedirs(join(self.datadir, 'pb'))
                zdump({'PBF1_i': PBdata['PBF1_i'], 'PBF2_i': PBdata['PBF2_i'], 'PBF1_q': PBdata['PBF1_q'], 'PBF2_q':PBdata['PBF2_q']}, FDPB)

            os.remove(WQ2list)

        return PBdata

    # mott cross section
    def _mott(self, E, theta):
        return self.alpha**2 * self.hc2 * np.cos(theta / 2.)**2 / (4 * E**2 * np.sin(theta / 2.)**4) # (1,17)

    # unpol cross section
    def __unpol(self, mott, F1, F2, nu, tan2thd2, M):
        return mott * (2. / M * F1 * tan2thd2 + F2 / nu) # (1,16)

    # pbosted cross section
    def xspb(self, E, E2, theta, Z, A, filename='', dumpfile=True):
        M = self.mp
        nu = E - E2
        sin2thd2, cos2thd2, tan2thd2 = self.trigo(theta)
        Q2 = 4 * E * E2 * sin2thd2
        W = np.sqrt(M**2 + 2 * M * nu - Q2)
        mott = self._mott(E, theta)
        pbdata = self.__pbf12(Q2, W, Z, A, filename, dumpfile)
        dxs_i = self.__unpol(mott, pbdata['PBF1_i'], pbdata['PBF2_i'], nu, tan2thd2, M)
        dxs_q = self.__unpol(mott, pbdata['PBF1_q'], pbdata['PBF2_q'], nu, tan2thd2, M)
        del pbdata
        return dxs_i + dxs_q

    # elastic form factor,W1 and W2 in (1,A3)
    # return GE,GM if GEM is true
    def __elform(self, Q2, Z, A, M, GEM=False):
        if Z == 1 and A == 1:
            # (6,1,tableII)
            Q2GeV = Q2 / 1e6
            pGE = [2.94, 3.04, -2.255, 2.002, -0.5338, 4.875e-2]
            pGM = [3, 1.39, 0.122, -8.34e-3, 4.25e-4, -7.79e-6]
            Ge, Gm = 1, 1
            for i in range(6):
                Ge += pGE[i] * Q2GeV**(i+1)
                Gm += pGM[i] * Q2GeV**(i+1)
            Ge, Gm = 1 / Ge, self.mup / Gm
            # sub_rtail.f
            #fd = lambda qms, a: 1. / (1 + qms / a**2)**2
            #q = np.sqrt(Q2GeV)
            #Gm = 1 + 0.35 * q + 2.44 * Q2GeV + 0.5 * q * Q2GeV
            #Gm = self.mup / (Gm + 1.04 * Q2GeV**2 + 0.34 * q * Q2GeV**2)
            #Ge = 1 / (1 + 0.62 * q + 0.68 * Q2GeV + 2.8 * q * Q2GeV + 0.84 * Q2GeV**2)
        elif Z == 2 and A == 4:
            a, b = 0.316, 0.675 # (4,below 10)
            Q2fm = Q2 / (self.hc**2) # fm^-2
            Ge = (1 - (a**2 * Q2fm)**6) * np.exp(-b**2 * Q2fm) * Z # (4,10), normalized by Z
            Gm = 0 # sub_rtail.f
        elif Z == 6 and A == 12:
            Q2fm = Q2 / (self.hc**2) # fm^-2
            Q2range1 = (Q2fm <= 3.2)
            Q2range2 = (Q2fm > 3.5)
            alpha = (Z - 2) / 3. # (5,17)
            Q2range3 = np.logical_not(np.logical_or(Q2range1, Q2range2))
            a = 1.64 * Q2range1.astype('float32') + 1.68 * Q2range2.astype('float32') # (5,below 18)
            a /= self.hc
            F = (1. - alpha / (2 * (2 + 3. * alpha)) * Q2 * a * a) * np.exp(-Q2 * a * a / 4.)
            F = np.sqrt(1e-5) * Q2range3.astype('float32') + F * np.logical_or(Q2range1, Q2range2).astype('float32') # (5,18)
            return Z * Z * F * F, 0 # (5,5)
        elif Z == 7 and A == 14:
            # sub_rtail.f
            Q2GeV = Q2 / 1.e6
            hcGeV = self.hc / 1.e3 # GeV fm
            B = 1.75 / hcGeV
            AP = 0.63 / hcGeV
            ALPHA = 0.44
            MU = 30
            X = 0.25 * Q2GeV * B * B
            D = 0.25 * Q2GeV * (AP * AP - B * B / A)
            C = 1. - 2. / 3. * (1. - ALPHA) * X
            FT2 = 2. / 3. * 2. * X / (B * B * (self.mp / 1e3)**2) * (MU * MU)/(Z * Z) * (C * C) * np.exp(-2.*(X + D))
            Q = 1.52 / (hcGeV**2)
            C = 1 - 10. / 21. * X
            FL2 = (C * C) * np.exp(-2. * (X + D)) + (Q2GeV * Q2GeV) / 180. * 10. * (Q * Q) / (Z * Z) * np.exp(-2. * (X + D))
            tau = Q2 / (4 * M**2) # (1,A7)
            Ge = Z * np.sqrt(tau * FT2 + (1 + tau) * FL2)
            Gm = Z * np.sqrt(0.5 * FT2 / tau)
        else:
            # Al, Cu, Au
            Q2fm = Q2 / (self.hc**2)
            b, c = 2.4, 1.07 * np.power(A, 1. / 3.)
            F = 1 / (1 + 1. / 6 * Q2fm * c**2) * np.exp(-1. / 6 * Q2fm * b**2) # (1,A18)
            return (Z * F)**2, 0
        if GEM:
            return Ge, Gm
        tau = Q2 / (4 * M**2) # (1,A7)
        W1 = tau * Gm**2
        W2 = (Ge**2 + tau * Gm**2) / (1 + tau) # (1,A6)
        return W2, W1

    # elastic xs, return F*xs # (1,A55)
    def xsel(self, E, theta, Z, A):
        M = self.mass(Z, A)
        eta, b, T, epsilon = self.__fetepbT(Z)
        sin2thd2, cos2thd2, tan2thd2 = self.trigo(theta)
        Eel = E / (1 + (2 * E / M) * sin2thd2) # elastic peak, (1,A4)
        mott = self._mott(E, theta)
        Q2 = 4 * E * Eel * sin2thd2
        W2, W1 = self.__elform(Q2, Z, A, M)
        spence = scipy.special.spence(cos2thd2) # (1,A48)
        xs = mott * Eel / E * (W2 + 2 * tan2thd2 * W1) # (1,A3)
        F = self.__fVF(Q2, E, Eel, b, T, spence)
        return F * xs

    # elastic asym, LT = 0 long, 1 tran. only for proton
    # hacked from sub_rtail.f
    def asymel(self, E, theta, LT):
        M = self.mp
        me = self.me
        sin2thd2, cos2thd2, tan2thd2 = self.trigo(theta)
        E2 = E / (1 + (2 * E / M) * sin2thd2) # elastic peak, (1,A4)
        me2 = me * me
        nu = E - E2
        M2 = M * M
        Q2 = 2 * M * nu
        S = 2 * M * np.sqrt(E * E + me2) # (7,S)
        lambdas = S * S - 4 * me2 * M2 # (7,6)
        GE, GM = self.__elform(Q2, 1, 1, M, True)
        #GE = 1.2742 / (1. + Q2GeV / 0.6394**2) - .2742 / (1. + Q2GeV / 1.582**2)
        #GM = (1.3262 / (1. + Q2GeV / 0.6397**2) - .3262 / (1. + Q2GeV / 1.3137**2)) * 2.7921 # GE GM from sub_rtails
        tau = Q2 / (4 * M2)
        F1 = 4 * tau * M2 * GM * GM # (7,11)
        F2 = 4 * M2 * (GE * GE + tau * GM * GM) / (1 + tau) # (7,11)
        F3 = -2 * M2 * GE * GM # (7,18)
        F4 = -M2 * GM * (GE - GM) / (1 + tau) # (7,18)
        X = S-Q2 # (7,below 20), then Sx = S - X = Q2
        Sx = Q2 # (7,above 20), Sx = S - X
        lambdaq = Sx * Sx + 4 * M2 * Q2 # (7,below 19)
        lambda0  = S * X * Q2 - me2 * lambdaq - M2 * Q2 * Q2 # (7,below 19,lambda)
        sqlq = np.sqrt(lambdaq)
        sql0 = np.sqrt(lambda0)
        sqls = np.sqrt(lambdas)
        ym = Q2 + 2 * me2
        a_s = S / (2 * me * sqls)
        b_s = 0
        c_s = -me / sqls
        if LT == 0:
            ae, be, ce = M / sqls, 0, -S / (2 * M * sqls)
        else:
            ae = (-S * X + 2 * M2 * ym) / (2 * sql0 * sqls)
            be = sqls / (2 * sql0)
            ce = -(S * Q2 + 2 * me2 * Sx) / (2 * sqls * sql0)
        apq = -Q2 * (ae - be) + ce * Sx
        apn = (Q2 + 4 * me2) * (ae + be) + ce * (S + X)
        dk2ks = a_s * ym + 2 * me2 * b_s + c_s * X
        dksp1 = a_s * S + b_s * X + c_s * 2 * M2
        dapks = 2 * (2 * me2 * (a_s * ae + b_s * be) + 2 * M2 * c_s * ce + ym * (a_s * be + b_s * ae) + S * (a_s * ce + c_s * ae) + X * (b_s * ce + c_s * be))
        thB1 = Q2-2 * me2 # (7,12)
        thB2 = (S * X - M2 * Q2)/(2 * M2) # (7,13)
        thB3 = (2. * (apq * dk2ks - dapks * Q2) * me) / M # (7,21)
        thB4 = apq * Q2 * me * (2. * dksp1 - dk2ks) / (M2 * M) # (7,21)
        sig0 = 2 * np.pi * self.alpha ** 2 / (S * S * Q2 * Q2)
        sig_unpol = sig0 * (thB1 * F1 + thB2 * F2) # (7,14h)
        sig_pol = sig0 * (thB3 * F3 + thB4 * F4)
        return sig_pol / sig_unpol

    # inelastic radiate correction
    # E, E2: energy of in, out electron (MeV), theta: scat angle (rad)
    def inelastic(self, xskey='', Z=1, A=1, pol=0, nodump=False):
        if pol == 1 and A != 1:
            print 'polarized radiative correction only work for proton!'

        MT = self.mass(Z, A)
        eta, b, T, epsilon = self.__fetepbT(Z)
        dE = 10 # (1,A83m),(2,3)
        Q2, nu, xs = self.__kins(xskey, Z, A)
        theta = self.data['theta']
        sin2thd2, cos2thd2, tan2thd2 = self.trigo(theta)
        spence = scipy.special.spence(cos2thd2) # (1,A48)

        # 1st term (1,A82)
        E2 = self.E - nu
        R = self.__fR(self.E, E2, MT, sin2thd2)
        tr = self.__fVtr(Q2, b, pol)
        F = self.__fVF(Q2, self.E, E2, b, T, spence, pol)
        xsrad_11 = np.power(R * dE / self.E, b * (self.tb + tr))
        xsrad_12 = np.power(dE / E2, b * (self.ta + tr))
        xsrad_13 = 1 - epsilon / dE / (1 - b * (T + 2 * tr))
        xsrad_1 = xsrad_11 * xsrad_12 * xsrad_13 * F * xs

        # 2nd/3rd term, integration
        def inte(Emin, Emax, Ea, Eb, rev):
            Nstep = 80
            Estep = (Emax - Emin) / float(Nstep)
            def inteone(i):
                #print 'inelastic radiate integral {} {}'.fromat(i, rev)
                Ep = Emin + i * Estep
                VQ2 = 4 * Ep * Eb * sin2thd2
                dEp = (Ep - Ea) * rev
                tr = self.__fVtr(VQ2, b, pol)
                VF = self.__fVF(VQ2, self.E, E2, b, T, spence, pol)

                if rev < 0:
                    #F = self.__fVF(VQ2, Ep, Eb, b, T, spence, pol)
                    dEpp = dEp / Ea
                    tar, tbr = self.tb + tr, self.ta + tr
                    #VR = self.__fR(Ep, Eb, MT, sin2thd2)
                    EI2 = np.power(dEp / (Eb * R), b * tbr)
                    if pol == 0:
                        xs = self.xspb(Ep, Eb, theta, Z, A, '{}{}'.format(i, rev))
                    else:
                        xs = self.choose_nuth(Ep - Eb, theta, [xskey], nodump=True)[xskey]
                else:
                    #F = self.__fVF(VQ2, Eb, Ep, b, T, spence, pol)
                    dEpp = dEp / Ep
                    tar, tbr = self.ta + tr, self.tb + tr
                    #VR = self.__fR(Eb, Ep, MT, sin2thd2)
                    EI2 = np.power(dEp * R / Eb, b * tbr)
                    if pol == 0:
                        xs = self.xspb(Eb, Ep, theta, Z, A, '{}{}'.format(i, rev))
                    else:
                        xs = self.choose_nuth(Eb - Ep, theta, [xskey], nodump=True)[xskey]

                EI1 = np.power(dEpp, b * tar)
                EI3 = self.__fphiv(dEpp) * b * tar / dEp + epsilon / (2 * dEp**2)
                xsin_inte = EI1 * EI2 * EI3 * xs * VF
                del Ep, VQ2, dEp, tr, VF, dEpp, tar, tbr, EI2, EI1, EI3, xs
                return xsin_inte

            pool = tp(processes = 4)
            intes = pool.map(inteone, range(Nstep + 1))
            #intes = []
            #for i in range(Nstep + 1):
            #    intes.append(inteone(i))
            interad = 0
            for i in range(Nstep):
                interad += Estep * (intes[i] + intes[i + 1]) / 2.
            return interad

        # (2,13,14)
        xsrad_2 = inte(E2 / (1 - 2 * E2 * sin2thd2 / MT), self.E - R * dE, self.E, E2, -1)
        xsrad_3 = inte(E2 + dE, self.E / (1 + 2 * self.E * sin2thd2 / MT), E2, self.E, 1)
        return xsrad_1 + xsrad_2 + xsrad_3

    # elastic tail, extern bremmsstrahlung, and peaking approximation for internal bremmstrahlung (unpol only)
    def elastic(self, xskey='', Z=1, A=1, pol=0):
        if pol == 1 and A != 1:
            print 'polarized radiative correction only work for proton!'
        LT = 1 if 'T' in xskey else 0

        MT = self.mass(Z, A)
        eta, b, T, epsilon = self.__fetepbT(Z)
        Q2, nu, tmp = self.__kins(xskey, Z, A)
        theta = self.data['theta']
        E2 = self.E - nu
        sin2thd2, cos2thd2, tan2thd2 = self.trigo(theta)
        omegas = self.E - E2 / (1 - (2 * E2 / MT) * sin2thd2) # (1,A50)
        omegap = self.E / (1 + (2 * self.E / MT) * sin2thd2) - E2 # (1,A51)
        nus = omegas / self.E # (1,A53)
        nup = omegap / (E2 + omegap) # (1,A53)
        xsel_s = self.xsel(self.E, theta, Z, A)
        xsel_sm = self.xsel(self.E - omegas, theta, Z, A)
        if pol == 1:
            # elastic differential cross section = asym * xs_unpol
            xsel_s *= self.asymel(self.E, theta, LT)
            xsel_sm *= self.asymel(self. E - omegas, theta, LT)
        bphos = b * self.__fphiv(nus) / omegas
        bphop = b * self.__fphiv(nup) / omegap

        term11 = (MT + 2 * (self.E - omegas) * sin2thd2) / (MT - 2 * E2 * sin2thd2)
        term12 = bphos * self.tb + epsilon / 2 / omegas**2
        term2 = bphop * self.ta + epsilon/ 2 / omegap**2
        xsb = term11 * xsel_sm * term12 + xsel_s * term2 # extern bremmsstrahlung (1,A49)
        tr = self.__fVtr(Q2,b)
        if pol == 0:
            xsp = term11 * xsel_sm * bphos * tr + xsel_s * bphop * tr # peaking approximation (1,A56)
        else:
            xsp = 0
        Fsoft = np.power(omegas / self.E, b * (self.tb + tr)) * np.power(omegap / (E2 + omegap), b * (self.ta + tr)) # Multiple-photon correction #(1,A58)
        return xsb * Fsoft, xsp * Fsoft

    # Q2, nu, xs
    def __kins(self, xskey='', Z=1, A=1):
        Q2 = self.data['Q2']
        nu = self.data['nu']
        if xskey == '':
            xs = self.data['dxs{}{}'.format(Z, A)]
        else:
            xs = self.data[xskey]
        return Q2, nu, xs

    # eta, b, epsilon, T from (1,A45-47,52)
    def __fetepbT(self, Z):
        logz13 = np.log(183 * np.power(Z, -1. / 3.))
        eta = np.log(1440 * np.power(Z, -2. / 3.)) / logz13 # (1,A46)
        b = 4. / 3. * (1 + 1 / 9. * ((Z + 1)/(Z + eta)) / logz13) # (1,A45)
        T = self.ta + self.tb # (1,A47)
        epsilon = np.pi * self.me / (2 * self.alpha) * T / ((Z + eta) * logz13) # (1,52)
        return eta, b, T, epsilon

    # R from (1,A83)
    def __fR(self, E, E2, M, sin2th):
        return (M + 2 * E * sin2th) / (M - 2 * E2 * sin2th) # (1,A83)

    # tr from (1,A57)
    def __fVtr(self, VQ2, b, pol=0):
        if pol == 1:
            return 0
        logQ2me = np.log(VQ2 / self.me**2)
        return 1. / b * self.alphapi * (logQ2me - 1) # (1,A57)

    # F from (1,A44)
    def __fVF(self, VQ2, VE, VE2, b, T, spence, pol=0):
        logQ2me = np.log(VQ2 / self.me**2)
        VF = self.__fVF_noQ2(b, T, spence, pol)
        if pol == 0:
            logEE2 = np.log(VE / VE2)**2
            VF = self.__fVF_2(logQ2me, logEE2, VF)
        return VF

    # stripped Q2, E, EF for F
    def __fVF_noQ2(self, b, T, spence, pol):
        VF = 1 + 0.5772 * b * T
        if pol == 0:
            VF = VF + self.alphapi * (np.pi**2 / 6. - spence) - 2 * self.alphapi * 14. / 9. # (1,A44)
        return VF

    # F from (1,A44) with input of 2nd vars
    def __fVF_2(self, logQ2me, logEE2, F_noQ2):
       return F_noQ2 + 2 * self.alphapi * 13. / 12. * logQ2me - self.alphapi / 2. * logEE2

    # phi function (1,A54)
    def __fphiv(self, v):
        return 1 - v + 0.75 * v * v
