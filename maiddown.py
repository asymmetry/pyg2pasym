#!/usr/bin/env python

import os, re, sys, urllib2
import numpy as np
from os.path import exists, join
from multiprocessing.dummy import Pool as tp

try:
    from tools import zload, zdump, Q2range, Wrange
except:
    from pyg2pasym import zload, zdump, Q2range, Wrange

class maiddown:
    def __init__(self):
        self.chans = [0, 1, 2, 3]
        self.Nchan = len(self.chans)
        self.theta = 6
        self.datatype = '2007tot'
        self.datadir = 'data'
        self.tmppkl = join(self.datadir, 'tmp/{chan}_{Q2}.pkl')
        self.allQ2, self.NQ2 = Q2range()
        self.allW, self.NW, self.Wmin, self.Wmax, self.Wstep = Wrange()
	self.Ndown = {}

    def down(self):
        url = self.geturl(self.datatype)
        if not exists(join(self.datadir, 'tmp')):
            os.makedirs(join(self.datadir, 'tmp'))

        def pooldown(par):
            k = str(par[0]) + str(par[1])
            if not self.Ndown.has_key(k):
                self.Ndown[k] = 0
            else:
                self.Ndown[k] += 1
            print self.Ndown[k], par[0], par[1]
            furl = url.format(CHANNEL = par[1] + 1, Q2 = par[0], WMIN = self.Wmin, WMAX = self.Wmax, STEP = self.Wstep, THETA = self.theta)
            try:
                web = urllib2.urlopen(furl)
                content = web.read()
            except Exception as err:
                print err
                return False
            zdump(content, self.tmppkl.format(chan = par[1], Q2 = par[0]))

        pars = []
        for Q2 in self.allQ2:
            for c in self.chans:
                fn = self.tmppkl.format(chan = c, Q2 = Q2)
                if exists(fn):
                    continue
                pars.append([Q2, c])
        if len(pars) < 1:
            return True
        else:
            print pars
            pool = tp(processes = 10)
            pool.map(pooldown, pars)
            return False

    def save(self):
        DMAID1D = {}
        #DMAID1D['info'] = ''
        keys = []
        ndata = 0
        while not self.down():
            pass
        for Q2 in self.allQ2:
            print Q2
            for cn in self.chans:
                fn = self.tmppkl.format(chan = cn, Q2 = Q2)
                content = zload(fn)
                for l in re.split('\n', content):
                    if len(l) < 1 or l[0] == '<':
                        continue
                    c = [x for x in re.split('\s', l) if len(x) > 0]
                    try:
                        c = [float(x) for x in c]
                        for i in range(len(c)):
                            DMAID1D[keys[i]][ndata] = c[i]
                        DMAID1D['Q2'][ndata] = Q2
                        DMAID1D['chan'][ndata] = cn
                        ndata += 1
                    except Exception as err:
                        #if Q2 == self.Q2min and chan == 0:
                        #    DMAID1D['info'] += '%s\n'%l
                        if c[0] == 'W':
                            keys = c
                            if not DMAID1D.has_key('Q2'):
                                DMAID1D['Q2'] = np.zeros(self.Nchan * self.NQ2 * self.NW, dtype = 'float32')
                                DMAID1D['chan'] = np.zeros(self.Nchan * self.NQ2 * self.NW, dtype = 'float32')
                            for i in range(len(c)):
                                if not DMAID1D.has_key(c[i]):
                                    DMAID1D[c[i]] = np.zeros(self.Nchan * self.NQ2 * self.NW, dtype = 'float32')

        zdump(DMAID1D, join(self.datadir, 'MAID1D{}.pdt'.format(self.datatype)))

    def geturl(self, c):
        url = {}
        url['2007tot'] = 'http://portal.kph.uni-mainz.de/cgi-bin/maid1?switch=217&param2={CHANNEL}&param50=2&value35={Q2}&value36={WMIN}&value41={STEP}&value42={WMAX}&param99=0&param11=1&param12=1&param13=1&param14=1&param15=1&param16=1&param17=1&param18=1&param19=1&param20=1&param21=1&param22=1&param23=1&param24=1&param25=1&param26=1&value11=1.0&value12=1.0&value13=1.0&value51=1.0&value52=1.0&value53=1.0&value54=1.0&value55=1.0&value56=1.0&value57=1.0&value58=1.0&value59=1.0&value60=1.0&value61=1.0&value62=1.0&value63=1.0&value64=1.0&value65=1.0&value66=1.0&value67=1.0&value68=1.0&value69=1.0&value70=1.0&value71=1.0&value72=1.0&value73=1.0&value74=1.0&value75=1.0&value76=1.0&value77=1.0&value78=1.0&value79=1.0&value80=1.0&value81=1.0&value82=1.0&value83=1.0&value84=1.0'
        url['2007diff'] = 'http://portal.kph.uni-mainz.de/cgi-bin/maid1?switch=211&param2={CHANNEL}&param50=2&value35={Q2}&value36={WMIN}&value37={THETA}&value41={STEP}&value42={WMAX}&param99=0&param11=1&param12=1&param13=1&param14=1&param15=1&param16=1&param17=1&param18=1&param19=1&param20=1&param21=1&param22=1&param23=1&param24=1&param25=1&param26=1&value11=1.0&value12=1.0&value13=1.0&value51=1.0&value52=1.0&value53=1.0&value54=1.0&value55=1.0&value56=1.0&value57=1.0&value58=1.0&value59=1.0&value60=1.0&value61=1.0&value62=1.0&value63=1.0&value64=1.0&value65=1.0&value66=1.0&value67=1.0&value68=1.0&value69=1.0&value70=1.0&value71=1.0&value72=1.0&value73=1.0&value74=1.0&value75=1.0&value76=1.0&value77=1.0&value78=1.0&value79=1.0&value80=1.0&value81=1.0&value82=1.0&value83=1.0&value84=1.0'
        url['chiraltot'] = 'http://portal.kph.uni-mainz.de/cgi-bin/maid1?switch=517&param2={CHANNEL}&param16=2&param11={Q2}&param12={WMIN}&param23={STEP}&param25={WMAX}&value9=-1.216&value8=-1.092&value20=4.337&value21=-4.260&valuee48=5.235&valuee49=0.925&valuee50=2.205&valuee51=6.629&valuee52=-4.103&valuee53=-2.654&valuee112=9.342&valuee67=-8.269&valuee68=-0.925&valuee69=-1.035&valuee71=-4.352&valuee72=10.593&valuee73=2.120&valuee113=-13.745&valuee70=3.910'
        url['chiraldiff'] = 'http://portal.kph.uni-mainz.de/cgi-bin/maid1?switch=511&param2={CHANNEL}&param50=2&param35={Q2}&param36={WMIN}&param37={THETA}&param41={STEP}&param42={WMAX}&value9=-1.216&value8=-1.092&value20=4.337&value21=-4.260&valuee48=5.235&valuee49=0.925&valuee50=2.205&valuee51=6.629&valuee52=-4.103&valuee53=-2.654&valuee112=9.342&valuee67=-8.269&valuee68=-0.925&valuee69=-1.035&valuee71=-4.352&valuee72=10.539&valuee73=2.120&valuee113=-13.745&valuee70=3.910'
        return url[c]

if __name__ == '__main__':
    d = maiddown()
    #d.datatype = 'chiraltot'
    d.save()
