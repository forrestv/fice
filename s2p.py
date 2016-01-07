from __future__ import division

import math
import cmath
import bisect

import numpy

import fice


def _S_to_Y(s, R):
    N = s.shape[0]
    assert s.shape == (N, N)
    sy = numpy.eye(N) * math.sqrt(1/R)
    return sy.dot(numpy.eye(N) - s).dot(numpy.linalg.inv(numpy.eye(N) + s)).dot(sy)

def _y_box(Y, vs): # list of (vp, vn)
    #assert Y(0).shape == (len(vs), len(vs))
    for i in xrange(len(vs)):
        for j in xrange(len(vs)):
            yield fice.VoltageControlledCurrentSource(lambda w, i=i, j=j: Y(w)[i,j], vs[i][0], vs[i][1], vs[j][0], vs[j][1])

def memoize(f):
    val = []
    def _(*args):
        if not val or val[0] != args:
            res = f(*args)
            val[:] = [args, res]
        return val[1]
    return _

class S2P(object):
    @classmethod
    def from_file(cls, filename):
        # implemented according to http://cp.literature.agilent.com/litweb/pdf/genesys200801/sim/linear_sim/sparams/touchstone_file_format.htm
        with open(filename) as f:
            line = f.readline().strip()
            
            while not line or line[0] == '!': line = f.readline().strip() # skip comments
            
            assert line[0] == '#'; line = line[1:].strip()
            freq_units_str, type_str, format_str, R_str, Rn_str = map(str.strip, line.split())
            freq_multiplier = dict(GHZ=1e9, MHZ=1e6, KHZ=1e3, HZ=1)[freq_units_str.upper()]
            assert type_str == 'S'
            assert format_str in ['MA', 'RI']
            assert R_str == 'R'
            ref_impedance = float(Rn_str)
            
            line = f.readline().strip()
            while not line or line[0] == '!': line = f.readline().strip() # skip comments
            
            res = {}
            while True:
                if format_str == 'MA':
                    freq, s11m, s11a, s21m, s21a, s12m, s12a, s22m, s22a = map(float, line.strip().split())
                    res[freq * freq_multiplier] = numpy.matrix([
                        [cmath.rect(s11m, math.radians(s11a)), cmath.rect(s12m, math.radians(s12a))],
                        [cmath.rect(s21m, math.radians(s21a)), cmath.rect(s22m, math.radians(s22a))],
                    ])
                elif format_str == 'RI':
                    freq, s11r, s11i, s21r, s21i, s12r, s12i, s22r, s22i = map(float, line.strip().split())
                    res[freq * freq_multiplier] = numpy.matrix([
                        [complex(s11r, s11i), complex(s12r, s12i)],
                        [complex(s21r, s21i), complex(s22r, s22i)],
                    ])
                else: assert False
                line = f.readline().strip()
                if not line or line[0] == '!': break
            return cls(res, ref_impedance)

    def __init__(self, freq_to_s, ref_impedance):
        self._freq_to_s = freq_to_s
        self._freqs = sorted(list(freq_to_s))
        self._ref_impedance = ref_impedance
    
    def get_frequencies(self):
        return self._freqs
    
    def get_closest_frequency(self, f):
        # XXX actually returns freq that is <= f
        return self._freqs[bisect.bisect_left(self._freqs, f)]
    
    def get_S(self, freq):
        i = bisect.bisect_right(self._freqs, freq)
        assert i >= 1
        fl, fr = self._freqs[i-1], self._freqs[i]
        assert fl <= freq <= fr, (fl, freq, fr)
        x = (freq - fl)/(fr - fl)
        return self._freq_to_s[fl] * (1-x) + self._freq_to_s[fr] * x
    
    def get_model(self, vs): # list of (vp, vn)
        assert len(vs) == 2
        return _y_box(memoize(lambda w: _S_to_Y(self.get_S(w/2/math.pi), self._ref_impedance)), vs)
