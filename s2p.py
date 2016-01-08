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
            
            line = f.readline().strip()
            while not line or line[0] == '!': line = f.readline().strip() # skip comments
            
            noise = {}
            while True:
                freq, fmin, Gm, Ga, rn = map(float, line.split())
                freq *= freq_multiplier
                #print freq, fmin, Gm, Ga, rn
                noise[freq] = 10**(fmin/10), cmath.rect(Gm, math.radians(Ga)), rn*ref_impedance
                #print freq, noise[freq]
                line = f.readline().strip()
                if not line or line[0] == '!': break
            
            return cls(res, noise, ref_impedance)

    def __init__(self, freq_to_s, freq_to_noise, ref_impedance):
        self._freq_to_s = freq_to_s
        self._freqs = sorted(list(freq_to_s))
        self._freq_to_noise = freq_to_noise
        self._noise_freqs = sorted(list(freq_to_noise))
        self._ref_impedance = ref_impedance
    
    def get_frequencies(self):
        return self._freqs
    
    def get_noise_frequencies(self):
        return self._noise_freqs
    
    def get_S(self, freq):
        i = bisect.bisect_right(self._freqs, freq)
        assert i >= 1
        fl, fr = self._freqs[i-1], self._freqs[i]
        assert fl <= freq <= fr, (fl, freq, fr)
        x = (freq - fl)/(fr - fl)
        return self._freq_to_s[fl] * (1-x) + self._freq_to_s[fr] * x
    
    def get_noise(self, freq):
        i = bisect.bisect_right(self._noise_freqs, freq)
        assert i >= 1
        fl, fr = self._noise_freqs[i-1], self._noise_freqs[i]
        assert fl <= freq <= fr, (fl, freq, fr)
        x = (freq - fl)/(fr - fl)
        return (
            self._freq_to_noise[fl][0] * (1-x) + self._freq_to_noise[fr][0] * x,
            self._freq_to_noise[fl][1] * (1-x) + self._freq_to_noise[fr][1] * x,
            self._freq_to_noise[fl][2] * (1-x) + self._freq_to_noise[fr][2] * x,
        )
    
    def get_noise2(self, freq):
        F_o, Gamma, R_n = self.get_noise(freq)
        Z = self._ref_impedance * (1 + Gamma)/(1 - Gamma)
        Y = 1/Z
        G_o, B_o = Y.real, Y.imag
        
        B_gamma = -B_o
        G_gamma = (F_o - 1)/2/R_n - G_o
        G_u = G_o**2*R_n - R_n * G_gamma**2
        Y_gamma = complex(G_gamma, B_gamma)
        T_0 = 290
        
        i2m = 4*fice.k_B*T_0*(abs(Y_gamma)**2 * R_n + G_u)
        e2m = 4*fice.k_B*T_0*R_n
        eiconjm = Y_gamma.conjugate() * e2m
        
        ieconjm = eiconjm.conjugate()
        a = math.sqrt(i2m)
        c = (ieconjm/a).conjugate()
        d = math.sqrt(e2m - abs(c)**2)
        
        return dict(a=a, c=c, d=d) # i = a n1, e = c n1 + d n2
    
    def get_model(self, vs): # list of (vp, vn)
        assert len(vs) == 2
        return _y_box(memoize(lambda w: _S_to_Y(self.get_S(w/2/math.pi), self._ref_impedance)), vs)
    
    def get_noisy_model(self, vs):
        assert len(vs) == 2
        
        n1 = fice.IndependentComplexRandom()
        n2 = fice.IndependentComplexRandom()
        
        np = memoize(lambda w: self.get_noise2(w/2/math.pi))
        
        tmp = fice.Net('tmp')
        res = []
        res.extend(self.get_model([(tmp, vs[0][1]), (vs[1][0], vs[1][1])]))
        res.append(CurrentNoiseSource(lambda w: {n1: np(w)['a']                }, tmp, vs[0][1]))
        res.append(VoltageNoiseSource(lambda w: {n1: np(w)['c'], n2: np(w)['d']}, tmp, vs[0][0]))
        return res

class CurrentNoiseSource(object):
    # current goes from net1 to net2
    def __init__(self, current, net1, net2):
        self.current = current
        self.net1 = net1
        self.net2 = net2
    def get_current_contributions(self, w): return {}
    def get_noise_contributions(self, w):
        return {
            nv: {self.net1.voltage: coeff, self.net2.voltage: -coeff}
        for nv, coeff in self.current(w).iteritems()}

class VoltageNoiseSource(object):
    # makes V(net2) - V(net1) = voltage
    def __init__(self, voltage, net1, net2):
        self.voltage = voltage
        self.net1 = net1
        self.net2 = net2
        self._fake_current_var = fice.Variable('_fake_current_var')
    def get_current_contributions(self, w):
        return {
            self.net2.voltage: ({self._fake_current_var: 1}, 0),
            self.net1.voltage: ({self._fake_current_var: -1}, 0),
            self._fake_current_var: ({self.net2.voltage: 1, self.net1.voltage: -1}, 0),
        }
    def get_noise_contributions(self, w):
        return {
            nv: {self._fake_current_var: -coeff}
        for nv, coeff in self.voltage(w).iteritems()}
