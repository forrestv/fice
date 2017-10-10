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
    return sy.dot(numpy.eye(N) - s).dot(fice.inv(numpy.eye(N) + s)).dot(sy)

def _Y_to_S(y, R):
    N = y.shape[0]
    assert y.shape == (N, N)
    sz = numpy.eye(N) * math.sqrt(R)
    return (numpy.eye(N) - sz.dot(y).dot(sz)).dot(numpy.linalg.inv(numpy.eye(N) + sz.dot(y).dot(sz)))

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
            lines = f.readlines()
        lines = [line[:line.index('!')] if '!' in line else line for line in lines]
        lines = map(str.strip, lines)
        lines = [line for line in lines if line]
        
        line = lines[0]
        assert line[0] == '#'; line = line[1:].strip()
        freq_units_str, type_str, format_str, R_str, Rn_str = line.split()
        freq_multiplier = dict(GHZ=1e9, MHZ=1e6, KHZ=1e3, HZ=1)[freq_units_str.upper()]
        assert type_str == 'S'
        parse_func = {
            'DB': lambda m_db, a: cmath.rect(10**(m_db/20), math.radians(a)),
            'MA': lambda m, a: cmath.rect(m, math.radians(a)),
            'RI': lambda r, i: complex(r, i),
        }[format_str]
        assert R_str == 'R'
        ref_impedance = float(Rn_str)
        
        lines = [map(float, line.split()) for line in lines[1:]]
        
        if len(lines[0]) == 3: ports = 1
        elif len(lines[0]) == 7: ports = 3
        elif len(lines[0]) == 9:
            if len(lines[1]) == 9:
                ports = 2
            else:
                total = 4
                for line in lines[1:]:
                    if len(line) % 2: break
                    total += len(line) // 2
                ports = int(round(math.sqrt(total)))
                assert total == ports**2
        else: assert False
        
        res = {}
        noise = {}
        if ports == 2:
            while lines:
                if len(lines[0]) != 9: break
                x = map(parse_func, lines[0][1::2], lines[0][2::2])
                res[lines[0][0] * freq_multiplier] = numpy.array([
                    [x[0], x[2]],
                    [x[1], x[3]],
                ])
                lines.pop(0)
            while lines:
                freq, fmin, Gm, Ga, rn = lines[0]
                noise[freq * freq_multiplier] = 10**(fmin/10), cmath.rect(Gm, math.radians(Ga)), rn*ref_impedance
                lines.pop(0)
        else:
            while lines:
                freq = lines[0][0]
                acc = [map(parse_func, lines[0][1::2], lines[0][2::2])]
                lines.pop(0)
                while lines and len(lines[0]) % 2 == 0:
                    acc.extend([map(parse_func, lines[0][::2], lines[0][1::2])])
                    lines.pop(0)
                res[freq * freq_multiplier] = numpy.array(acc).reshape((ports, ports))
        
        return cls(ports, res, noise, ref_impedance)

    def __init__(self, ports, freq_to_s, freq_to_noise, ref_impedance):
        self._ports = ports
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
        if self._freqs[0] * (1 - 1e-6) <= freq < self._freqs[0] * (1 + 1e-6):
            freq = self._freqs[0] * (1 + 1e-6)
        elif self._freqs[-1] * (1 - 1e-6) < freq <= self._freqs[-1] * (1 + 1e-6):
            freq = self._freqs[-1] * (1 - 1e-6)
        i = bisect.bisect_right(self._freqs, freq)
        assert i >= 1
        fl, fr = self._freqs[i-1], self._freqs[i]
        assert fl <= freq <= fr, (fl, freq, fr)
        x = (freq - fl)/(fr - fl)
        return self._freq_to_s[fl] * (1-x) + self._freq_to_s[fr] * x
    
    def get_noise(self, freq):
        if self._noise_freqs[0] * (1 - 1e-6) <= freq < self._noise_freqs[0] * (1 + 1e-6):
            freq = self._noise_freqs[0] * (1 + 1e-6)
        elif self._noise_freqs[-1] * (1 - 1e-6) < freq <= self._noise_freqs[-1] * (1 + 1e-6):
            freq = self._noise_freqs[-1] * (1 - 1e-6)
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
        min_R_n = (F_o-1)/4/G_o # otherwise, G_u is negative
        if R_n < min_R_n:
            print 'warning: pushing impossible noise parameters into valid region (R_n from %f to %f)' % (R_n, min_R_n)
            R_n = min_R_n
        
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
        if e2m - abs(c)**2 < 0: # should never be violated with perfect arithmetic due to above correction
            d = 0
        else:
            d = math.sqrt(e2m - abs(c)**2)
        
        return dict(a=a, c=c, d=d) # i = a n1, e = c n1 + d n2
    
    def get_model(self, vs): # list of (vp, vn)
        assert len(vs) == self._ports
        return _y_box(memoize(lambda w: _S_to_Y(self.get_S(w/2/math.pi), self._ref_impedance)), vs)
    
    def get_noisy_model(self, vs):
        assert len(vs) == self._ports == 2
        
        n1 = fice.IndependentComplexRandom()
        n2 = fice.IndependentComplexRandom()
        
        np = memoize(lambda w: self.get_noise2(w/2/math.pi))
        
        tmp = fice.Net('tmp')
        res = []
        res.extend(self.get_model([(tmp, vs[0][1]), (vs[1][0], vs[1][1])]))
        res.append(CurrentNoiseSource(lambda w: [(n1, np(w)['a'])                  ], tmp, vs[0][1]))
        res.append(VoltageNoiseSource(lambda w: [(n1, np(w)['c']), (n2, np(w)['d'])], tmp, vs[0][0]))
        return res

class CurrentNoiseSource(object):
    # current goes from net1 to net2
    def __init__(self, current, net1, net2):
        self.current = current
        self.net1 = net1
        self.net2 = net2
    def get_current_contributions(self, w): return []
    def get_noise_contributions(self, w):
        return [
            (nv, [(self.net1.voltage, coeff), (self.net2.voltage, -coeff)])
        for nv, coeff in self.current(w)]

class VoltageNoiseSource(object):
    # makes V(net2) - V(net1) = voltage
    def __init__(self, voltage, net1, net2):
        self.voltage = voltage
        self.net1 = net1
        self.net2 = net2
        self._fake_current_var = fice.Variable('_fake_current_var')
    def get_current_contributions(self, w):
        return [
            (self.net2.voltage, [(self._fake_current_var, 1)], 0),
            (self.net1.voltage, [(self._fake_current_var, -1)], 0),
            (self._fake_current_var, [(self.net2.voltage, 1), (self.net1.voltage, -1)], 0),
        ]
    def get_noise_contributions(self, w):
        return [
            (nv, [(self._fake_current_var, --coeff)])
        for nv, coeff in self.voltage(w)]
