from __future__ import division

import cmath
import math
import collections

import numpy

k_B = 1.3806485279e-23

class Variable(object):
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return 'Variable(%r)' % (self.name,)

class Net(object):
    def __init__(self, name):
        self.name = name
        self.voltage = Variable((name, 'voltage'))
    def __repr__(self):
        return 'Net(%r)' % (self.name,)

class IndependentComplexRandom(object): pass

class Impedor(object):
    def __init__(self):
        self._current_when_short = Variable((self, '_current_when_short'))
    def get_current_contributions(self, w):
        # return list of (destination node, [(variable, coefficient)], constant)
        # convention is to return currents going into the node
        Z = self.get_impedance(w)
        if Z == 0:
            return [
                (self.net2.voltage, [(self._current_when_short, 1)], 0),
                (self.net1.voltage, [(self._current_when_short, -1)], 0),
                (self._current_when_short, [(self.net2.voltage, 1), (self.net1.voltage, -1)], 0),
            ]
        if Z is None: # inf
            return []
        return [
            (self.net1.voltage, [(self.net1.voltage, -1/Z), (self.net2.voltage, +1/Z)], 0),
            (self.net2.voltage, [(self.net2.voltage, -1/Z), (self.net1.voltage, +1/Z)], 0),
        ]

class Admittor(object):
    def __init__(self, admittance):
        self.get_admittance = admittance
    def get_current_contributions(self, w):
        Y = self.get_impedance(w)
        return [
            (self.net1.voltage, [(self.net1.voltage, -Y), (self.net2.voltage, +Y)], 0),
            (self.net2.voltage, [(self.net2.voltage, -Y), (self.net1.voltage, +Y)], 0),
        ]

class Resistor(Impedor):
    def __init__(self, resistance, net1, net2, temperature=290):
        Impedor.__init__(self)
        self.resistance = resistance
        self.net1 = net1
        self.net2 = net2
        self._nv = IndependentComplexRandom()
        self.temperature = temperature
    def get_impedance(self, w): return self.resistance
    def get_noise_contributions(self, w): # map from variables -> (map from nodes to injected current)
        if self.resistance == 0: return {} # not a hack
        Y = 1/self.resistance
        x = math.sqrt(4 * k_B * self.temperature * Y.real)
        return [
            (self._nv, [(self.net1.voltage, x), (self.net2.voltage, -x)]),
        ]

class Capacitor(Impedor):
    def __init__(self, capacitance, net1, net2):
        Impedor.__init__(self)
        self.capacitance = capacitance
        self.net1 = net1
        self.net2 = net2
    def get_impedance(self, w): return 1/(1j * w * self.capacitance) if w * self.capacitance != 0 else None
    def get_noise_contributions(self, w): return []

class Inductor(Impedor):
    def __init__(self, inductance, net1, net2):
        Impedor.__init__(self)
        self.inductance = inductance
        self.net1 = net1
        self.net2 = net2
    def get_impedance(self, w): return 1j * w * self.inductance
    def get_noise_contributions(self, w): return []

class CurrentSource(object):
    # current source from net1 to net2
    def __init__(self, current, net1, net2):
        self.current = current
        self.net1 = net1
        self.net2 = net2
    def voltage(self, res):
        return res[self.net2.voltage] - res[self.net1.voltage]
    def get_current_contributions(self, w):
        c = self.current(w)
        return [
            (self.net1.voltage, [], -c),
            (self.net2.voltage, [], c),
        ]

class VoltageControlledCurrentSource(object):
    # generates current from net1 to net2 equal to gain * (V(sense1) - V(sense2))
    def __init__(self, gain, net1, net2, sense1, sense2):
        self.gain = gain
        self.net1 = net1
        self.net2 = net2
        self.sense1 = sense1
        self.sense2 = sense2
    def get_current_contributions(self, w):
        g = self.gain(w)
        return [
            (self.net1.voltage, [(self.sense1.voltage, -g), (self.sense2.voltage, g)], 0),
            (self.net2.voltage, [(self.sense1.voltage, g), (self.sense2.voltage, -g)], 0),
        ]
    def get_noise_contributions(self, w): return []

class VoltageSource(object):
    # makes V(net2) - V(net1) = voltage
    def __init__(self, voltage, net1, net2):
        self.voltage = voltage
        self.net1 = net1
        self.net2 = net2
        self._fake_current_var = Variable('_fake_current_var')
    def current(self, res):
        return res[self._fake_current_var]
    def get_current_contributions(self, w):
        return [
            (self.net2.voltage, [(self._fake_current_var, 1)], 0),
            (self.net1.voltage, [(self._fake_current_var, -1)], 0),
            (self._fake_current_var, [(self.net2.voltage, 1), (self.net1.voltage, -1)], -self.voltage(w)),
        ]
    def get_noise_contributions(self, w): return []

def exp(x):
    if isinstance(x, D):
        if isinstance(x._v, complex):
            res = cmath.exp(x._v)
        else:
            res = math.exp(x._v)
        return D(res, scale_dict(x._d, res))
    else:
        if isinstance(x, complex):
            return cmath.exp(x)
        else:
            return math.exp(x)
cosh = lambda x: (exp(x)+exp(-x))/2
sinh = lambda x: (exp(x)-exp(-x))/2
tanh = lambda x: sinh(x)/cosh(x)
coth = lambda x: cosh(x)/sinh(x)
sech = lambda x: 1/cosh(x)
csch = lambda x: 1/sinh(x)

def sqrt(x):
    if isinstance(x, D):
        res = math.sqrt(x._v)
        return D(res, scale_dict(x._d, 1/(2*res)))
    else:
        return math.sqrt(x)

def log10(x):
    if isinstance(x, D):
        return D(math.log10(x._v), scale_dict(x._d, 1/x._v / math.log(10)))
    else:
        return math.log10(x)

def angle(x, deg=0):
    fact = 180/math.pi if deg else 1
    if isinstance(x, D):
        res = cmath.phase(x._v)
        denom = x._v.real**2 + x._v.imag**2
        return D(res, map_values(x._d, lambda val: (-x._v.imag*val.real + x._v.real*val.imag)/denom)) * fact
    else:
        return cmath.phase(x) * fact    

class TransmissionLine(object):
    def __init__(self, characteristic_impedance, attenuation_per_second, length_in_seconds, net1, gnd, net2):
        self.net1 = net1
        self.gnd = gnd
        self.net2 = net2
        self._waver_var = Variable('_waver_var')
        self._wavel_var = Variable('_wavel_var')
        self._Z_0 = characteristic_impedance
        self._attenuation_per_second = attenuation_per_second
        self._length_in_seconds = length_in_seconds
    def get_current_contributions(self, w):
        k = exp(self._attenuation_per_second(w) * self._length_in_seconds) * exp(1j*w*self._length_in_seconds)
        return [
            (self.net1.voltage, [(self._wavel_var, 1), (self._waver_var, -1)], 0),
            (self.net2.voltage, [(self._waver_var, 1/k), (self._wavel_var, -k)], 0),
            (self.gnd.voltage, [(self._wavel_var, -1), (self._waver_var, 1), (self._waver_var, -1/k), (self._wavel_var, k)], 0),
            (self._waver_var, [(self._waver_var, self._Z_0), (self.net1.voltage, -1), (self.gnd.voltage, 1), (self._wavel_var, self._Z_0)], 0),
            (self._wavel_var, [(self._waver_var, self._Z_0/k), (self.net2.voltage, -1), (self.gnd.voltage, 1), (self._wavel_var, self._Z_0*k)], 0),
        ]

def CoupledTransmissionLine(even_impedance, odd_impedance, even_attenuation_per_second, odd_attenuation_per_second, even_length_in_seconds, odd_length_in_seconds):
    '''
    p1 -- p2
    p4 -- p3
    '''
    from fice import s2p
    Z_0 = 50 # sqrt(even_impedance * odd_impedance)
    def get_S(w):
        gamma_l_e = (even_attenuation_per_second(w) + 1j*w) * even_length_in_seconds
        gamma_l_o = ( odd_attenuation_per_second(w) + 1j*w) *  odd_length_in_seconds
        D_e = 2 * even_impedance * Z_0 * cosh(gamma_l_e) + (even_impedance**2 + Z_0**2) * sinh(gamma_l_e)
        D_o = 2 *  odd_impedance * Z_0 * cosh(gamma_l_o) + ( odd_impedance**2 + Z_0**2) * sinh(gamma_l_o)
        X_e = (even_impedance**2 - Z_0**2) * sinh(gamma_l_e) / (2 * D_e)
        X_o = ( odd_impedance**2 - Z_0**2) * sinh(gamma_l_o) / (2 * D_o)
        Y_e = even_impedance * Z_0 / D_e
        Y_o =  odd_impedance * Z_0 / D_o
        S = numpy.array([
            [X_e+X_o, Y_e+Y_o, Y_e-Y_o, X_e-X_o],
            [Y_e+Y_o, X_e+X_o, X_e-X_o, Y_e-Y_o],
            [Y_e-Y_o, X_e-X_o, X_e+X_o, Y_e+Y_o],
            [X_e-X_o, Y_e-Y_o, Y_e+Y_o, X_e+X_o],
        ])
        return S
    def get_model(vs): # list of (vp, vn)
        assert len(vs) == 4
        return s2p._y_box(s2p.memoize(lambda w: s2p._S_to_Y(get_S(w), Z_0)), vs)
    return get_model

def MultiTransmissionLine(L_per_meter, C_per_meter, length_in_meters):
    N = L_per_meter.shape[0]
    assert L_per_meter.shape == (N, N)
    assert C_per_meter.shape == (N, N)
    
    YZ = C_per_meter.dot(L_per_meter) # real is times -w^2
    T = numpy.linalg.eig(YZ)[1] # T_I in text
    Tsum = numpy.sum(T, axis=0)
    Tinv = numpy.linalg.inv(T)
    gamma = numpy.sqrt(numpy.diagonal(Tinv.dot(YZ).dot(T))) # real is times jw
    assert numpy.allclose(YZ, T.dot(numpy.diag(gamma**2)).dot(Tinv))
    Zc = L_per_meter.dot(T).dot(numpy.diag(1/gamma)).dot(Tinv) # 7.34d
    assert numpy.allclose(Zc, numpy.linalg.inv(C_per_meter).dot(T).dot(numpy.diag(gamma)).dot(Tinv))
    Zc_dot_T = Zc.dot(T)
    class MultiTransmissionLineInstance(object):
        def __init__(self, left_ref_net, left_nets, right_ref_net, right_nets):
            assert len(left_nets) == N
            assert len(right_nets) == N
            self._left_ref_net = left_ref_net.voltage
            self._left_nets = [x.voltage for x in left_nets]
            self._right_ref_net = right_ref_net.voltage
            self._right_nets = [x.voltage for x in right_nets]
            self._waver_vars = [Variable('_waver_var%i'%i) for i in xrange(N)]
            self._wavel_vars = [Variable('_wavel_var%i'%i) for i in xrange(N)]
        def get_current_contributions(self, w):
            # k = exp(gamma * z)
            k = numpy.exp(gamma * (1j * w * length_in_meters))
            
            # waver = I+
            # wavel = I-
            
            # left going currents; remember to negate I from that given in 7.86b
            for i, Trow in enumerate(T):
                yield self._left_nets[i], zip(self._waver_vars, -Trow) + zip(self._wavel_vars, Trow), 0
            yield self._left_ref_net, zip(self._waver_vars, Tsum) + zip(self._wavel_vars, -Tsum), 0
            
            # right going currents
            for i, Trow in enumerate(T):
                yield self._right_nets[i], zip(self._waver_vars, Trow/k) + zip(self._wavel_vars, -Trow*k), 0
            yield self._right_ref_net, zip(self._waver_vars, -Tsum/k) + zip(self._wavel_vars, Tsum*k), 0
            
            for i, (var, net, Zc_dot_T_row) in enumerate(zip(self._wavel_vars, self._left_nets, Zc_dot_T)):
                yield var, zip(self._waver_vars, Zc_dot_T_row) + zip(self._wavel_vars, Zc_dot_T_row) + [(net, -1), (self._left_ref_net, 1)], 0
            
            for i, (var, net, Zc_dot_T_row) in enumerate(zip(self._waver_vars, self._right_nets, Zc_dot_T)):
                yield var, zip(self._waver_vars, Zc_dot_T_row/k) + zip(self._wavel_vars, Zc_dot_T_row*k) + [(net, -1), (self._right_ref_net, 1)], 0
    return MultiTransmissionLineInstance

class Ground(object):
    # a 1 ohm resistor to 0 volt reference
    # DO NOT USE MULTIPLE IN ONE CIRCUIT
    def __init__(self, net): self.net = net
    def get_current_contributions(self, w):
        return [
            (self.net.voltage, [(self.net.voltage, -1)], 0),
        ]
    def get_noise_contributions(self, w): return []

def add_dicts(a, b, add_func=lambda a, b: a + b):
    res = dict(a)
    for k, v in b.iteritems():
        res[k] = add_func(res[k], v) if k in res else v
    return res
def map_values(d, f):
    return {k: f(v) for k, v in d.iteritems()}
def scale_dict(d, k):
    return map_values(d, lambda x: k*x)

def inv(A):
    A_v = [[D.v(entry) for entry in row] for row in A]
    A_d = collections.defaultdict(lambda: numpy.zeros(A.shape, dtype=A.dtype))
    for i, row in enumerate(A):
        for j, entry in enumerate(row):
            for k, v in D.d(entry).iteritems():
                A_d[k][i, j] = v
    
    A_inv = numpy.linalg.inv(A_v)
    A_inv2 = A_inv.copy()
    for k in A_d:
        A_inv2 = A_inv2 + [[D(0, {k: entry}) for entry in row] for row in -A_inv.dot(A_d[k]).dot(A_inv)]
    return A_inv2

def do_nodal(objects, w=0):
    gnd, = [obj.net for obj in objects if isinstance(obj, Ground)]
    
    equations = {}
    var_list = set()
    for obj in objects:
        for dest, coeff, const in obj.get_current_contributions(w):
            assert isinstance(dest, Variable), dest
            x = equations.get(dest, ({}, 0))
            equations[dest] = x[0], x[1] + const
            for k, v in coeff:
                x[0][k] = x[0].get(k, 0) + v
    
    equations = list(equations.iteritems())
    net_list = [k for k, v in equations]
    
    doing_differentiation = any(isinstance(const, D) or any(isinstance(x, D) for x in coeff.itervalues()) for dest, (coeff, const) in equations)
    if doing_differentiation:
        A = numpy.zeros((len(equations), len(equations)), dtype=complex)
        b = numpy.zeros((len(equations)), dtype=complex)
        A_d = collections.defaultdict(lambda: numpy.zeros((len(equations), len(equations)), dtype=complex))
        b_d = collections.defaultdict(lambda: numpy.zeros((len(equations)), dtype=complex))
        for dest, (coeff, const) in equations:
            for k, v in coeff.iteritems():
                A[net_list.index(dest), net_list.index(k)] = D.v(v)
                for k2, v2 in D.d(v).iteritems():
                    A_d[k2][net_list.index(dest), net_list.index(k)] = v2
            b[net_list.index(dest)] = D.v(-const)
            for k2, v2 in D.d(-const).iteritems():
                b_d[k2][net_list.index(dest)] = v2
        A_inv = numpy.linalg.inv(A)
        x = list(A_inv.dot(b))
        for k in set(A_d)|set(b_d):
            A_inv_d_b = -A_inv.dot(A_d[k].dot(A_inv.dot(b)))
            x_d = A_inv_d_b + A_inv.dot(b_d[k])
            x = [a+D(0, {k:v}) for a, v in zip(x, x_d)]
    else:
        A = numpy.zeros((len(equations), len(equations)), dtype=complex)
        b = numpy.zeros((len(equations)), dtype=complex)
        for dest, (coeff, const) in equations:
            for k, v in coeff.iteritems():
                A[net_list.index(dest), net_list.index(k)] = v
            b[net_list.index(dest)] = -const
        x = numpy.linalg.solve(A, b)
    
    res = dict(zip(net_list, x))
    assert (res[gnd.voltage]*res[gnd.voltage].conjugate()).real < 1e-12
    return res

class D(object):
    def __init__(self, value, dvalue):
        # these asserts are a bit slow
        assert isinstance(value, (int, float, complex, long))
        assert isinstance(dvalue, dict)
        assert all(isinstance(k, Variable) for k in dvalue)
        assert all(isinstance(v, (int, float, complex, long)) for v in dvalue.itervalues())
        self._v = value
        self._d = dvalue
    def __neg__(self):
        return -1*self
    def __add__(self, other):
        if not isinstance(other, D): other = D(other, {})
        return D(self._v+other._v, add_dicts(self._d, other._d))
    def __radd__(self, other):
        return self+other
    def __sub__(self, other):
        return self + (-1*other)
    def __rsub__(self, other):
        return (-1*self) + other
    def __mul__(self, other):
        if not isinstance(other, D): other = D(other, {})
        return D(self._v*other._v, add_dicts(scale_dict(self._d, other._v), scale_dict(other._d, self._v)))
    def __rmul__(self, other):
        return self*other
    def __truediv__(self, other):
        if not isinstance(other, D): other = D(other, {})
        return D(self._v/other._v, add_dicts(scale_dict(self._d, 1/other._v), scale_dict(other._d, -self._v/other._v**2)))
    def __rtruediv__(self, other):
        if not isinstance(other, D): other = D(other, {})
        return other/self
    def __pow__(self, other):
        if not isinstance(other, D): other = D(other, {})
        res = self._v ** other._v
        res_d = scale_dict(self._d, self._v ** (other._v-1) * other._v)
        if other._d:
            res_d = add_dicts(res_d, scale_dict(other._d, math.log(self._v)*res))
        return D(res, res_d)
    def __rpow__(self, other):
        if not isinstance(other, D): other = D(other, {})
        return other ** self
    def __mod__(self, other):
        return self.__class__(self._v % other, self._d)
    def __lt__(self, other): return self._v < other
    def __le__(self, other): return self._v <= other
    def __eq__(self, other): return self._v == other
    def __ne__(self, other): return self._v != other
    def __ge__(self, other): return self._v >= other
    def __gt__(self, other): return self._v > other
    @classmethod
    def v(cls, x):
        return (cls(0, {})+x)._v
    @classmethod
    def d(cls, x):
        return (cls(0, {})+x)._d
    @property
    def real(self):
        return self.__class__(self._v.real, map_values(self._d, lambda x: x.real))
    @property
    def imag(self):
        return self.__class__(self._v.imag, map_values(self._d, lambda x: x.imag))
    def conjugate(self):
        return self.__class__(self._v.conjugate(), map_values(self._d, lambda x: x.conjugate()))    
    def __abs__(self):
        return sqrt((self * self.conjugate()).real)
    def __repr__(self):
        return self.__class__.__name__ + repr((self._v, self._d))

def do_noise(objects, w=0): # returns map var -> variance(var)
    noise_vars = set()
    
    for obj in objects:
        noise_vars.update(nvar for nvar, dests in obj.get_noise_contributions(w))
    
    res = {}
    
    for noise_var in noise_vars:
        class Surrogate(object):
            def __init__(self, inner):
                self._inner = inner
            def get_current_contributions(self, w2):
                assert w2 is None
                
                x = [(dest, coeff, 0) for dest, coeff, const in self._inner.get_current_contributions(w)]
                for nvar, dests in self._inner.get_noise_contributions(w):
                    assert isinstance(dests, list)
                    if nvar is noise_var:
                        for node, coeff in dests:
                            x.append((node, [], coeff))
                return x
        objects2 = map(lambda o: Surrogate(o) if not isinstance(o, Ground) else o, objects)
        res2 = do_nodal(objects2, None)
        for k, v in res2.iteritems():
            res[k] = res.get(k, 0) + (v.real*v.real + v.imag*v.imag)
    
    return res
