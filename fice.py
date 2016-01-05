from __future__ import division

import cmath
import math

import numpy


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

class Impedor(object):
    def __init__(self):
        self._current_when_short = Variable((self, '_current_when_short'))
    def get_current_contributions(self, w):
        # return list of (destination node, {net voltage coefficients}, constant)
        # convention is to return currents going into the node
        Z = self.get_impedance(w)
        if Z == 0:
            return {
                self.net2.voltage: ({self._current_when_short: 1}, 0),
                self.net1.voltage: ({self._current_when_short: -1}, 0),
                self._current_when_short: ({self.net2.voltage: 1, self.net1.voltage: -1}, 0),
            }
        if Z is None: # inf
            return {}
        return {
            self.net1.voltage: ({self.net1.voltage: -1/Z, self.net2.voltage: +1/Z}, 0),
            self.net2.voltage: ({self.net2.voltage: -1/Z, self.net1.voltage: +1/Z}, 0),
        }

class Resistor(Impedor):
    def __init__(self, resistance, net1, net2):
        Impedor.__init__(self)
        self.resistance = resistance
        self.net1 = net1
        self.net2 = net2
    def get_impedance(self, w): return self.resistance

class Capacitor(Impedor):
    def __init__(self, capacitance, net1, net2):
        Impedor.__init__(self)
        self.capacitance = capacitance
        self.net1 = net1
        self.net2 = net2
    def get_impedance(self, w): return 1/(1j * w * self.capacitance) if w else None

class Inductor(Impedor):
    def __init__(self, inductance, net1, net2):
        Impedor.__init__(self)
        self.inductance = inductance
        self.net1 = net1
        self.net2 = net2
    def get_impedance(self, w): return 1j * w * self.inductance

class CurrentSource(object):
    # current source from net1 to net2
    def __init__(self, current, net1, net2):
        self.current = current
        self.net1 = net1
        self.net2 = net2
    def voltage(self, res):
        return res[self.net2.voltage] - res[self.net1.voltage]
    def get_current_contributions(self, w):
        # return list of (destination node, {net voltage coefficients}, constant)
        c = self.current(w)
        return {
            self.net1.voltage: ({}, -c),
            self.net2.voltage: ({}, c),
        }

class VoltageControlledCurrentSource(object):
    # generates current from net1 to net2 equal to gain * (V(sense1) - V(sense2))
    def __init__(self, gain, net1, net2, sense1, sense2):
        self.gain = gain
        self.net1 = net1
        self.net2 = net2
        self.sense1 = sense1
        self.sense2 = sense2
    def get_current_contributions(self, w):
        # return list of (destination node, {net voltage coefficients}, constant)
        return {
            self.net1.voltage: ({self.sense1.voltage: -self.gain, self.sense2.voltage: self.gain}, 0),
            self.net2.voltage: ({self.sense1.voltage: self.gain, self.sense2.voltage: -self.gain}, 0),
        }

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
        return {
            self.net2.voltage: ({self._fake_current_var: 1}, 0),
            self.net1.voltage: ({self._fake_current_var: -1}, 0),
            self._fake_current_var: ({self.net2.voltage: 1, self.net1.voltage: -1}, -self.voltage(w)),
        }

class TransmissionLine(object):
    def __init__(self, characteristic_impedance, attenuation_per_second, length_in_seconds, net1, gnd, net2):
        self.net1 = net1
        self.gnd = gnd
        self.net2 = net2
        self._waver_var = Variable('_wave1_var')
        self._wavel_var = Variable('_wave2_var')
        self._Z_0 = characteristic_impedance
        self._attenuation_per_second = attenuation_per_second
        self._length_in_seconds = length_in_seconds
    def get_current_contributions(self, w):
        k = math.exp(-self._attenuation_per_second(w) * self._length_in_seconds) * cmath.exp(-1j*w*self._length_in_seconds)
        return {
            self.net1.voltage: ({self._wavel_var: 2/self._Z_0, self.net1.voltage: -1/self._Z_0, self.gnd.voltage: 1/self._Z_0}, 0),
            self.net2.voltage: ({self._waver_var: 2/self._Z_0, self.net2.voltage: -1/self._Z_0, self.gnd.voltage: 1/self._Z_0}, 0),
            self.gnd.voltage: ({self._wavel_var: -2/self._Z_0, self._waver_var: -2/self._Z_0, self.net1.voltage: 1/self._Z_0, self.net2.voltage: 1/self._Z_0, self.gnd.voltage: -2/self._Z_0}, 0),
            self._waver_var: ({self._waver_var: -1, self.net1.voltage: k, self.gnd.voltage: -k, self._wavel_var: -k}, 0),
            self._wavel_var: ({self._wavel_var: -1, self.net2.voltage: k, self.gnd.voltage: -k, self._waver_var: -k}, 0),
        }

class Ground(object):
    # a 1 ohm resistor to 0 volt reference
    # DO NOT USE MULTIPLE IN ONE CIRCUIT
    def __init__(self, net): self.net = net
    def get_current_contributions(self, w):
        return {
            self.net.voltage: ({self.net.voltage: -1}, 0),
        }

def add_dicts(a, b, add_func=lambda a, b: a + b):
    res = dict(a)
    for k, v in b.iteritems():
        res[k] = add_func(res[k], v) if k in res else v
    return res

def do_nodal(objects, w=0):
    equations = {}
    var_list = set()
    for obj in objects:
        for dest, (coeff, const) in obj.get_current_contributions(w).iteritems():
            x = equations.get(dest, ({}, 0))
            equations[dest] = add_dicts(x[0], coeff), x[1] + const
    
    equations = list(equations.iteritems())
    net_list = [k for k, v in equations]
    A = numpy.zeros((len(equations), len(equations)), dtype=complex)
    b = numpy.zeros((len(equations), 1), dtype=complex)
    for dest, (coeff, const) in equations:
        for k, v in coeff.iteritems():
            A[net_list.index(dest), net_list.index(k)] = v
        b[net_list.index(dest)] = -const
    x = numpy.linalg.solve(A, b)
    
    return dict(zip(net_list, map(lambda a: complex(a.real, a.imag), x)))
