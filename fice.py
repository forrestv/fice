from __future__ import division

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
    def get_current_contributions(self, frequency):
        # return list of (destination node, {net voltage coefficients}, constant)
        Z = self.get_impedance(frequency)
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
    def get_other_constraints(self):
        return []

class Resistor(Impedor):
    def __init__(self, resistance, net1, net2):
        Impedor.__init__(self)
        self.resistance = resistance
        self.net1 = net1
        self.net2 = net2
    def get_impedance(self, frequency): return self.resistance

class Capacitor(Impedor):
    def __init__(self, capacitance, net1, net2):
        Impedor.__init__(self)
        self.capacitance = capacitance
        self.net1 = net1
        self.net2 = net2
    def get_impedance(self, frequency): return 1/(1j * frequency * self.capacitance) if frequency else None

class Inductor(Impedor):
    def __init__(self, inductance, net1, net2):
        Impedor.__init__(self)
        self.inductance = inductance
        self.net1 = net1
        self.net2 = net2
    def get_impedance(self, frequency): return 1j * frequency * self.inductance

class CurrentSource(object):
    def __init__(self, current, net1, net2):
        self.current = current
        self.net1 = net1
        self.net2 = net2
    def voltage(self, res):
        return res[self.net2.voltage] - res[self.net1.voltage]
    def get_current_contributions(self, frequency):
        # return list of (destination node, {net voltage coefficients}, constant)
        return {
            self.net1.voltage: ({}, -self.current),
            self.net2.voltage: ({}, self.current),
        }

class VoltageSource(object):
    def __init__(self, voltage, net1, net2):
        self.voltage = voltage
        self.net1 = net1
        self.net2 = net2
        self._fake_current_var = Variable('_fake_current_var')
    def current(self, res):
        return res[self._fake_current_var]
    def get_current_contributions(self, frequency):
        return {
            self.net2.voltage: ({self._fake_current_var: 1}, 0),
            self.net1.voltage: ({self._fake_current_var: -1}, 0),
            self._fake_current_var: ({self.net2.voltage: 1, self.net1.voltage: -1}, -self.voltage),
        }

class Ground(object):
    # a 1 ohm resistor to 0 volt reference
    def __init__(self, net): self.net = net
    def get_current_contributions(self, frequency):
        return {
            self.net.voltage: ({self.net.voltage: -1}, 0),
        }

def add_dicts(a, b, add_func=lambda a, b: a + b):
    res = dict(a)
    for k, v in b.iteritems():
        res[k] = add_func(res[k], v) if k in res else v
    return res

def do_nodal(objects, freq=0):
    equations = {}
    var_list = set()
    for obj in objects:
        for dest, (coeff, const) in obj.get_current_contributions(freq).iteritems():
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

'''v1 = Net('v1')
v2 = Net('v2')
gnd = Net('gnd')
res = do_nodal([
    Ground(gnd),
    VoltageSource(1, gnd, v1),
    Resistor(10, v1, v2),
    Capacitor(10, v2, gnd),
], 10)
for net, v in res.iteritems():
    print net.name, v'''

if 0:
    def f(x):
        v1 = Net('v1')
        v2 = Net('v2')
        gnd = Net('gnd')
        res = do_nodal([
            Ground(gnd),
            VoltageSource(1, gnd, v1),
            Resistor(10, v1, v2),
            Capacitor(.002, v2, gnd),
        ], x)
        return abs(res[v2])
    for x in xrange(100):
        xx = 2**(x/12)
        print str(xx).rjust(30), 20*math.log10(f(xx))

if __name__ == '__main__':
    H, W = 2, 5
    nets = {(i, j): Net('%s,%s' % (i, j)) for i in xrange(H) for j in xrange(W)}
    objs = []
    for i in xrange(H-1):
        for j in xrange(W):
            objs.append(Resistor(1, nets[i, j], nets[i+1, j]))
    for i in xrange(H):
        for j in xrange(W-1):
            objs.append(Resistor(1, nets[i, j], nets[i, j+1]))
    objs.append(Ground(nets[0, 0]))
    objs.append(VoltageSource(1, nets[0, 0], nets[H-1, W-1]))
    res = do_nodal(objs)
    for i in xrange(H):
        for j in xrange(W):
            print '%.02f' % (res[nets[i, j]].real,),
        print
