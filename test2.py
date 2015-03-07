from __future__ import division

import math

import fice


def f(x):
    v1 = fice.Net('v1')
    v2 = fice.Net('v2')
    gnd = fice.Net('gnd')
    res = fice.do_nodal([
        fice.Ground(gnd),
        fice.VoltageSource(lambda w: 1, gnd, v1),
        fice.Resistor(10, v1, v2),
        fice.Capacitor(.002, v2, gnd),
    ], x)
    return abs(res[v2.voltage])

for x in xrange(100):
    xx = 2**(x/12)
    print str(xx).rjust(30), 20*math.log10(f(xx))
