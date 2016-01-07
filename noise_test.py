from __future__ import division

import math

import fice


gnd = fice.Net('gnd')
x = fice.Net('x')

objs = []
objs.append(fice.Ground(gnd))
objs.append(fice.Resistor(100, x, gnd))
objs.append(fice.Resistor(100, x, gnd, temperature=0))

res = fice.do_nodal(objs, 2*math.pi*100)
print res[x.voltage]

res = fice.do_noise(objs, 2*math.pi*100)
print res[x.voltage]

print 4 * fice.k_B * 273 * 50
