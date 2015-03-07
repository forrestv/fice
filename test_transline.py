from __future__ import division

import math

import fice


v1 = fice.Net('v1')
v2 = fice.Net('v2')
v22 = fice.Net('v22')
v3 = fice.Net('v3')
gnd = fice.Net('gnd')

objs = [
    fice.VoltageSource(lambda w: 1, gnd, v1),
    fice.Resistor(0, v1, v2),
    fice.TransmissionLine(50, lambda w: 0, 0.2, v2, gnd, v22),
    fice.TransmissionLine(50, lambda w: 0, 0.3, v22, gnd, v3),
    fice.Resistor(50, v3, gnd),
    #fice.VoltageSource(lambda w: 0, gnd, v3),
    fice.Ground(gnd),
]

res = fice.do_nodal(objs, w=1*2*math.pi)

print res[gnd.voltage]
print res[v1.voltage]
print res[v2.voltage]
print res[v3.voltage]
print 'current', objs[0].current(res)
print 'waver', res[objs[3]._waver_var]
print 'wavel', res[objs[3]._wavel_var]
