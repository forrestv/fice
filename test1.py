from __future__ import division

import fice


v1 = fice.Net('v1')
v2 = fice.Net('v2')
gnd = fice.Net('gnd')
res = fice.do_nodal([
    fice.Ground(gnd),
    fice.VoltageSource(lambda w: 1, gnd, v1),
    fice.Resistor(10, v1, v2),
    fice.Capacitor(10, v2, gnd),
], 10)
for var, value in res.iteritems():
    print var.name, '=', value
