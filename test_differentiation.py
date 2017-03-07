from __future__ import division

import math

import fice


v1 = fice.Net('v1')
v2 = fice.Net('v2')
v22 = fice.Net('v22')
v3 = fice.Net('v3')
gnd = fice.Net('gnd')

theta = fice.Variable('theta')
theta2 = fice.Variable('theta2')

objs = [
    fice.VoltageSource(lambda w: 1, gnd, v1),
    fice.Resistor(0, v1, v2),
    fice.TransmissionLine(fice.D(50, {theta: 1}), lambda w: 0, 0.2, v2, gnd, v22),
    fice.TransmissionLine(fice.D(50, {theta2: 1}), lambda w: 0, 0.3, v22, gnd, v3),
    fice.Resistor(50, v3, gnd),
    #fice.VoltageSource(lambda w: 0, gnd, v3),
    fice.Ground(gnd),
]

res = fice.do_nodal(objs, w=1*2*math.pi)
objs[2]._Z_0 = 50
objs[3]._Z_0 = 50
res2 = fice.do_nodal(objs, w=1*2*math.pi)
objs[2]._Z_0 = 50.1
res3 = fice.do_nodal(objs, w=1*2*math.pi)
objs[2]._Z_0 = 50
objs[3]._Z_0 = 50.1
res4 = fice.do_nodal(objs, w=1*2*math.pi)

for k in res:
    a, b = fice.D.d(res[k]).get(theta, 0), (res3[k]-res2[k])/0.1
    print '%.5f %.5f %.5f %.5f' % (a.real, a.imag, b.real, b.imag), k
print
for k in res:
    a, b = fice.D.d(res[k])[theta2], (res4[k]-res2[k])/0.1
    print '%.5f %.5f %.5f %.5f' % (a.real, a.imag, b.real, b.imag), k
