import fice

v1 = fice.Net('v1')
v2 = fice.Net('v2')
gnd = fice.Net('gnd')

objs = [
    fice.VoltageSource(lambda w: 1, gnd, v1),
    fice.CurrentSource(lambda w: 5, v1, v2),
    fice.Resistor(1, v2, gnd),
    fice.Resistor(1, v1, gnd),
    fice.Ground(gnd),
]

res = fice.do_nodal(objs)

print res[v1.voltage]
print res[v2.voltage]
print objs[0].current(res)
print objs[1].voltage(res)
