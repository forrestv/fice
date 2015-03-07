from __future__ import division

import fice


H, W = 2, 5
nets = {(i, j): fice.Net('%s,%s' % (i, j)) for i in xrange(H) for j in xrange(W)}
objs = []
for i in xrange(H-1):
    for j in xrange(W):
        objs.append(fice.Resistor(1, nets[i, j], nets[i+1, j]))
for i in xrange(H):
    for j in xrange(W-1):
        objs.append(fice.Resistor(1, nets[i, j], nets[i, j+1]))
objs.append(fice.Ground(nets[0, 0]))
objs.append(fice.VoltageSource(lambda w: 1, nets[0, 0], nets[H-1, W-1]))
res = fice.do_nodal(objs)
for i in xrange(H):
    for j in xrange(W):
        print '%.02f' % (res[nets[i, j].voltage].real,),
    print
