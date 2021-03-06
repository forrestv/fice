from __future__ import division

import math

import numpy

import fice
import s2p
from autoee import util

def power_ratio_to_db(x):
    if x == 0:
        return -1e30000
    else:
        return 10*math.log10(x)

def _rig(f, r, drive_voltages, r_temperatures):
    N = len(r)
    assert len(drive_voltages) == N
    
    objects = []
    gnd = fice.Net('gnd'); objects.append(fice.Ground(gnd))
    vs = [fice.Net('v%i'%i) for i in xrange(N)]
    for i in xrange(N):
        objects.append(fice.VoltageSource(lambda w, i=i: drive_voltages[i], gnd, vs[i]))
    ports = [fice.Net('port%i' % i) for i in xrange(N)]
    for i, (port, R) in enumerate(zip(ports, r)):
        objects.append(fice.Resistor(R, vs[i], port, temperature=r_temperatures[i]))
    objects.extend(f(gnd, *ports))
    return objects, ports

def ana(f, w, r=[50, 50], drives=None, return_res=False, dtype=complex): # r is list of reference impendances per port
    N = len(r)
    
    ret = numpy.full((N, N), numpy.nan, dtype=dtype)
    ress = []

    for drive_index in (xrange(N) if drives is None else drives):
        DRIVE_VOLTAGE = 2
        objects, ports = _rig(f, r, [DRIVE_VOLTAGE if i == drive_index else 0 for i in xrange(N)], [290 for i in xrange(N)])
        res = fice.do_nodal(objects, w)
        for sense_index in xrange(N):
            ret[sense_index][drive_index] = 2 * (res[ports[sense_index].voltage]/fice.sqrt(r[sense_index])) / (DRIVE_VOLTAGE/fice.sqrt(r[drive_index])) - (1 if sense_index == drive_index else 0)
        ress.append(res)
    
    return (ret, ress) if return_res else ret

def noise(f, w, r):
    assert len(r) == 2

    @util.chain(util.flatten)
    def modified_f(gnd, p1, p2, temp=290):
        tmp = fice.Net('tmp')
        tmp2 = fice.Net('tmp2')
        yield f(gnd, p1, tmp)
        Y = s2p._S_to_Y(numpy.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]), r[1])
        yield s2p._y_box(lambda w: Y, [(tmp, gnd), (p2, gnd), (tmp2, gnd)])
        yield fice.Resistor(r[1], tmp2, gnd, temperature=temp)

    objects, ports = _rig(lambda *args: modified_f(*args, temp=290), r, [0, 0], [290, 0])
    output_noise = fice.do_noise(objects, w)[ports[1].voltage]
    
    objects, ports = _rig(lambda *args: modified_f(*args, temp=290), r, [0, 0], [0, 0])
    output_noise_added = fice.do_noise(objects, w)[ports[1].voltage]
    
    res = output_noise / (output_noise - output_noise_added)
    #print power_ratio_to_db(res)
    
    return res
    

def print_S(S):
    for i, row in enumerate(S):
        for j, entry in enumerate(row):
            print 'S%i_%i' % (i, j), entry, '=', abs(entry)**2, '=', power_ratio_to_db(abs(entry)**2), 'dB'
