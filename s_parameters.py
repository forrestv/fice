from __future__ import division

import math

import numpy

import fice

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

def ana(f, w, r=[50, 50]): # r is list of reference impendances per port
    N = len(r)
    
    ret = numpy.zeros((N, N), dtype=complex)
    
    for drive_index in xrange(N):
        DRIVE_VOLTAGE = 1
        objects, ports = _rig(f, r, [DRIVE_VOLTAGE if i == drive_index else 0 for i in xrange(N)], [290 for i in xrange(N)])
        res = fice.do_nodal(objects, w)
        for sense_index in xrange(N):
            ret[sense_index][drive_index] = 2 * (res[ports[sense_index].voltage]/math.sqrt(r[sense_index])) / (DRIVE_VOLTAGE/math.sqrt(r[drive_index])) - (1 if sense_index == drive_index else 0)
    
    return ret

def noise(f, w, r):
    assert len(r) == 2
    
    objects, ports = _rig(f, r, [0, 0], [290, 0])
    output_noise = fice.do_noise(objects, w)[ports[1].voltage]
    
    objects, ports = _rig(f, r, [0, 0], [0, 0])
    output_noise_added = fice.do_noise(objects, w)[ports[1].voltage]
    
    res = output_noise / (output_noise - output_noise_added)
    #print power_ratio_to_db(res)
    
    return res
    

def print_S(S):
    for i, row in enumerate(S):
        for j, entry in enumerate(row):
            print 'S%i_%i' % (i, j), entry, '=', abs(entry)**2, '=', power_ratio_to_db(abs(entry)**2), 'dB'
