from __future__ import division

import math

import numpy

import fice

def power_ratio_to_db(x):
    if x == 0:
        return -1e30000
    else:
        return 10*math.log10(x)

def ana(f, w, r=[50, 50]): # r is list of reference impendances per port
    N = len(r)
    
    ret = numpy.zeros((N, N), dtype=complex)
    
    for drive_index in xrange(N):
        objects = []
        gnd = fice.Net('gnd'); objects.append(fice.Ground(gnd))
        DRIVE_VOLTAGE = 1
        v = fice.Net('v'); objects.append(fice.VoltageSource(lambda w: DRIVE_VOLTAGE, gnd, v))
        ports = [fice.Net('port%i' % i) for i in xrange(N)]
        for i, (port, R) in enumerate(zip(ports, r)):
            objects.append(fice.Resistor(R, v if i == drive_index else gnd, port))
        objects.extend(f(gnd, *ports))
        res = fice.do_nodal(objects, w)
        for sense_index in xrange(N):
            ret[sense_index][drive_index] = 2 * (res[ports[sense_index].voltage]/math.sqrt(r[sense_index])) / (DRIVE_VOLTAGE/math.sqrt(r[drive_index])) - (1 if sense_index == drive_index else 0)
    
    return ret

def print_S(S):
    for i, row in enumerate(S):
        for j, entry in enumerate(row):
            print 'S%i_%i' % (i, j), entry, '=', abs(entry)**2, '=', power_ratio_to_db(abs(entry)**2), 'dB'
