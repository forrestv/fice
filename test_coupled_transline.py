from __future__ import division

import math

import numpy
from matplotlib import pyplot

import fice
from fice import s_parameters
from autoee import util

@util.chain(util.flatten)
def f(gnd, p1, p2, p3, p4):
    yield fice.CoupledTransmissionLine(
        even_impedance=120.71,
        odd_impedance=20.71,
        even_attenuation_per_second=lambda w: 0,
        odd_attenuation_per_second=lambda w: 0,
        even_length_in_seconds=1/10e9/4,
        odd_length_in_seconds=1/10e9/4,
    )([
        (gnd, p1),
        (gnd, p2),
        (gnd, p3),
        (gnd, p4),
    ])

c = 299792458

@util.chain(util.flatten)
def f2(gnd, p1, p2, p3, p4):
    ze, zo = 120.71, 20.71
    c11, c12 = (zo + ze)/(2*c*zo*ze), (zo - ze)/(2*c*zo*ze)
    C = numpy.array([[c11, c12], [c12, c11]])
    L = numpy.linalg.inv(C)/c**2
    yield fice.MultiTransmissionLine(L, C, 1/10e9/4 * c)(gnd, [p1, p4], gnd, [p2, p3])

ws = 2*math.pi*numpy.linspace(1e9, 19e9, 1000)

for w in ws:
    assert numpy.allclose(
        s_parameters.ana(f, w, [50]*4),
        s_parameters.ana(f2, w, [50]*4),
    )

if 0:
    Ses = [s_parameters.ana(f, w, [50]*4) for w in ws]
    
    for i in xrange(4):
        for j in xrange(4):
            pyplot.plot(ws/(2*math.pi), [10*math.log10(abs(S[i, j])**2) for S in Ses])
    pyplot.show()
