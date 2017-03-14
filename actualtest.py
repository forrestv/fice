from __future__ import division

import math

import numpy
from matplotlib import pyplot

import fice
from fice import s_parameters
from autoee import util

c = 299792458

@apply
def test_transline():
    Z0 = 50
    LENGTH_IN_SECONDS = 1e-6

    @util.chain(util.flatten)
    def f1(gnd, p1, p2):
        yield fice.TransmissionLine(Z0, lambda w: 0, LENGTH_IN_SECONDS, p1, gnd, p2)

    @util.chain(util.flatten)
    def f2(gnd, p1, p2):
        C = numpy.array([[1/c/Z0]])
        L = numpy.array([[Z0/c]])
        yield fice.MultiTransmissionLine(L, C, LENGTH_IN_SECONDS * c)(gnd, [p1], gnd, [p2])

    for f in [f1, f2]:
        ws = 2*math.pi*numpy.linspace(0, 1e9, 1000)

        for w in ws:
            S = s_parameters.ana(f, w, [Z0]*2)
            x = numpy.exp(-LENGTH_IN_SECONDS * 1j * w)
            assert numpy.allclose(S, numpy.array([[0, x], [x, 0]]))

@apply
def test_coupled_transline():
    EVEN_IMPEDANCE = 120.71
    ODD_IMPEDANCE = 20.71
    LENGTH_IN_SECONDS = 1/10e9/4
    
    @util.chain(util.flatten)
    def f1(gnd, p1, p2, p3, p4):
        yield fice.CoupledTransmissionLine(
            even_impedance=EVEN_IMPEDANCE,
            odd_impedance=ODD_IMPEDANCE,
            even_attenuation_per_second=lambda w: 0,
            odd_attenuation_per_second=lambda w: 0,
            even_length_in_seconds=LENGTH_IN_SECONDS,
            odd_length_in_seconds=LENGTH_IN_SECONDS,
        )([
            (gnd, p1),
            (gnd, p2),
            (gnd, p3),
            (gnd, p4),
        ])
    
    @util.chain(util.flatten)
    def f2(gnd, p1, p2, p3, p4):
        ze, zo = EVEN_IMPEDANCE, ODD_IMPEDANCE
        c11, c12 = (zo + ze)/(2*c*zo*ze), (zo - ze)/(2*c*zo*ze)
        C = numpy.array([[c11, c12], [c12, c11]])
        L = numpy.linalg.inv(C)/c**2
        yield fice.MultiTransmissionLine(L, C, LENGTH_IN_SECONDS * c)(gnd, [p1, p4], gnd, [p2, p3])

    ws = 2*math.pi*numpy.linspace(1e9, 19e9, 1000)

    for w in ws:
        assert numpy.allclose(
            s_parameters.ana(f1, w, [50]*4),
            s_parameters.ana(f2, w, [50]*4),
        )
