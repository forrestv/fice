from __future__ import division

import math
import cmath
import sys

import numpy

import fice
from fice import s2p, s_parameters

d1 = s2p.S2P.from_file('BGB707L7ESD_3V06M0.s2p')
for freq in d1.get_noise_frequencies():
    F_o, Gamma, R_n = d1.get_noise(freq)
    Z = d1._ref_impedance * (1 + Gamma)/(1 - Gamma)
    def f(gnd, port1, port2):
        x = fice.Net('x')
        yield fice.Resistor(1j*Z.imag, x, port1)
        for x in d1.get_noisy_model([(x, gnd), (port2, gnd)]): yield x
    rs = [Z.real, 50]
    #S = s_parameters.ana(f, 2*math.pi*freq, rs)
    N = s_parameters.noise(f, 2*math.pi*freq, rs)
    assert abs(N - F_o) / F_o < 1e-6
    print freq/1e9, N, F_o
