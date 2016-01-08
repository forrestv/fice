from __future__ import division

import math

import fice
from fice import s_parameters

Z0 = 50

def f(gnd, p1, p2):
    PAD_LOSS = 10 # dB
    PAD_A = 10**(-PAD_LOSS/20)
    PAD_ZS = Z0
    PAD_RA = PAD_RB = PAD_ZS * (1 - PAD_A) / (1 + PAD_A)
    PAD_RC = (PAD_ZS**2 - PAD_RB**2)/2/PAD_RB
    #print PAD_RA, PAD_RC
    x = fice.Net('pad')
    yield fice.Resistor(PAD_RA, p1, x)
    yield fice.Resistor(PAD_RC, x, gnd)
    yield fice.Resistor(PAD_RB, x, p2)

print s_parameters.power_ratio_to_db(s_parameters.noise(f, None, [Z0, Z0]))
