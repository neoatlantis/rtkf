#!/usr/bin/env python3

import math
from math import pi, sqrt, sin, cos, tan, asin, acos, atan
from scipy.optimize import newton
import numpy as np
from numpy import array

from ..prop import Properties
from ..prop.units import *
from .choked import getChokingPoint



def maximalMassflow(self, point0_total, A):
    chokingPoint = getChokingPoint(self, point0_total)
    return float(chokingPoint.rho) * A * float(chokingPoint.a)



def estimateMassflow(calculation, speedratio=0.5):
    """Requires a calculation input, with all geometry loaded or calculated."""

    volute, nozzle_row, rotor, diffuser = (
        calculation.volute,
        calculation.nozzle_row,
        calculation.rotor,
        calculation.diffuser
    )

    point0_total = calculation.workingPoint.properties # inlet properties

    alpha3 = nozzle_row.alpha3

    # expand p0 to such a speed that the local speed is 0.5 sound speed

    h0 = float(point0_total.h)
    s = point0_total.s
    point = None
    def iterate(c):
        nonlocal point
        h = h0 - c**2/2
        point = calculation.queryProperties(
            Enthalpy(h).units(),
            s
        )
        return c - float(point.a) * speedratio 
    c = newton(iterate, 0)
    rho = float(point.rho)

    r3, b3 = float(nozzle_row.r3), float(nozzle_row.b3)
    A = 2 * pi * r3 * b3

    m_dot = rho * c * A * sin(alpha3)
    
    return m_dot
