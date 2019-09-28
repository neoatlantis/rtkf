#!/usr/bin/env python3
from scipy.optimize import newton
from ..prop import Properties
from ..prop.units import *


class ChokedException(Exception):

    def __init__(self, massflow, p_max=None, p_min=None, where=""):
        Exception.__init__(
            self,
             "Calculation found one component choked.\n" + 
             (" Massflow=%f" % massflow) +
             ((" Pressure %f <= p <= %f" % (p_min, p_max)) if (p_min) else "") +
             (" at %s" % where)
        )

        self.where = where
        self.m_dot = massflow
        self.p_min = p_min
        self.p_max = p_max



def getChokingPoint(self, point_total):
    """Expand the point with total conditions isentropic, and get another point
    that yields c=a.
      This condition yields a point with pressure and velocity, giving a limit
    on the minimum of achievable discharge pressure in nozzle row and rotor.
    """
    h0 = float(point_total.h)
    s = point_total.s
    point = None
    def iterate(c):
        nonlocal point
        h = h0 - c**2/2
        point = self.queryProperties(
            Enthalpy(h).units(),
            s
        )
        return c - float(point.a)
    c = newton(iterate, 0)
    return point