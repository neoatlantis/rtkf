#!/usr/bin/env python3

import math
from math import pi, sqrt, sin, cos, tan, asin, acos, atan
from numpy import array

from ._report import CalculationReport
from ..prop import Properties
from ..prop.units import *


average = lambda l: sum(l) / len(l)



def findY(self,
    c_theta_u, r_u, b_u,    # speed and geometry of upstream
    r1, b1, r3, b3,         # geometry of this passage
    alpha1,                 # flow angle of incoming stream
    p1_t, p1, p3_t, p3
):

    # TODO: boundary layer!
    Y_p = 0 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    A_u = 2 * pi * r_u * b_u
    A1 = 2 * pi * r1 * b1

    Y_in = ((A1/A_u-1)*sin(alpha1)) ** 2 * (p1_t-p1) / (p3_t-p3)        #(9-82)

    return Y_p + Y_in                                                   #(9-83)



def solve(self,
    c_m_u, c_theta_u,   # upstream speed
    r_u, b_u,           # upstream outlet geometry
    r1, b1, r3, b3,     # geometry for this passage
    
    p_ut, h_ut,         # upstream total pressure and enthalpy
    m_dot               # mass flow
):
    medium = self.parent.workingPoint.medium

    c_theta = lambda r: c_theta_u * r_u / r                             #(9-81)

    h1_t = h_ut
    rho1 = rho_u
    A_u = 2 * pi * r_u * b_u
    A1 = 2 * pi * r1 * b1
    A3 = 2 * pi * r3 * b3

    c_m1 = A_u / A1 * c_m_u # speed at station 1, as uncompressible flow
    c_theta1 = c_theta(r1)
    c_theta3 = c_theta(r3)

    h1 = h1_t - 0.5 * (c_m1**2 + c_theta1**2)
    point1 = Properties(medium).query(
        Enthalpy(h1).units(),
        Density(rho1).units()
    )
    point1_t = Properties(medium).query(
        Enthalpy(h1_t).units(),
        point1.s
    )

    alpha1 = atan(c_m1 / c_theta1)
    p1_t = float(point1_t.p)
    p1 = float(point1.p)
    point3, point3_t = None

    def iteration(p3, c3m):
        nonlocal point3, point3_t
        h3_t = h1_t
        h3 = h3_t - 0.5 * (c_theta3**2 + c3m**2)
        point3 = Properties(medium).query(
            Enthalpy(h3).units(),
            Pressure(p3).units()
        )
        point3_t = Properties(medium).query(
            point3.s,
            Enthalpy(h3_t).units()
        )
        p3_t = float(point3_t.h)
        
        Y1 = (p_ut - p3_t) / (p3_t - p3)
        Y2 = findY(self, 
            c_theta_u=c_theta_u, r_u=r_u, b_u=b_u,
            r1=r1, b1=b1, r3=r3, b3=b3,
            alpha1=alpha1,
            p1_t=p1_t, p1=p1, p3_t=p3_t, p3=p3
        )

        e1 = Y1 - Y2
        e2 = m_dot - float(point3.rho)*A3*c3m
        return array([e1, e2])











class VanelessPassageCalculation(CalculationReport):
    pass 