#!/usr/bin/env python3

import math
from math import pi, sqrt, sin, cos, tan, asin, acos, atan, copysign
from scipy.optimize import newton, fminbound, fsolve
from scipy import integrate
from numpy import array

from ._report import CalculationReport
from ..prop import Properties
from ..prop.units import *

average = lambda l: sum(l) / len(l)
averageR = lambda r1, r2: sqrt(0.5 * (r1**2 + r2**2))
deg2rad = lambda deg: deg / 180 * pi
rad2deg = lambda rad: rad / pi * 180



def findDelta2(self):
    # ---- INPUTS
    p1_total, p1 = self.p1_total, self.p1
    p2_ideal = self.p2_ideal
    b1, r1, z1, A1 = self.b1, self.r1, self.z1, self.A1
    r2, z2, A2 = self.r2, self.z2, self.A2
    phi_1, phi_2 = 0, 0 #!!!!!! for CFTurbo cyclindal diffuser only
    # -------------------------------------------------------------------------

    self.comment("Find delta2 - fractional blockage")

    d = sqrt( (z2 - z1) ** 2 + (r2 - r1) ** 2)
    delta_phi = abs(phi_2 - phi_1)
    
    if delta_phi > 0:
        L = d * delta_phi / ( 2 * sin(delta_phi / 2) )      #(9-88)
    else:
        L = d

    theta_c = atan(b1 * (A2 / A1 - 1) / (2 * L))

    p_vr = (p1_total - p1) / (p1_total - p2_ideal)          #(9-90)
    D = (sqrt(p_vr) + 1) ** 2 / 4                           #(9-92)

    self.set("p_vr", p_vr, "difuser velocity-pressure ratio")
    self.set("D", D, "diffusion factor")

    K_theta = 2 * theta_c / deg2rad(11)
    if K_theta < 1: K_theta = 1                             #(9-93)

    K1 = 0.005 + (K_theta) / 5                              #(9-94)
    K2 = 2 * theta_c * (\
        1 - 2*rad2deg(theta_c)/(22*K_theta)
    ) / (125*K_theta)                                       #(9-95)

    delta_2 = (K1 + K2 * (D-1)) * L * A1 / (A2 * b1)        #(9-91)

    return delta_2 # Diffuser fractional blockage, used for mass balance





class DiffuserCalculation(CalculationReport):

    RELATIVE_ERROR = 1e-6

    def __init__(self, calc):
        self.parent = calc
        self.queryProperties = self.parent.queryProperties
        self.config = calc.config.diffuser # input configuration

        CalculationReport.__init__(self) # recording begins
        try:
            self.override(self.config)
        except:
            self.__findGeometry()
        self.__deriveGeometry()

    def __findPreviousElementOutput(self):
        """Copy calculation results from upstreaming element."""
        self.comment("Copy calculation results from previous element.")
        upstreamingElement = self.parent.rotor
        mapping = {
            # last element => current element
            "m_dot"         : "m_dot",
            "p3_total"      : "p1_total",
            "p3"            : "p1",
            "rho3"          : "rho1",
            "h3"            : "h1",
            "h3_total"      : "h1_total",
            "s3"            : "s1",
            "c3"            : "c1",
            "c3_m"          : "c1_m",
            "c3_theta"      : "c1_theta",
        }
        for lastEntry in sorted(mapping.keys()):
            self.set(
                mapping[lastEntry], # key assigned: new name
                getattr(upstreamingElement, lastEntry)) # value: from old name

    def skip(self):
        """Skip calculation, copy rotor output as final."""
        upstreamingElement = self.parent.rotor
        mapping = {
            # last element => current element
            "m_dot"         : "m_dot",
            "p3_total"      : "p2_total",
            "p3"            : "p2",
            "rho3"          : "rho2",
            "h3"            : "h2",
            "h3_total"      : "h2_total",
            "s3"            : "s2",
            "c3"            : "c2",
            "c3_m"          : "c2_m",
            "c3_theta"      : "c2_theta",
        }
        for lastEntry in sorted(mapping.keys()):
            self.set(
                mapping[lastEntry], # key assigned: new name
                getattr(upstreamingElement, lastEntry)) # value: from old name

    def __findGeometry(self):
        MainDimensionsElement = self.config.MainDimensionsElement
        MerInlet = MainDimensionsElement.MerData["MerInlet"]
        MerOutlet = MainDimensionsElement.MerData["MerOutlet"]

        z1hub = MerInlet.Hub.x + MerInlet.OffsetHub.x
        z1shr = MerInlet.Shr.x + MerInlet.OffsetShr.x
        self.set("z1", average([z1hub, z1shr]))

        z2hub = MerOutlet.Hub.x
        z2shr = MerOutlet.Shr.x
        self.set("z2", average([z2hub, z2shr]))

        r1hub = MerInlet.Hub.y + MerInlet.OffsetHub.y
        r1shr = MerInlet.Shr.y + MerInlet.OffsetShr.y
        self.set("r1", averageR(*[r1hub, r1shr]))

        r2hub = MerOutlet.Hub.y
        r2shr = MerOutlet.Shr.y
        self.set("r2", averageR(*[r2hub, r2shr]))

        self.set("b1", r1shr - r1hub)
        self.set("b2", r2shr - r2hub)



    def __deriveGeometry(self):
        self.set("A1", 2*pi*self.r1*self.b1)
        self.set("A2", 2*pi*self.r2*self.b2)


    def __p2_ideal(self):
        """Find the ideal discharge static pressure"""
        c2_theta = self.c1_theta * self.r1 / self.r2
        self.set("c2_theta", c2_theta)

        point2_ideal = None
        def iterate(c2_m_ideal):
            nonlocal point2_ideal
            h2_total = float(self.h1_total)
            h2 = h2_total - 0.5 * (c2_m_ideal**2 + c2_theta**2)
            s2 = self.s1
            point2_ideal = self.queryProperties(
                Enthalpy(h2).units(),
                Entropy(s2).units()
            )
            rho2 = float(point2_ideal.rho)
            return c2_m_ideal - (self.m_dot / (rho2 * self.A2))
        c2_m_ideal = newton(iterate, 0)
        self.set("p2_ideal", point2_ideal.p, "ideal discharge static pressure")
        self.set("c2_m_ideal", c2_m_ideal)
        return point2_ideal


    def __p2(self):
        point2 = None
        def iterate(c2_m):
            nonlocal point2
            h2_total = float(self.h1_total)
            h2 = h2_total - 0.5 * (c2_m**2 + float(self.c2_theta)**2)
            s2 = self.s1
            point2 = self.queryProperties(
                Enthalpy(h2).units(),
                Entropy(s2).units()
            )
            rho2 = float(point2.rho)
            A2_blocked = self.A2 * (1-self.delta2)
            return c2_m - (self.m_dot / (rho2 * A2_blocked))
        c2_m = newton(iterate, 0)
        self.set("p2", point2.p, "actual discharge static pressure")
        return point2


    def solve(self):
        self.__findPreviousElementOutput()

        point2_ideal = self.__p2_ideal()

        delta2 = findDelta2(self)
        self.set("delta2", delta2, "fractional discharge blockage")

        point2_withdelta = self.__p2()


        c2_m = 0
        c2_m_last = 0
        while True:
            c2 = sqrt(c2_m**2 + self.c2_theta**2)
            h2 = float(self.h1_total) - c2**2/2
            point2 = self.queryProperties(
                Enthalpy(h2).units(),
                point2_withdelta.p
            )
            rho2 = float(point2.rho)
            c2_m_last = c2_m
            c2_m = self.m_dot / (rho2 * self.A2)

            if abs(c2_m_last/c2_m - 1) < self.RELATIVE_ERROR: break



        self.set("c2_m", c2_m)
        self.set("c2", c2)
        self.set("h2", point2.h)
        self.set("p2", point2.p)
        self.set("s2", point2.s)
        self.set("T2", point2.T)

        point2_total = self.queryProperties(
            self.h1_total,
            point2.s
        )

        self.set("h2_total", point2_total.h)
        self.set("p2_total", point2_total.p)
        self.set("T2_total", point2_total.T)


        print("Diffuser p2", self.p2)