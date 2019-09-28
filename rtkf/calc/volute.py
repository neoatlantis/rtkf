#!/usr/bin/env python3

import enum
import math
from math import pi, sqrt
from scipy.optimize import newton

from ._report import CalculationReport
from ..prop import Properties
from ..prop.units import *
from .multivar_newton import newton_multi
from .estimate_massflow import maximalMassflow
from .choked import ChokedException
from .boundary_layer import BoundaryLayerCalculation
from logging import *


class DesignRule:

    PFLEIDERER = "constant angular momentum (Pfleiderer)"
    STEPANOFF = "constant angular speed (Stepanoff)"
    RADIUS_CURVE = "curve of radius"
    AREA_CURVE = "curve of area"



class VoluteCalculation(CalculationReport):

    """Calculations for the volute part.

    Input
    :calc - An instance of calc.Calculation, the main calculation object.
    """

    RELATIVE_ERROR = 1e-6

    def __init__(self, calc):
        self.parent = calc
        self.config = calc.config.volute # input configuration
        self.queryProperties = self.parent.queryProperties

        #print(self.config.SpiralCasing)
        #print(self.config.SpiralCasingBC)

        CalculationReport.__init__(self)
        try:
            self.override(self.config)
        except:
            self.__findGeometry()
        self.__deriveGeometry()

    def __findGeometry(self):
        """Find important geometry variables. All variables are notated
        according to H. Aungier's book __Turbine Aerodynamics__."""

        self.comment("Geometry parameters.\n")

        self.comment("Design rule:")
        self.designRule = self.__findDesignRule()

        r3 = (
            self.config.SpiralCasingBC.MerData["MerInlet"].Hub.y +
            self.config.SpiralCasingBC.MerData["MerInlet"].Shr.y
        ) / 2
        A1 = self.__A1()
        r1 = r3 + (A1 / pi) ** 0.5
        
        self.comment("At inlet:")
        self.A1 = A1
        self.set("r1", r1, "approx. radius, assuming circle cross-section!")


        self.comment("At outlet:")
        self.r3 = r3
        self.set("b3", abs(
            self.config.SpiralCasingBC.MerData["MerInlet"].Hub.x -
            self.config.SpiralCasingBC.MerData["MerInlet"].Shr.x
        ), "width")
        


    def __deriveGeometry(self):
        self.A2 = self.A1 / 2
        self.r2 = self.r3 + (self.A2 / pi)**0.5
        self.A3 = 2 * pi * self.r3 * self.b3


    def __findAerodynamicalProperties(self):

        self.comment("\nAerodynamical parameters at inlet.")

        self.set(
            "m_dot_max",
            maximalMassflow(
                self, self.parent.workingPoint.properties, self.A1),
            "maximal massflow that may enter"
        )

        if self.parent.workingPoint.m_dot > self.m_dot_max:
            raise ChokedException(self.m_dot_max, "volute1") 

        self.set("m_dot", self.parent.workingPoint.m_dot)
        self.set("p1_total", self.parent.workingPoint.p_total_inlet, "total")
        self.set("h1_total", self.parent.workingPoint.properties.h)

        self.__findInletStaticProperties()
        self.__findMiddleStaticProperties()


    ###########################################################################

    def __massflowBalanceForStation1And2(self, massflow, area):
        workingPoint = self.parent.workingPoint
        point0_total = workingPoint.properties
        point1 = None

        c1 = 0
        s1 = float(point0_total.s)
        rho1 = float(point0_total.rho)
        p1_total = float(point0_total.p)

        c1_last = 0

        iterations = 0
        while True:
            iterations += 1
            if iterations > 100:
                raise Exception("Massflow balance @ volute station 1 or 2 not convergence.")

            point1_total = self.queryProperties(
                point0_total.h,
                Pressure(p1_total).units()
            )

            c1_last = c1
            c1 = massflow / (float(rho1) * area)
            if abs(c1_last / c1 - 1) < self.RELATIVE_ERROR:
                break

            h1 = float(point0_total.h) - c1**2/2
            point1 = self.queryProperties(
                point1_total.s,
                Enthalpy(h1).units()
            )
            rho1 = float(point1.rho)

        return point1, c1

    def __findInletStaticProperties(self):
        self.comment("Finding static conditions at inlet via iteration.")
        print("*************************************")
        point1, c1 = self.__massflowBalanceForStation1And2(
            massflow = self.m_dot,
            area = float(self.A1)
        )
        self.set("s1", point1.s)
        self.set("c1", c1)
        self.set("h1", point1.h)
        self.set("p1", point1.p)
        self.set("rho1", point1.rho)
        self.set("a1", point1.a, "sound speed at inlet")


    def __findMiddleStaticProperties(self):
        point2, c2 = self.__massflowBalanceForStation1And2(
            massflow = 0.5 * self.m_dot,
            area = float(self.A2)
        )
        self.set("c2", c2)




    ###########################################################################

    def __A1(self):
        SpiralCasing = self.config.SpiralCasing
        if self.designRule in [DesignRule.PFLEIDERER, DesignRule.STEPANOFF]:
            self.comment("Velocity based design rules. Area from volume flow.")
            try:
                assert SpiralCasing.AutoUpdate
            except Exception as e:
                print("Error: design rules not set as auto updated.")
                raise e
            return SpiralCasing.Qi / SpiralCasing.cu
        else:
            lastValue = SpiralCasing.BezierProgr.Points[-1].y
            meaning = SpiralCasing.ProgressionType
            if self.designRule == DesignRule.AREA_CURVE:
                self.comment("Geometry based design rules. Area given.")
                return lastValue
            elif self.designRule == DesignRule.RADIUS_CURVE:
                self.comment("Geometry based design rules. Radius given.")
                self.comment("Warning: assume volute cross-section is circle!")
                return pi * (lastValue ** 2)

    def __c3_theta(self):
        if self.designRule == DesignRule.PFLEIDERER:
            return Speed(self.c1 * self.r1 / self.r3).units()
        elif self.designRule == DesignRule.STEPANOFF:
            return Speed(self.c1).units()
        else:
            raise Exception("""Not implemented: geometry based definition.
                Must know A2 & r2.
                Please parse the curve definition A=A(phi) or r=r(phi),
                then code with A2 and r2.
            """)
            m2_dot = self.m_dot / 2
            c2 = m2_dot / (float(self.rho) * A2)
            return Speed(c2 * r2 / self.r3).units()

    def __Y_theta(self, c3_m=0):
        # c3_m is not known and must be solved by iteration.
        c1, r1 = self.c1, self.r1
        c3_theta, r3 = self.c3_theta, self.r3
        c3 = (c3_theta ** 2 + c3_m ** 2) ** 0.5
        Y_theta = (
            (r1 * c1 / r3 - c3_theta) /
            c3
        ) ** 2
        return Y_theta

    ###########################################################################

    def solve(self):
        
        self.__findAerodynamicalProperties()

        c3_theta = self.c2 * self.r2 / self.r3
        self.set("c3_theta", c3_theta, "angular speed at volute outlet")

        boundaryLayer = BoundaryLayerCalculation(
            b_w=sqrt(self.A1/pi)
        )

        # Find solution with iteration

        point1_total = self.parent.workingPoint.properties
        p1_total = float(point1_total.p)

        s3 = float(point1_total.s)
        p3_total = float(point1_total.p)
        rho3 = float(point1_total.rho)
        c3_m_last = 0
        c3_m = 0

        iterations = 0
        while True:
            iterations += 1
            if iterations > 100:
                raise Exception("Massflow balance @ volute station 3 not convergence.")

            c3_m_last = c3_m
            c3_m = self.m_dot / (float(rho3) * float(self.A3))

            if abs(c3_m_last / c3_m - 1) < self.RELATIVE_ERROR:
                break

            c3 = sqrt(c3_m**2 + self.c3_theta**2)

            point3 = self.queryProperties(
                Enthalpy(float(self.h1_total) - c3**2/2).units(),
                Entropy(s3).units()
            )
            p3 = float(point3.p)
            rho3 = float(point3.rho)


            Y_theta = self.__Y_theta(c3_m)
            Y = Y_theta # + Y_p

            p3_total = (p1_total + Y * p3) / (1+Y)
            point3_total = self.queryProperties(
                Pressure(p3_total).units(),
                self.h1_total
            )
            s3 = float(point3_total.s)

            # next iteration ready


        self.set("Y", Y)
        self.set("Y_theta", Y_theta)


        self.set("c3_m", c3_m)
        self.set("c3", (float(c3_m)**2 + float(c3_theta)**2) ** 0.5)

        self.set("rho3", point3.rho)
        self.set("p3", point3.p)
        self.set("p3_total", point3_total.p)
        self.set("h3", point3.h)
        self.set("h3_total", point3_total.h)
        self.set("s3", point3.s, "static entropy")
        self.set("s3_total", point3_total.s, "total, = static entropy")


        info("c1=%f, c2=%f, c3=%f" % (self.c1, self.c2, self.c3))
        info("c3_theta=%f, c3_m=%f" % (self.c3_theta, self.c3_m))
        info("Volute calculation done.")