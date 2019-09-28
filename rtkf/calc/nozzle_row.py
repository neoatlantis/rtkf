#!/usr/bin/env python3

from logging import *

import math
from math import pi, sqrt, sin, cos, tan, asin, acos, atan, copysign
from scipy.optimize import *
from scipy import integrate
from numpy import array

from ._report import CalculationReport
from ..prop import Properties
from ..prop.units import *
from .curve import DiscreteCurve, BezierCurve, MiddleOf2Curves
from .parametric_curve_distance import *
from .multivar_newton import newton_multi 
from .choked import *
from .boundary_layer import BoundaryLayerCalculation


average = lambda l: sum(l) / len(l)

def findDeltaC(self, c3_theta):
    r3, r1 = self.r3, self.r1
    c1_theta = float(self.c1_theta)
    L, N = self.L, self.N

    deltaC = pi*abs(r3*c3_theta-r1*c1_theta) / (L*N)
    return deltaC




class NozzleRowCalculation(CalculationReport):

    ITERATIONS_MAX = 100
    RELATIVE_ERROR = 5e-4

    def __init__(self, calc):
        self.parent = calc
        self.config = calc.config.nozzleRow # input configuration
        self.queryProperties = self.parent.queryProperties

        CalculationReport.__init__(self) # recording begins
        try:
            self.override(self.config)
        except:
            self.__findGeometry()
        self.__deriveGeometry()

    def __findGeometry(self):
        self.comment("Geometry parameters.\n")

        MainDimensionsElement = self.config.MainDimensionsElement
        BladeProperties = self.config.BladeProperties

        MainBladeMeanlines = self.config.geometry.Blade["Main"].Meanlines
        mainCurve3dHub = MainBladeMeanlines[0]
        mainCurve3dMid = MainBladeMeanlines[int(len(MainBladeMeanlines)/2)]
        mainCurve3dShr = MainBladeMeanlines[-1]

        hubRtoL = DiscreteCurve([np.array([e.L, e.r]) for e in mainCurve3dHub])
        hubZtoL = DiscreteCurve([np.array([e.L, e.z]) for e in mainCurve3dHub])
        midRtoL = DiscreteCurve([np.array([e.L, e.r]) for e in mainCurve3dMid])
        midZtoL = DiscreteCurve([np.array([e.L, e.z]) for e in mainCurve3dMid])
        shrRtoL = DiscreteCurve([np.array([e.L, e.r]) for e in mainCurve3dShr])
        shrZtoL = DiscreteCurve([np.array([e.L, e.z]) for e in mainCurve3dShr])

        MerInlet = MainDimensionsElement.MerData["MerInlet"]
        MerOutlet = MainDimensionsElement.MerData["MerOutlet"]

        self.comment("Assuming nozzle row is rotation-symmetric on axis x=0!")

        r1 = max(*[e.r for e in mainCurve3dMid])
        b1 = abs(MerInlet.Shr.x - MerInlet.Hub.x)
        r3 = min(*[e.r for e in mainCurve3dMid])
        b3 = abs(MerOutlet.Shr.x - MerOutlet.Hub.x)


        BladeValues = BladeProperties.TReadWriteArray_TBladeProps[
            "BladeValues"]
        TBladePropsStator_dict = BladeValues.TBladePropsStator
        MainBlade = TBladePropsStator_dict[
            list(TBladePropsStator_dict.keys())[0]]


        t_b1 = average(MainBlade.s1)
        t_b3 = average(MainBlade.s2)

        beta1 = getattr(mainCurve3dMid[0], "ß") / 180 * pi
        beta2 = getattr(mainCurve3dMid[int(len(mainCurve3dMid)/2)], "ß") / 180 * pi
        beta3 = getattr(mainCurve3dMid[-1], "ß") / 180 * pi

        N = BladeProperties.nBl
        t_b = 0.5 * (t_b1+t_b3)

        self.set("N", N, "blades count")
        self.set("t_b", t_b, "blade thickness")

        self.comment("At inlet:")
        self.set("r1", r1, "radius")
        self.set("b1", b1, "width")
        self.set(
            "beta1", beta1, "blade camberling angle",
            formatter=self.formatters.ANGLE)

        self.set("beta2", beta2, formatter=self.formatters.ANGLE)

        self.comment("At outlet:")
        self.set("r3", r3)
        self.set("b3", b3)
        self.set("beta3", beta3, formatter=self.formatters.ANGLE)

        self.comment("Calculated:")
        self.set("L", mainCurve3dMid[-1].L, "flow path length")

        o, r_th, b_th = self.__findThroat(
            mainCurve3dHub, mainCurve3dMid, mainCurve3dShr)
        self.set("o", o, "blade-to-blade width of throat")
        self.set("r_th", r_th, "radius at throat")
        self.set("b_th", b_th, "throat width")        


    def __deriveGeometry(self):
        self.comment("Derived geometry:")

        r1, r3, beta1, beta3 = self.r1, self.r3, self.beta1, self.beta3
        b1, b3, b_th = self.b1, self.b3, self.b_th
        o = self.o
        beta1, beta3 = self.beta1, self.beta3
        t_b = self.t_b
        N = self.N

        self.set("pitch1", self.__pitch(r1), "blade pitch")
        self.set("pitch3", self.__pitch(r3))

        #self.set("A1",
        #    b1 * (2*pi * r1 * sin(beta1) - t_b * N), "passage area")
        #self.set("A3",
        #    b3 * (2 * pi * r3 * sin(beta3) - t_b * N), "passage area")

        alpha_os = asin(self.b_th * self.o / self.pitch3 / self.b3)
        self.set(
            "alpha3", atan(self.r3 / self.r_th * tan(alpha_os)),
            "subsonic discharge angle",
            formatter=self.formatters.ANGLE)

        self.__findOptimiumNozzleIncidenceAngle()


    def __pitch(self, r):
        return 2 * pi * r / self.N

    def __findThroat(self, curve3dHub, curve3dMid, curve3dShr):
        delta = 2 * pi / self.N

        i1, i2, omid = findThroatFromCurve3d(delta, curve3dMid)

        # TODO substract omid with blade thickness!

        M_th = curve3dMid[i1].M
        r_th = (curve3dMid[i1].r + curve3dMid[i2].r)/2
        t = M_th / curve3dMid[-1].M

        hubRtoM = DiscreteCurve([np.array([e.M, e.r]) for e in curve3dHub])
        hubZtoM = DiscreteCurve([np.array([e.M, e.z]) for e in curve3dHub])
        shrRtoM = DiscreteCurve([np.array([e.M, e.r]) for e in curve3dShr])
        shrZtoM = DiscreteCurve([np.array([e.M, e.z]) for e in curve3dShr])

        hubR, hubZ = hubRtoM(t)[1], hubZtoM(t)[1]
        shrR, shrZ = shrRtoM(t)[1], shrZtoM(t)[1]

        bmid = sqrt((shrR-hubR)**2 + (shrZ-hubZ)**2)
        return omid, r_th, bmid

    def __findOptimiumNozzleIncidenceAngle(self):
        """Find the optimium nozzle incidence angle."""
        t_2 = self.t_b

         # R.H.Aungier, Turbine Aerodynamics, Eq.(9-35)
        i_asterisk = (
            3.6 * (pi/180) * sqrt(10 * t_2 / self.L) + 
            abs(self.beta3 - self.beta1) / 3.4
        ) * sqrt(self.L / self.pitch3) - 0.5 * abs(self.beta3 - self.beta1)

        sgn = lambda i: copysign(1, i)
        alpha_asterisk = self.beta1 - i_asterisk * sgn(self.beta3 - self.beta1)
        self.comment("Optimium incidence angle.")
        self.set(
            "alpha_asterisk", alpha_asterisk, formatter=self.formatters.ANGLE)


    def __findPreviousElementOutput(self):
        """Copy calculation results from upstreaming element."""
        self.comment("Copy calculation results from previous element.")
        workingPoint = self.parent.workingPoint
        point1_total = workingPoint.properties
        upstreamingElement = self.parent.volute

        if hasattr(upstreamingElement, "SKIP"):
            info("Volute Calculation skipped, use working point properties.")
            alpha_asterisk = self.alpha_asterisk

            self.set("m_dot", workingPoint.m_dot)
            self.set("p1_total", point1_total.p)
            self.set("h1_total", point1_total.h)
            self.set("s1_total", point1_total.s)
            self.set("s1", point1_total.s)
            A1 = 2 * pi * self.r1 * self.b1
            
            iterations = 0
            rho1 = float(point1_total.rho)
            c1_theta = 0
            c1_m = 0
            while True:
                iterations += 1
                if iterations > 100:
                    raise Exception("Massflow balance @ stator not convergence.")

                c1_m_last = c1_m
                c1_m = self.m_dot / (rho1 * A1)
                if abs(c1_m_last / c1_m - 1) < self.RELATIVE_ERROR:
                    break

                c1 = c1_m / sin(self.beta1) #alpha_asterisk)

                h1 = float(point1_total.h) - c1**2/2
                point1 = self.queryProperties(
                    point1_total.s,
                    Enthalpy(h1).units()
                )
                rho1 = float(point1.rho)
            self.set("p1", point1.p)
            self.set("rho1", point1.rho)
            self.set("c1", c1)
            self.set("c1_m", c1_m)
            self.set("c1_theta", sqrt(c1**2-c1_m**2))
            info("Nozzle row entry mass balance done.")


        else:
            mapping = {
                # last element => current element
                "m_dot": "m_dot",
                "p3_total": "p1_total",
                "p3": "p1",
                "rho3": "rho1",
                "h3_total": "h1_total",
                "s3": "s1",
                "c3": "c1",
                "c3_m": "c1_m",
                "c3_theta": "c1_theta",
            }

            for lastEntry in sorted(mapping.keys()):
                self.set(
                    mapping[lastEntry],
                    getattr(upstreamingElement, lastEntry))
            self.set(
                "alpha1",
                atan(float(self.c1_m)/float(self.c1_theta)),
                formatter=self.formatters.ANGLE
            )


    def solve(self):
        """Solve for given component. If p3_after is given, supersonic solution
        will use this value as static pressure after shockwave. Otherwise,
        a ChokedException will be thrown with suggested choking massflow and
        range(maximum and minimum of p3_after)."""
        
        self.__findPreviousElementOutput()

        boundaryLayer = BoundaryLayerCalculation(
            b_w=self.b3,                        # end-wall width
            b_b=self.pitch3*sin(self.beta3)     # blade-to-blade width
        )

        alpha1 = atan(float(self.c1_m) / float(self.c1_theta))
        p1_total = float(self.p1_total)
        p1 = float(self.p1)

        def calcY_inc(p3_total, p3, c3):    
            Y_inc = (sin(alpha1 - self.alpha_asterisk)**2) * (
                    (p1_total - p1) / (float(p3_total) - float(p3))
                )
            return Y_inc


        ## ---- MASS BALANCE PROCEDURE ----

        point1_total = self.queryProperties(
            Enthalpy(self.h1_total).units(),
            Pressure(self.p1_total).units(),
        )
        s3 = float(point1_total.s)
        c3 = 0
        A3_m = 2 * pi * self.r3 * self.b3 * sin(self.alpha3) # *c3*rho3=>m_dot
        Delta = 0 
        rho3 = None

        point3_total = self.queryProperties(
            point1_total.h,       # always
            Entropy(s3).units(),  # will shift due to iteration
        )

        c3_last, rho3_last = 0, float(point3_total.rho)
        chokingSuspected = False
        choked, chokedMassflow, chokedPressure = False, 0, 0
        iterations = 0
        while True:
            iterations += 1
            if iterations > self.ITERATIONS_MAX:
                critical("Iteration not finding convergence. Exit now.")
                exit()

            if rho3 is None: 
                rho3 = float(point3_total.rho)

            A3_m_reduced = A3_m * (1-Delta)

            # record results from last iteration
            c3_last = c3
            rho3_last = rho3

            if chokingSuspected:
                # Exit the iteration when the current speed is near enough
                # to the sound speed.
                # This is important! If not and go into another iteration with
                # new c3 calculated based on this situation, it would go over
                # the limit and cause numerical error.
                if abs(c3 / a3 - 1) < self.RELATIVE_ERROR:
                    # Exit the iteration, thermodynamical parameters are now
                    # correct, mass flow may not.
                    warning("Choking occured!")
                    choked = True
                    #chokedMassflow = rho3 * A3_m_reduced * a3
                    chokedPoint = getChokingPoint(self, point3_total)
                    chokedMassflow = \
                        self.N * self.b_th * self.o * (1-Delta) * \
                        float(chokedPoint.rho) * float(chokedPoint.a)
                    chokedMassflow = min(chokedMassflow, self.m_dot)
                    chokedPressure = float(point3.p)
                    break

                c3 = self.m_dot / (rho3 * A3_m_reduced)
                if c3 > a3:
                    warning("Limiting c3=%f < a3=%f" % (c3, a3)) 
                    c3 = a3
            else:
                c3 = self.m_dot / (rho3 * A3_m_reduced)

            c3_theta = c3 * cos(self.alpha3)
            c3_m = c3 * sin(self.alpha3)

            c2 = (c3_m + self.c1_m) / (2 * sin(self.beta2))

            debug("c3_m=%f, c3_theta=%f, c3=%f" % (c3_m, c3_theta, c3))
            debug("s3=%f" % s3)

            h3 = float(point3_total.h) - c3**2/2

            point3 = self.queryProperties(
                Enthalpy(h3).units(),
                Entropy(s3).units()
            )
            p3 = float(point3.p)
            rho3 = float(point3.rho)
            a3 = float(point3.a)

            Y_inc = calcY_inc(
                p3_total=float(point3_total.p),
                p3=float(point3.p),
                c3=c3
            )

            deltaC = findDeltaC(self, c3_theta)

            Delta, Y_p = boundaryLayer.clear()\
                .setBladeLoadingSpeedDifference(deltaC)\
                .setLastDelta(Delta)\
                .addEndwall(
                    mu=point3.mu, L=self.L,
                    u1=self.c1, u3=c3,
                    rho1=self.rho1, rho3=point3.rho
                )\
                .addEndwall(
                    mu=point3.mu, L=self.L,
                    u1=self.c1, u3=c3,
                    rho1=self.rho1, rho3=point3.rho
                )\
                .addBlade2Blade(
                    mu=point3.mu, L=self.L,
                    u1=self.c1, u2=c2, u3=c3,
                    rho1=self.rho1, rho3=point3.rho,
                    role="suction"
                )\
                .addBlade2Blade(
                    mu=point3.mu, L=self.L,
                    u1=self.c1, u2=c2, u3=c3,
                    rho1=self.rho1, rho3=point3.rho,
                    role="pressure"
                )\
                ()

            Y = Y_p + Y_inc
            p3_total= (p1_total + Y*p3) / (1+Y)

            debug("Y(%f)=Y_p(%f)+Y_inc(%f)" % (Y, Y_p, Y_inc))

            point3_total = self.queryProperties(
                Enthalpy(self.h1_total).units(),
                Pressure(p3_total).units()
            )
            s3 = float(point3_total.s)


            e_m = abs(self.m_dot / (rho3 * A3_m_reduced * c3) - 1)

            debug("Rel. error = %E" % e_m)

            if e_m < self.RELATIVE_ERROR:
                break

            # check for \partial{rho3*c3_m}/\partial{c3_m}
            if (
                (rho3_last*c3_last-rho3*c3) * 
                (c3_last - c3)
            ) < 0:
                chokingSuspected = True
                warning("Choked solution suspected.")

            info("Iteration #%d done.\n" % iterations)

        # ---- END OF ITERATION ----

        self.comment("Iteration error: %s%%" % (e_m*100))
        print("e_m", e_m)

        # ---- SUPERSONIC SOLUTION CORRECTION ----

        self.set("choked", choked)
        self.set("Y", Y)
        self.set("Y_inc", Y_inc)
        self.set("Y_p", Y_p)
        self.set(
            "Delta", Delta, "percentage of boundary layer",
            formatter=self.formatters.PERCENT
        )

        self.set("s3", point3.s)
        self.set("h3_total", point3_total.h)
        self.set("p3_total", point3_total.p)
        self.set("T3_total", point3_total.T)


        if choked:
            info("Nozzle row is choked.")
            debug("Choked massflow: %f <= %f" % (chokedMassflow, self.m_dot))

            self.set("h3_choked", point3.h)
            self.set("p3_choked", point3.p)
            self.set("m_dot", chokedMassflow, "choked mass flow")



            point3_after = None
            def calc(c1_m_ring):
                """The minimum pressure is given by assuming the discharge
                ring passage has got a meridian velocity == a(must still
                consider circum. velocity!)."""
                nonlocal point3_after
                h = float(point3_total.h) - (c3_theta**2+c1_m_ring**2)/2
                point3_after = self.queryProperties(
                    Enthalpy(h).units(),
                    point3_total.s
                )
                return float(point3_after.a) - c1_m_ring
            c1_m_ring = newton(calc, 0)

            raise ChokedException(
                massflow=chokedMassflow,
                p_max=chokedPressure,
                p_min=point3_after.p,
                where="nozzle_row"
            )

            # Choked solution terminates here. After pressure between p_max
            # and p_min is required to continue calculation, which could be
            # decided by external analysis program and then given to
            # self.setChokedAfterPressure method to resume calculation.

        info("Nozzle row is NOT choked.")
        self.finalize(c3, self.alpha3, point3)



    def setChokedAfterPressure(self, p3_after):
        assert self.choked == True

        # do a corrected solution with shockwave expansion.
        point3_total = self.queryProperties(self.s3, self.p3_total)
        point3_after = self.queryProperties(
            Pressure(p3_after).units(),
            self.s3
        )
        A3_ring = (1-self.Delta) * 2*pi * self.r3 * self.b3
        c3 = sqrt(2*(float(point3_total.h) - float(point3_after.h)))
        alpha3 = asin(self.m_dot / (
            A3_ring * float(point3_after.rho) * c3
        ))
        self.set(
            "alpha3", alpha3, "supersonic discharge angle !",
            formatter=self.formatters.ANGLE
        )
        self.finalize(c3, alpha3, point3_after)


    def finalize(self, c3, alpha3, point3):
        c3_theta = c3 * cos(alpha3)
        c3_m = c3 * sin(alpha3)

        self.set("c3_theta", Speed(c3_theta).units())
        self.set("c3_m", Speed(c3_m).units())
        self.set("c3", Speed(c3).units())
        self.set("h3", point3.h)
        self.set("p3", point3.p)
        self.set("rho3", point3.rho)
        self.set("a3", point3.a, "sonic speed")

        print("\n")
        print("********** End of Nozzle Row Solution ***********")
        print("\n")
