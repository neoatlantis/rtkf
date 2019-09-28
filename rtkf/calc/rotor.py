#!/usr/bin/env python3

import sys
from logging import *

import math
from math import pi, sqrt, sin, cos, tan, asin, acos, atan
from scipy.optimize import *
from scipy import integrate
import numpy as np
from numpy import array

from ._report import CalculationReport
from ..prop import Properties
from ..prop.units import *
from .curve import DiscreteCurve, BezierCurve, MiddleOf2Curves
from .parametric_curve_distance import *
from .multivar_newton import newton_multi 
from .choked import *
from .boundary_layer import BoundaryLayerCalculation
from .disk_friction import DiskFrictionLoss


average = lambda l: sum(l) / len(l)
cot = lambda x: 1/tan(x)
acot = lambda x: atan(1/x)




def findSigma(self,
    r1, r3, phi1, beta1, N, N_SB
):
    """Calculate the slip factor."""
    sigma = 1 - abs(sin(phi1)) * sqrt(sin(beta1)) / ((N+N_SB) ** 0.7)   #(9-56)

    sigma0 = sin(19 / 180*pi + beta1/5)                                 #(9-58)
    epsilon_lim = (sigma-sigma0) / (1-sigma0)                           #(9-57)

    epsilon = r3 / r1
    if epsilon > epsilon_lim:
        xi = ((epsilon-epsilon_lim)/(1-epsilon_lim)) ** sqrt(beta1/10)  #(9-60)
        sigma_cor = sigma * (1-xi)                                      #(9-59)
        return sigma_cor
    else:
        return sigma


def findDeltaW(self, c3_theta):
    r3, r1 = self.r3, self.r1
    c1_theta = float(self.c1_theta)
    L, N, N_SB, F_SB = self.L, self.N, self.N_SB, self.F_SB

    deltaW = pi*abs(r3*c3_theta-r1*c1_theta) / (
        L * (N + N_SB*F_SB)
    )  # (9-50)
    return deltaW



def findY_CL(self,
    rho1, rho2, rho3, r1, r2, r3, b1, b2, b3, L, N, N_SB, F_SB, delta_c,
    m_dot, c1_theta, c3_theta,
    p3_total_rotate, p3
):
    """Loss coefficient of clearance"""
    weighted_average = lambda a,b,c: (a + 2*b + c) / 4                  #(9-3)
    #rho_ave = weighted_average(rho1, rho2, rho3)
    rb_ave = weighted_average(r1*b1, r2*b2, r3*b3)
    rho_ave = 1 # !!!!!?????
    #print(rho_ave)

    deltaP = m_dot * abs(r1 * c1_theta - r3 * c3_theta) / (
        rho_ave * rb_ave * L * (N + N_SB*F_SB)
    )                                                                   #(9-63)
    u_CL = sqrt(2 * deltaP / rho_ave)                                   #(9-64)
    m_dot_CL = 0.816 * rho_ave * u_CL * L * (N + N_SB*F_SB) * delta_c   #(9-65)
    Y_CL = m_dot_CL * deltaP / (m_dot * (p3_total_rotate-p3))           #(9-66)
    return Y_CL



"""def findY_Q(self,
    Q1, Q3,
    u1, c1_theta, u2, c2_theta, s1, s3, T1_total_rotate,
    rho1_total_rotate, p3_total_rotate, p3
):
    #Moisture loss coefficient
    deltaHQ = (
        (1-(Q1+Q3)/2) * abs(u1*c1_theta - u2*c2_theta) + 
        T1_total_rotate*(s3-s1)
    )
    Y_Q = deltaHQ * rho1_total_rotate / (p3_total_rotate - p3)
    return Y_Q"""



def findY(self,
    # DYNAMICS
    rho2, c2_theta,
    w3_theta, w3_m, c3_theta,

    # THERMODYNAMICAL PARAMETERS
    point1,
    point1_total_rotate,
    point1_total_absolute,

    point3,
    point3_total_rotate,
    point3_total_absolute,
):
    """Calculation the sum of all loss coefficients."""

    # TODO following values must be found

    N_SB = self.N_SB       # number of separate blades
    F_SB = self.F_SB       # fractional length of separate blades
    delta_d = getattr(self, "delta_d", 0)  # clearance from disk back to wall
    delta_c = getattr(self, "delta_c", 0)  # clearance from blade tip to wall

    # geometry

    beta1 = self.beta1
    b1, b2, b3 = self.b1, self.b2, self.b3
    r1, r2, r3 = self.r1, self.r2, self.r3
    b_th, r_th, o = self.b_th, self.r_th, self.o
    phi1, phi3, m3, L = self.phi1, self.phi3, self.m3, self.L
    N = self.N
    alpha1 = self.alpha1
    omega = self.omega
    c1_m, c1_theta = float(self.c1_m), float(self.c1_theta)
    m_dot = self.m_dot

    # thermodynamical values at in- and outlet

    u1, u2, u3 = float(self.u1), float(self.u2), float(self.u3)
    rho1, rho3 = float(point1.rho), float(point3.rho)
    rho2 = (rho1+rho3)/2
    s1, s3 = float(point1.s), float(point3.s)
    rho1_total_rotate = float(point1_total_rotate.rho)
    p1_total_rotate = float(point1_total_rotate.p)
    p1 = float(point1.p)
    p3_total_rotate = float(point1_total_rotate.p)
    p3 = float(point1.p)
    Q1, Q3 = float(point1.Q), float(point3.Q)
    T1_total_rotate = float(point1_total_rotate.T)


    # TODO calculate boundary thickness delta_b, (9-51)
    Y_p = 0 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #theta_b = theta_b * (1 + N_SB * F_SB / N)                           #(9-52)

    #print(r3*c3_theta, r1*c1_theta)

    deltaW = findDeltaW(self, c3_theta)

    sigma = findSigma(self,
        r1=r1, r3=r3, phi1=phi1, beta1=beta1, N=N, N_SB=N_SB) 
    #c1_theta_asterisk = sigma * (u1 - c1_m * cot(beta1))                #(9-53)
    c1_theta_asterisk = (1/sigma) * (u1 - c1_m * cot(beta1))             #(9-53)
    self.set("c1_theta_asterisk", c1_theta_asterisk, "optimal circ. speed")

    alpha_1_asterisk = atan(c1_m/c1_theta_asterisk)                     #(9-54)

    # Incidence loss coefficient
    Y_inc = (sin(alpha1-alpha_1_asterisk)**2) * (
        (p1_total_rotate-p1) / (p3_total_rotate-p3)
    )                                                                   #(9-55)

    # Blade load loss coefficient
    w3 = sqrt(w3_m**2 + w3_theta**2)
    Y_BL = 1/24 * (2*deltaW/w3) ** 2                                    #(9-61)
    #              ^-- Definitions of deltaW is different from Aungier's book
    #              "Centrifugal Compressors" and "Turbine Aerodynamics", where
    #              it seems the deltaW given latter in Fig.(9-5) doubles the
    #              definition in the earlier book. However (9-61) still uses
    #              the old form, so the deltaW here should also be doubled.

    # Hub-to-shroud load loss coefficient
    beta2 = 0.5 * (self.beta1 + self.beta3)
    w2 = (float(self.c1_m) + float(w3_m)) / 2 / sin(beta2)
    kappa_m = abs(phi3 - phi1) / m3

    c3_m = float(w3_m)
    sin_alpha3 = c3_m / sqrt(c3_m**2 + c3_theta**2)
    Y_HS = 1/6 * (kappa_m * b3 * w2 / (w3 * sin_alpha3)) ** 2          #(9-62)
    #Y_HS = 0# TODO INCLUDE BLADE LOADING EFFECT!!!!

    # Clearance loss, (9-63) to (9-66)
    Y_CL = findY_CL(self,
        rho1, rho2, rho3, r1, r2, r3, b1, b2, b3, L, N, N_SB, F_SB, delta_c,
        m_dot, c1_theta, c3_theta,
        p3_total_rotate, p3)

    # Moisure loss
    """Y_Q = findY_Q(self,
        Q1=Q1, Q3=Q3,
        u1=u1, c1_theta=c1_theta, u2=u2, c2_theta=c2_theta,
        s1=s1, s3=s3, T1_total_rotate=T1_total_rotate,
        rho1_total_rotate=rho1_total_rotate, 
        p3_total_rotate=p3_total_rotate,
        p3=p3
    )"""
    Y_Q = 0


    sumY = sum([Y_inc, Y_BL, Y_HS, Y_CL])

    print("Y_inc=%f Y_BL=%f Y_HS=%f Y_CL=%f   ->   Y - Y_p=%f" % (
        Y_inc, Y_BL, Y_HS, Y_CL, sumY
    ))

    return sumY





class RotorCalculation(CalculationReport):

    RELATIVE_ERROR = 1e-3
    ITERATIONS_MAX = 100

    def __init__(self, calc):
        self.parent = calc
        self.config = calc.config.rotor # input configuration
        self.queryProperties = self.parent.queryProperties

        CalculationReport.__init__(self) # recording begins
        try:
            self.override(self.config)
        except:
            self.__findGeometry()
        self.__deriveGeometry()

    def __findGeometry(self):
        """
        Find following properties

        : beta1,
        : b1, b2, b3, r1, r2, r3,
        : b_th, o,            # throat
        : phi1, phi3, m3, L,  # angles and lengthes of flow path
        : N,                  # number of blades
        N_SB,               # number of separate blades
        F_SB,               # fractional length of separate blades
        delta_d,            # clearance from disk back to wall
        delta_c,            # clearance from blade tip to wall
        """

        bladeProperties = self.config.BladeProperties
        meridian = self.config.Meridian

        TgrMer_RImpTurb = meridian.TgrMer_RImpTurb

        geoLeadingEdge = TgrMer_RImpTurb.Bezier4MerLE["GeoLeadingEdge"]
        geoTrailingEdge = TgrMer_RImpTurb.Bezier4MerTE["GeoTrailingEdge"]
        geoHub = TgrMer_RImpTurb.TCCurveMer_RImpTurb["GeoHub"]
        geoShroud = TgrMer_RImpTurb.TCCurveMer_RImpTurb["GeoShroud"]

        BladeValues = bladeProperties.TReadWriteArray_TBladeProps[
            "BladeValues"]
        TBladeProps_RImpTurb_dict = BladeValues.TBladeProps_RImpTurb
        MainBlade = TBladeProps_RImpTurb_dict[
            list(TBladeProps_RImpTurb_dict.keys())[0]]


        MainBladeMeanlines = self.config.geometry.Blade["Main"].Meanlines
        mainCurve3dHub = MainBladeMeanlines[0]
        mainCurve3dMid = MainBladeMeanlines[int(len(MainBladeMeanlines)/2)]
        mainCurve3dShr = MainBladeMeanlines[-1]

        distance = lambda array: sqrt(
            (array[0].x - array[1].x) ** 2 +
            (array[0].y - array[1].y) ** 2
        )

        pointsArray = lambda array: [np.array([e.x, e.y]) for e in array]

        self.set("N", bladeProperties.nBl, "count of blades")
        self.set("b1", distance(geoLeadingEdge.Points))
        self.set("b3", distance(geoTrailingEdge.Points))

        # ---- find betas

        beta1 = getattr(mainCurve3dMid[0], "ß") / 180 * pi
        beta3 = getattr(mainCurve3dMid[-1], "ß") / 180 * pi
        self.set("beta1", beta1, formatter=self.formatters.ANGLE)
        self.set("beta3", beta3, formatter=self.formatters.ANGLE)

        # ---- find out meridian curve r(z) for middle line

        meridianMid = DiscreteCurve([
            np.array([e.z, e.r])
            for e in mainCurve3dMid
        ])
        meridianHub = DiscreteCurve([
            np.array([e.z, e.r])
            for e in self.config.geometry.Meridian.Hub
        ])
        meridianShr = DiscreteCurve([
            np.array([e.z, e.r])
            for e in self.config.geometry.Meridian.Shroud
        ])

        hubRtoL = DiscreteCurve([np.array([e.L, e.r]) for e in mainCurve3dHub])
        hubZtoL = DiscreteCurve([np.array([e.L, e.z]) for e in mainCurve3dHub])
        shrRtoL = DiscreteCurve([np.array([e.L, e.r]) for e in mainCurve3dShr])
        shrZtoL = DiscreteCurve([np.array([e.L, e.z]) for e in mainCurve3dShr])

        self.set("r1", max(*[e.r for e in mainCurve3dMid]))
        self.set("r3", min(*[e.r for e in mainCurve3dMid]))

        self.set("L", self.__3dPathLength(mainCurve3dMid), "flow path length")
        self.set("m3", meridianMid.length, "mer. length of flow path")

        self.set("r2", (hubRtoL(0.5)[1] + shrRtoL(0.5)[1])/2)
        self.set("b2", sqrt(
            (hubRtoL(0.5)[1] - shrRtoL(0.5)[1]) ** 2 +
            (hubZtoL(0.5)[1] - shrZtoL(0.5)[1]) ** 2
        ))

        vec1 = meridianMid(0.05) - meridianMid(0.0)
        vec3 = meridianMid(1.0) - meridianMid(0.95)
        phi1 = acos(vec1[0] / norm(vec1))
        phi3 = acos(vec3[0] / norm(vec3))

        self.set("phi1", phi1, formatter=self.formatters.ANGLE)
        self.set("phi3", phi3, formatter=self.formatters.ANGLE)

        # ---- flow path related

        o, r_th, b_th = self.__findThroat(
            mainCurve3dHub, mainCurve3dMid, mainCurve3dShr)
        self.set("o", o, "blade-to-blade width of throat")
        self.set("b_th", b_th, "throat width")
        self.set("r_th", r_th, "throat radius")


        

    def __deriveGeometry(self):
        # ---- geometry regarding discharge flow

        self.comment("Geometry of relative discharge flow.")

        pitch = lambda r: 2 * pi * r / self.N
        pitch3 = pitch(self.r3)
        kappa_m = abs(self.phi3-self.phi1) / self.m3                    #(9-47)
        alpha_os_rotate = asin(self.b_th * self.o / (pitch3 * self.b3)) #(9-48)
        #r_th = r3 - o ** 2 * sin(phi3) / (2 * pitch3)                  #(9-50)
        alpha3_rotate = atan(self.r3 / self.r_th * tan(alpha_os_rotate))#(9-49)

        self.set("pitch1", pitch(self.r1))
        self.set("pitch3", pitch3)
        self.set("kappa_m", kappa_m, "average mean line curvature")
        self.set(
            "alpha_os_rotate", alpha_os_rotate,
            formatter=self.formatters.ANGLE)

        self.set(
            "alpha3_rotate", alpha3_rotate,
            formatter=self.formatters.ANGLE
        )


    def __beta(self, r):
        """Linear interpolation for beta from r1 to r3"""
        return (self.beta3 - self.beta1) / (self.r3 - self.r1) * (r - self.r1)\
            + self.beta1

    def __3dPathLength(self, curve3d):
        p0 = curve3d[0]
        l = 0
        for p1 in curve3d[1:]:
            v = norm([p1.x-p0.x, p1.y-p0.y, p1.z-p0.z])
            l += norm(v)
            p0 = p1
        return l

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

    def __findPreviousElementOutput(self):
        """Copy calculation results from upstreaming element."""
        self.comment("Copy calculation results from previous element.")
        upstreamingElement = self.parent.nozzle_row
        mapping = {
            # last element => current element
            "m_dot": "m_dot",
            "p3_total": "p1_total",
            "p3": "p1",
            "rho3": "rho1",
            "h3_total": "h1_total",
            "h3": "h1",
            "s3": "s1",
            "c3": "c1",
            "c3_m": "c1_m",
            "c3_theta": "c1_theta",
            "T3_total": "T1_total",
        }
        for lastEntry in sorted(mapping.keys()):
            self.set(
                mapping[lastEntry],
                getattr(upstreamingElement, lastEntry))
        self.set(
            "alpha1",
            atan(self.c1_m / self.c1_theta),
            formatter=self.formatters.ANGLE)
        self.set("omega", self.parent.workingPoint.omega, "rotating speed, radians/s")

    def solve(self):
        """Seeks a solution for rotor. If p3_override is given, in choking
        conditions this will be the final discharge pressure."""

        self.__findPreviousElementOutput()

        b_w = self.b3
        b_b = self.pitch3 * sin(self.beta3)
        b_avg = b_w * b_b * 2 / (b_w + b_b) 


        #boundaryLayer = BoundaryLayerCalculation(
        #    b_w=self.b3,                        # end-wall width
        #    b_b=self.pitch3*sin(self.beta3)     # blade-to-blade width
        #)
        boundaryLayer = BoundaryLayerCalculation(
            b_w=b_avg,                        # end-wall width
            b_b=b_avg     # blade-to-blade width
        )

        # ---- Known dynamics

        self.comment("Known dynamic values.")

        self.set("u1", self.omega * self.r1, "m/s")
        self.set("u2", self.omega * self.r2, "m/s")
        self.set("u3", self.omega * self.r3, "m/s")
        debug("u1=%f u2=%f u3=%f" % (self.u1, self.u2, self.u3))

        # ---- find out all thermodynamical conditions at station 1 (inlet)

        point1_total_absolute = self.queryProperties(
            Entropy(self.s1).units(),
            Enthalpy(self.h1_total).units()
        )
        point1 = self.queryProperties(
            Entropy(self.s1).units(),
            Enthalpy(self.h1).units()
        )
        
        w1_theta = float(self.c1_theta) - float(self.u1)  # TODO  double check this
        w1_m = float(self.c1_m)
        w1 = sqrt(w1_theta**2+w1_m**2)

        self.set("c1_theta", self.c1_theta)
        self.set("w1_theta", w1_theta)
        self.set("w1_m", w1_m)
        self.set("w1", w1)

        debug("m_dot=%f" % self.m_dot)
        debug("c1=%f c1_m=%f c1_theta=%f" %(self.c1, self.c1_m, self.c1_theta))
        debug("w1=%f w1_m=%f w1_theta=%f" %(self.w1, self.w1_m, self.w1_theta))

        h1_total_rotate = self.h1 + w1**2/2
        point1_total_rotate = self.queryProperties(
            Enthalpy(h1_total_rotate).units(),
            Entropy(self.s1).units()
        )

        # ---- Disk Friction Loss

        deltaH_DF = DiskFrictionLoss(
            m_dot=float(self.m_dot),
            rho=float(point1.rho),
            omega=self.omega,
            r=self.r1,
            delta=self.delta_d,
            mu=float(point1.mu)
        )

        # ---- prepare for iteration

        def calcY(
            rho2, c2_theta,
            w3_theta, w3_m, c3_theta,
            point3, point3_total_absolute, point3_total_rotate,

        ):
            return findY(self,
                rho2=rho2, c2_theta=c2_theta, c3_theta=c3_theta,
                w3_theta=w3_theta, w3_m=w3_m, 
                point1=point1, point1_total_rotate=point1_total_rotate,
                point1_total_absolute=point1_total_absolute,
                point3=point3, point3_total_rotate=point1_total_rotate,
                point3_total_absolute=point3_total_absolute
            )

        #I = float(point1_total_absolute.h) - float(self.u1) * float(self.c1_theta)
        #self.set("I", I, "rothalpy")

        # lets now figure out the ideal discharge condition
        point3_total_rotate_ideal = self.queryProperties(
            Enthalpy(
                h1_total_rotate + self.u3**2/2 - self.u1**2/2
            ).units(),
            point1.s
        )
        p3_total_rotate_ideal = float(point3_total_rotate_ideal.p)
        #print("p3_total_rotate_ideal", p3_total_rotate_ideal)

        # initial guess
        point3_rotate_ideal = self.queryProperties(
            Enthalpy(
                float(point3_total_rotate_ideal.h) - float(self.c1_m)**2/2
            ).units(),
            point3_total_rotate_ideal.s
        )
        p3_guess = float(point3_rotate_ideal.p)

        self.set("p3_total_rotate_ideal", p3_total_rotate_ideal)
        self.set("p3_guess", p3_guess, "initial guess of p3")

        """Y = None
        c3, c3_m, c3_theta, w3_theta = None, None, None, None
        point3, point3_total_rotate, point3_total_absolute = None, None, None
        isChokingSolution, m_dot_choking_rotate = False, None
        """

        print("\n********** Begin of Rotor Solution **********\n")

        # ---- MASS BALANCE PROCEDURE ----

        u1, u2, u3, r3, b3 = self.u1, self.u2, self.u3, self.r3, self.b3
        c1_theta = float(self.c1_theta)
        c1_m = float(self.c1_m)

        s3 = float(point1.s)
        w3 = 0
        A3_m = 2 * pi * self.r3 * self.b3 * sin(self.alpha3_rotate)
        Delta = 0 
        rho3 = None

        p3_total_rotate = float(point3_total_rotate_ideal.p)

        w3_last = 0
        rho3_last = float(point3_total_rotate_ideal.rho)
        rho3 = rho3_last
        chokingSuspected = False
        choked, chokedMassflow, chokedPressure = False, 0, 0
        iterations = 0
        while True:
            iterations += 1
            if iterations > self.ITERATIONS_MAX:
                critical("Iteration not finding convergence. Exit now.")
                print(self.parent.workingPoint.variation)
                exit()

            if rho3 is None:
                rho3 = float(point3_total_rotate_ideal.rho)

            A3_m_reduced = A3_m * (1 - Delta)

            w3_last = w3
            rho3_last = rho3

            if chokingSuspected:
                # Exit the iteration when the current speed is near enough
                # to the sound speed.
                # This is important! If not and go into another iteration with
                # new c3 calculated based on this situation, it would go over
                # the limit and cause numerical error.
                if abs(w3 / a3 - 1) < self.RELATIVE_ERROR:
                    warning("Choking occured! Rel. err=%e" % abs(w3/a3-1))
                    choked = True
                    chokedPoint = getChokingPoint(self, point3_total_rotate)
                    chokedMassflow = \
                        self.N * self.b_th * self.o * (1-Delta) * \
                        float(chokedPoint.rho) * float(chokedPoint.a)
                    print(Delta)
                    chokedPressure = float(point3.p)
                    break

                w3 = self.m_dot / (rho3 * A3_m_reduced)
                if abs(w3) > a3:
                    warning("Limiting w3=%f < a3=%f" % (w3, a3)) 
                    w3 = a3
            else:
                w3 = self.m_dot / (rho3 * A3_m_reduced)

            debug("A3_m_reduced: %f, Delta: %f" % (A3_m_reduced, Delta))

            w3_theta = w3 * cos(self.alpha3_rotate)
            w3_m = w3 * sin(self.alpha3_rotate)
            debug("w3=%f w3_m=%f w3_theta=%f alpha3'=%f" %
                (w3, w3_m, w3_theta, self.alpha3_rotate*180/pi))

            print(rho3_last, w3_last, rho3, w3)
            print("!!!!!!!!!!!!!!")
            debug("choking check %f" % (
                (rho3_last*w3_last-rho3*w3) * 
                (w3_last - w3)
            ))
            if (
                (rho3_last*w3_last-rho3*w3) * 
                (w3_last - w3)
            ) < 0:
                chokingSuspected = True
                warning("Choked solution suspected.")

            c3_m = w3_m
            c3_theta = self.u3 - w3_theta # TODO DOUBLE CHECK
            c3 = sqrt(c3_m**2 + c3_theta**2)

            debug("c3=%f c3_m=%f c3_theta=%f" % (c3, c3_m, c3_theta))
            debug("u3=%f" % u3)

            c2_theta = 0.5 * (float(self.c1_theta) + c3_theta)
            # TODO use mass balance to achieve c2_theta

            work = c1_theta * u1 - c3_theta * u3 
            debug("Circu. work = %f" % work)
            debug("------------> %f x %f - %f x %f = %f" % (
                c1_theta, u1, c3_theta, u3, work
            ))


            h3_total_absolute = float(point1_total_absolute.h) - work
            h3 = h3_total_absolute - c3**2/2 + deltaH_DF
            h3_total_rotate = h3 + w3**2/2

            debug("h3*=%f, h3*rel=%f, h3=%f" % (
                h3_total_absolute, h3_total_rotate, h3
            ))

            # determine calculated point3 static and total conditions            

            try:
                point3_total_rotate = self.queryProperties(
                    Pressure(p3_total_rotate).units(),
                    Enthalpy(h3_total_rotate).units())

                s3_old = s3
                s3 = point3_total_rotate.s

                debug("s3: %f -> %f" % (s3_old, s3))

                point3 = self.queryProperties(
                    s3,
                    Enthalpy(h3).units()
                )
                p3 = float(point3.p)
                rho3_last = rho3
                rho3 = float(point3.rho)
                debug("New density: %f (last density: %f)" % (rho3, rho3_last))
                a3 = float(point3.a)
            except:
                print(self)
                exit()

            point3_total_absolute = self.queryProperties(
                s3,
                Enthalpy(h3_total_absolute).units())

            # now p3, rho3, a3 updated. calculate loss coefficients.

            deltaW = findDeltaW(self, c3_theta)

            debug("deltaW=%f" % deltaW)
            
            beta2 = 0.5 * (self.beta1 + self.beta3)
            w2 = 0.5 * (self.w1_m + w3_m) / sin(beta2)

            Delta, Y_p = boundaryLayer.clear()\
                .setSplitBlades(N=self.N, N_SB=self.N_SB, F_SB=self.F_SB)\
                .setBladeLoadingSpeedDifference(deltaW)\
                .setLastDelta(Delta)\
                .addEndwall(
                    mu=point3.mu, L=self.L,
                    u1=self.c1, u3=c3,
                    rho1=self.rho1, rho3=point3.rho
                )\
                .addBlade2Blade(
                    mu=point3.mu, L=self.L,
                    u1=self.w1, u2=w2, u3=w3,
                    rho1=self.rho1, rho3=point3.rho,
                    role="other"
                )\
                .addBlade2Blade(
                    mu=point3.mu, L=self.L,
                    u1=self.w1, u2=w2, u3=w3,
                    rho1=self.rho1, rho3=point3.rho,
                    role="suction"
                )\
                .addBlade2Blade(
                    mu=point3.mu, L=self.L,
                    u1=self.w1, u2=w2, u3=w3,
                    rho1=self.rho1, rho3=point3.rho,
                    role="pressure"
                )\
                ()
            #Delta, Y_p = 0,0 ######################
            

            Y = Y_p + calcY(
                rho2=float(point3.rho + point1.rho) / 2,
                c2_theta=c2_theta,
                w3_theta=w3_theta,
                w3_m=w3_m,
                c3_theta=c3_theta,
                point3=point3,
                point3_total_absolute=point3_total_absolute,
                point3_total_rotate=point3_total_rotate
            )


            p3_total_rotate = (p3_total_rotate_ideal + Y * p3) / (1+Y)

            
            e_m = abs(self.m_dot / (rho3 * A3_m_reduced * w3) - 1)
            debug("Rel. error = %E" % e_m)

            if e_m < self.RELATIVE_ERROR:
                break

            # check for \partial{rho3*c3_m}/\partial{c3_m}
            debug("choking check %f" % (
                (rho3_last*w3_last-rho3*w3) * 
                (w3_last - w3)
            ))
            if (
                (rho3_last*w3_last-rho3*w3) * 
                (w3_last - w3)
            ) < 0:
                chokingSuspected = True
                warning("Choked solution suspected.")

            info("Iteration #%d done.\n" % iterations)
            

        # ---- END OF MASS BALANCE PROCEDURE
        
        self.set("choked", choked)
        self.set("Y", Y)
        self.set(
            "Delta", Delta, "percentage of boundary layer",
            formatter=self.formatters.PERCENT
        )
        self.set("work", work, "circumferential work")
        self.set("s3", point3.s)
        self.set("h3_total_rotate", point3_total_rotate.h)
        self.set("h3_total", point3_total_absolute.h)
        self.set("p3_total", point3_total_absolute.p)
        self.set("T3_total", point3_total_absolute.T)

        
        if choked:
            info("Rotor is choked.")
            debug("Choked massflow = %f" % chokedMassflow)
            
            self.set("h3_choked", point3.h)
            self.set("p3_choked", point3.p)
            self.set("m_dot", chokedMassflow, "choked mass flow")
            
            point3_after = None
            def calc(c1_m_ring):
                """The minimum pressure is given by assuming the discharge
                ring passage has got a meridian velocity == a(must still
                consider circum. velocity!)."""
                nonlocal point3_after
                h = float(point3_total_absolute.h) -\
                    (u3**2+c1_m_ring**2)/2  # w3_theta=0 -> c3_theta=u3
                point3_after = self.queryProperties(
                    Enthalpy(h).units(),
                    point3_total_rotate.s
                )
                return float(point3_after.a) - c1_m_ring
            c1_m_ring = newton(calc, 0)

            raise ChokedException(
                massflow=chokedMassflow,
                p_max=chokedPressure,
                p_min=point3_after.p,
                where="rotor"
            )

        info("Rotor is NOT choked.")
        self.finalize(w3, self.alpha3_rotate, point3)

    def setChokedAfterPressure(self, p3_after):
        assert self.choked == True

        A3_ring = (1-self.Delta) * 2*pi * self.r3 * self.b3

        # do a corrected solution with shockwave expansion.
        point3_total_rotate = self.queryProperties(
            self.s3, self.h3_total_rotate
        )
        point3_total_absolute = self.queryProperties(
            self.s3, self.h3_total
        )
        point3_after = self.queryProperties(
            Pressure(p3_after).units(), self.s3
        )

        w3 = sqrt(2*(float(point3_total_rotate.h) - float(point3_after.h)))
        alpha3_rotate = asin(self.m_dot / (
            A3_ring * float(point3_after.rho) * w3
        ))
        self.set(
            "alpha3_rotate", alpha3_rotate,
            "relative supersonic discharge angle !",
            formatter=self.formatters.ANGLE
        )

        self.finalize(w3, alpha3_rotate, point3_after)
       

    def finalize(self, w3, alpha3_rotate, point3):
        w3_m = w3 * sin(alpha3_rotate)
        w3_theta = w3 * cos(alpha3_rotate)
        c3_m = w3_m
        c3_theta = self.u3 - w3_theta ## CHECK THIS!
        c3 = sqrt(c3_m**2 + c3_theta**2)
        
        self.set("w3_theta", w3_theta)
        self.set("w3_m", w3_m)
        self.set("w3", w3)
        
        self.set("c3", c3)
        self.set("c3_m", c3_m)
        self.set("c3_theta", c3_theta)
        
        self.set("p3", point3.p)
        self.set("h3", point3.h)
        self.set("rho3", point3.rho)
        self.set("a3", point3.a, "sonic speed")

        point3_ideal = self.queryProperties(
            self.parent.workingPoint.properties.s,
            point3.p
        )

        self.set("h3_ideal", point3_ideal.h, "ideal entropy at outlet")
        self.set(
            "dh_total2static",
            (
                self.parent.workingPoint.properties.h -
                point3_ideal.h
            )
        )
        self.set(
            "eta_total2static",
            (self.work / float(self.dh_total2static)),
            formatter=self.formatters.PERCENT
        )
        self.set(
            "speed_ratio",
            self.u1 / sqrt(self.dh_total2static),
            "u1/c0",
            formatter=self.formatters.PERCENT
        )

        print("\n")
        print("********** End of Rotor Solution ***********")
        print("\n")

        
