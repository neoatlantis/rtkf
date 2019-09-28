#!/usr/bin/env python3

from enum import Enum
from math import log, sqrt
from scipy.optimize import *



class FlowState:
    LAMINAR   = "LAMINAR" 
    TRANSIENT = "TRANSIENT"
    TURBULENT = "TURBULENT"

def __cf_s(Re_d):
    def f(x):
        return 1 / sqrt(4*x) + 2*log(2.51/(Re_d*sqrt(4*x)), 10)
    return newton(f, x0=1e-4)

def boundaryLayerThickness(
    mu,                 # shear viscosity / dynamic viscosity
    d,                  # width of flow passage
    L,                  # flow length
    u1, u2, u3,         # relative speed of flow to component
    rho1, rho2, rho3,   # density
    e=0                 # roughness peak-to-valley, default: hydraulic smooth!
):
    #print("... calc boundary layer thickness: u1=%f u2=%f u3=%f" % (
    #    u1, u2, u3
    #))

    Re_d = float(rho3) * float(u3) * float(d) / float(mu)
    rho_ave = 0.25 * float(rho1 + 2*rho2 + rho3)
    H = 9 / 7 # 1.2857

    # 1. find cf for skin friction coefficient
    cf_l = 16 / Re_d  
    if Re_d < 2000: # laminar
        cf = cf_l
        flowState = FlowState.LAMINAR 
    else:           # transient or turbulent
        # Find cf_t (t: turbulent)
        Re_e = (Re_d - 2000) * e / d
        cf_s = __cf_s(Re_d)
        if Re_e <= 60:  # hydraulic smooth
            cf_t = cf_s
        else:
            cf_r = (
                1 / (-4 * log(e / (3.71*d), 10))
            )**2
            cf_t = cf_s + (cf_r - cf_s) * (1 - 60 / Re_e)

        if Re_d < 4000: # transient
            cf = cf_l + (cf_t - cf_l) * (Re_d / 2000 - 1)
            flowState = FlowState.TRANSIENT
        else: # fully turbulent
            cf = cf_t
            flowState = FlowState.TURBULENT

    # 2. momentum loss thickness
    n=5
    theta = cf * rho_ave * (L/8/rho3) * ((u1/u3)**n + 2*(u2/u3)**n + 1)

    # 3. displacement thickness
    delta_asterisk = H * theta

    # 4. boundary thickness
    delta = theta * H * (H+1)/(H-1)

    #print("... Re_d=%f cf=%f" % (Re_d, cf)) #, rho_ave, u1, u2, u3)
    #print("... delta=%f delta*=%f delta**=%f" % (
    #    delta, delta_asterisk, theta
    #))

    return {
        "Re_d": Re_d,
        "theta": theta,                     # momentum loss thickness
        "delta_asterisk": delta_asterisk,   # boundary displacement thickness
        "delta": delta,                     # boundary (speed) thickness
        "c_f": cf,
        "flow_state": flowState,
    }



class BoundaryLayerCalculation:

    def __init__(self, b_w=0, b_b=0):
        self.b_w = b_w
        self.b_b = b_b
        self.clear()

    def clear(self):
        self.theta_w = [] # momentum thickness (divided by b_w), end-wall
        self.theta_b = [] #                    (           b_b), relative
        self.delta_asterisk_w = [] # defect thickness (divided by b_w), end-wall
        self.delta_asterisk_b = [] #                  (           b_b), relative
        self.delta_w = [] # boundary speed thickness
        self.delta_b = []
        self.adjust_split_blades = 0
        self.du2 = 0
        self.lastDelta = 0
        return self

    def __calcCoeff(self, b, **args):
        for e in args:
            args[e] = float(args[e]) if args[e] is not None else args[e]
        args["u3"] /= (1 - self.lastDelta)
        boundaryLayer = boundaryLayerThickness(
            mu=args["mu"],
            d=b,
            L=args["L"],
            u1=args["u1"], u3=args["u3"],
            u2=(
                args["u2"]\
                if args["u2"] is not None\
                else (args["u1"] + args["u3"])/2
            ) + args["du2"],
            rho1=args["rho1"], rho3=args["rho3"],
            rho2=args["rho2"]\
                if args["rho2"] is not None else (args["rho1"]+args["rho3"])/2,
            e=0
        )
        delta = boundaryLayer["delta"]
        theta = boundaryLayer["theta"]
        delta_asterisk = boundaryLayer["delta_asterisk"]
        return (theta / b, delta_asterisk / b, delta)

    def setSplitBlades(self, N, N_SB, F_SB):
        self.adjust_split_blades = N_SB * F_SB / N
        return self

    def setBladeLoadingSpeedDifference(self, du2):
        self.du2 = abs(du2)
        return self

    def setLastDelta(self, lastDelta):
        self.lastDelta = lastDelta
        return self

    def addEndwall(self, mu, L, u1, u3, rho1, rho3, u2=None, rho2=None):
        #print('endwall')
        theta_w, delta_asterisk_w, delta_w = self.__calcCoeff(
            b=self.b_w,
            mu=mu, L=L, u1=u1, u2=u2, u3=u3, rho1=rho1, rho2=rho2, rho3=rho3,
            du2=0
        )
        self.theta_w.append(theta_w)
        self.delta_asterisk_w.append(delta_asterisk_w)
        self.delta_w.append(delta_w)
        #print("BOUNDARY LAYER: WALL :: u1=%s  u2=%s  u3=%s" % (
        #    u1, u2, u3
        #))
        print("BOUNDARY LAYER: WALL :: delta=%f delta*=%f delta**=%f" % (
            delta_w, delta_asterisk_w, theta_w
        ))
        return self

    def addBlade2Blade(self, mu, L, u1, u3, rho1, rho3, role, u2=None, rho2=None):
        assert role in ["suction", "pressure", "other"]
        du2 = {
            "suction":  self.du2,
            "pressure": -self.du2,
            "other": 0
        }[role]
        #print(role)
        #print("****", float(u1),  float(u3), u3)
        theta_b, delta_asterisk_b, delta_b = self.__calcCoeff(
            b=self.b_b,
            mu=mu, L=L, u1=u1, u2=u2, u3=u3, rho1=rho1, rho2=rho2, rho3=rho3,
            du2=du2
        )
        self.theta_b.append(theta_b)
        self.delta_asterisk_b.append(delta_asterisk_b)
        self.delta_b.append(delta_b)
        #print("BOUNDARY LAYER: BLADE : u1=%f  u2=%f  u3=%f" % (
        #    u1, u2, u3
        #))
        #print("BOUNDARY LAYER: BLADE : delta=%f delta*=%f delta**=%f" % (
        #    delta_b, delta_asterisk_b, theta_b
        #))
        return self

    def __call__(self):
        sum_theta_w = sum(self.theta_w)
        sum_theta_b = sum(self.theta_b)
        sum_delta_asterisk_w = sum(self.delta_asterisk_w)
        sum_delta_asterisk_b = sum(self.delta_asterisk_b)
        sum_delta_w = sum(self.delta_w)
        sum_delta_b = sum(self.delta_b)

        sum_theta_b *= (1 + self.adjust_split_blades)

        #print("***", sum_delta_asterisk_w, sum_theta_w)

        if sum_delta_w / self.b_w > 1:
            sum_delta_asterisk_w *= self.b_w / sum_delta_w
            sum_theta_w *= self.b_w / sum_delta_w
        if sum_delta_b / self.b_b > 1:
            sum_delta_asterisk_b *= self.b_b / sum_delta_b
            sum_theta_b *= self.b_b / sum_delta_b

        #print("***", sum_delta_asterisk_w, sum_theta_w)


        Theta = 1 - (1-sum_theta_w) * (1-sum_theta_b)
        Delta = 1 - (1-sum_delta_asterisk_w) * (1-sum_delta_asterisk_b)


        if Theta > 1:
            Delta = Theta - 1
            Theta = 1

        Y_p = (2*Theta + Delta**2) / (1 - Delta)**2
        #print("BOUNDARY LAYER: Delta=%f  Theta=%f  Y_p=%f" % (
        #    Delta, Theta, Y_p))

        return Delta, Y_p










if __name__ == "__main__":
    # tests
    print(__cf_s(1), __cf_s(4000), __cf_s(1e8))

    print(boundaryLayerThickness(
        mu=18.6e-6,
        d=0.01,
        L=0.1,
        u1=10, u2=10, u3=10,
        rho1=1.225, rho2=1.225, rho3=1.225
    ))