#!/usr/bin/env python3

from math import pi, log, sqrt

def diskFrictionEnergy(C_M, rho, omega, r):
    return 0.5 * C_M * rho * omega**3 * r**5 # for one side disc

def calcRe(rho, omega, r, mu):
    return rho * omega * r**2 / mu

calcC_M1 = lambda Re, delta2r: 2 * pi / (delta2r * Re)
calcC_M2 = lambda Re, delta2r: 3.7 * delta2r**0.1 / sqrt(Re)
calcC_M3 = lambda Re, delta2r: 0.08 / (delta2r**(1/6) * Re**(1/4))
calcC_M4 = lambda Re, delta2r: 0.102 * delta2r**0.1 / Re**0.2

def calcC_Mr(e2r, delta2r):
    if e2r > 0:
        return (1/(3.8 * log(1/e2r, 10) - 2.4 * delta2r**0.25))**2
    return 0

calcRe_r = lambda e2r: 1100/e2r - 6e6


def calcC_M(rho, omega, mu, e, r, delta):
    delta2r = delta / r
    e2r = e / r
    Re = calcRe(rho, omega, r, mu)    

    C_Ms = sum([
        calc(Re, delta2r)
        for calc in [calcC_M1, calcC_M2, calcC_M3, calcC_M4] 
    ])
    C_M3 = calcC_M3(Re, delta2r)

    C_Mr = calcC_Mr(e2r, delta2r)
    C_M = C_Ms

    if e > 0:
        Re_s = 1100 * e2r**(-0.4) / sqrt(C_M3)
        Re_r = calcRe_r(e2r)

        C_M += (C_Mr - C_Ms) * log(Re/Re_s) / log(Re_r/Re_s)

    
    return C_M


def DiskFrictionLoss(m_dot, rho, omega, r, delta, mu, e=0):
    if delta == 0:
        return 0
        
    C_M = calcC_M(rho, omega, mu, e, r, delta)
    H_DF = diskFrictionEnergy(C_M, rho, omega, r)

    return H_DF / m_dot


if __name__ == "__main__":
    print(DiskFrictionLoss(
        m_dot=0.3,
        rho=2,
        omega=1800 * 2 * pi,
        r=0.05,
        delta=0.005,
        mu=18.5e-6
    ))




