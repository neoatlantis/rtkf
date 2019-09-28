#!/usr/bin/env python3

import numpy as np
from numpy import tan, sin, cos, pi, arctan
from numpy.linalg import norm
from scipy.optimize import newton, fminbound, fsolve
from scipy import integrate, interpolate


def findThroatFromCurve3d(delta, curve3d):
    n = len(curve3d)
    assert n > 50

    cos_delta = cos(delta)
    sin_delta = sin(delta)

    getXYZ1 = lambda i: np.array([curve3d[i].x, curve3d[i].y, curve3d[i].z])
    def getXYZ2(i):
        x1, y1, z1 = getXYZ1(i)
        x2 = cos_delta * x1 - sin_delta * y1
        y2 = sin_delta * x1 + cos_delta * y1
        return np.array([x2, y2, z1])

    minimalFinder = []

    for i1 in range(0, n-1):
        v_t1 = getXYZ1(i1+1) - getXYZ1(i1)

        lastDot1, lastDot2 = None, None
        for i2 in range(0, n-1):
            v_t2 = getXYZ2(i2+1) - getXYZ2(i2)
            v_12 = getXYZ2(i2) - getXYZ1(i1)

            dot1 = np.dot(v_t1, v_12)
            dot2 = np.dot(v_12, v_t2)

            if lastDot1 is not None and lastDot2 is not None:
                if dot1 * lastDot1 <=0 and dot2 * lastDot2 <= 0:
                    # found
                    return i1, i2, norm(v_12)
            lastDot1 = dot1
            lastDot2 = dot2
            minimalFinder.append((i1, i2, norm(v_12)))

    minimalFinder.sort(key=lambda e: e[2])
    return minimalFinder[0]
