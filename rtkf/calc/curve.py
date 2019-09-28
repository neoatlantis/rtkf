#!/usr/bin/env python3

import numpy as np
from numpy.linalg import norm
from scipy.special import comb
from scipy import interpolate
from scipy.optimize import newton



CURVE_PRECISION = 50

separate = lambda points: (
    [e[0] for e in points], 
    [e[1] for e in points]
)

class DiscreteCurve:

    def __init__(self, points=None):
        self.__length = None
        if points:
            cx, cy = separate(points)
            self._curve = interpolate.interp1d(cx, cy)
            self._cxrange = cx[0], cx[-1]

    def __call__(self, t):
        if hasattr(self, "_curve"):
            cxrange = self._cxrange
            x = (cxrange[1] - cxrange[0]) * t + cxrange[0]
            y = self._curve(x)
            return np.array([x,y])
        else:
            raise NotImplementedError("Must override this.")

    def __iter__(self):
        for i in np.linspace(0, 1, CURVE_PRECISION):
            yield self.__call__(i)

    @property
    def length(self):
        if self.__length is not None:
            return self.__length
        self.__length = self.lengthOfT(1.0)
        return self.__length

    def lengthOfT(self, t):
        assert 0 <= t <= 1
        points = [self.__call__(i) for i in np.linspace(0, t, CURVE_PRECISION)]
        l = 0
        for i in range(0, len(points)-1):
            l += norm(points[i+1] - points[i])
        return l

    def tOfLength(self, l):
        if l < 0 or l > self.__length:
            raise Exception("Required length exceeds curve total length.")
        def solve(t):
            return l - self.lengthOfT(t)
        return newton(solve, 0.0)
        


class BezierCurve(DiscreteCurve):

    def __init__(self, points):
        DiscreteCurve.__init__(self)

        n = len(points)
        self.n = n

        pascalCoeffs = np.array([
            comb(n-1, i) for i in range(0, n)])

        self.xPoints = np.array([p[0] for p in points]) * pascalCoeffs
        self.yPoints = np.array([p[1] for p in points]) * pascalCoeffs

    def __call__(self, t):
        n = self.n
        poly = np.array([ 
            ((1-t)**(n-1-i)) * (t**i)
            for i in range(0, n)
        ])

        x = np.dot(self.xPoints, poly)
        y = np.dot(self.yPoints, poly)
        return np.array([x,y])



class MiddleOf2Curves(DiscreteCurve):
    """Find the middle of 2 curves. Each curve is a Y(X) function,
    but X ranges can differ. The idea is regarding both starting
    and ending point as correlated and find a Y(X) middle curve
    between, which is parameterized by t in [0,1]."""

    def __init__(self, curve1, curve2, func=lambda a,b: (a+b)/2):
        DiscreteCurve.__init__(self)
        self.curve1 = DiscreteCurve(curve1)
        self.curve2 = DiscreteCurve(curve2)
        self.func = func

    def __call__(self, t):
        p1 = self.curve1(t)
        p2 = self.curve2(t)
        return self.func(p1, p2)



if __name__ == "__main__":

    b1 = BezierCurve([
        [0, 0.02],
        [0.00106, 0.0149],
        [0.002726, 0.00727],
        [0.00554, 0.00658],
        [0.014, 0.0045],
    ])

    b2 = BezierCurve([
        [0.006955, 0.02],
        [0.0076098, 0.0184],
        [0.008592, 0.016],
        [0.009409, 0.016],
        [0.01186, 0.016],
    ])

    bm = MiddleOf2Curves(b1, b2)

    print("bm_length", bm.length)
    print("b1_length", b1.length)
    print("b2_length", b2.length)
    exit()

    a = list(b1)
    b = list(b2)
    m = list(bm)

    for i in range(0, len(m)):
        print(*a[i], *b[i], *m[i]) 

