#!/usr/bin/env python3

import numpy
from numpy.linalg import inv

def newton_multi(func, x0, etol=1e-3, dx=1, itermax=1000, stepdiv=1):
    f0 = func(x0)
    n = len(f0)
    if n != len(x0):
        raise Exception("func() must return a vector of same length as x0")
    def can_exit(v):
        for e in v:
            if abs(e) > etol: return False
        return True
    x1 = lambda i, dx: numpy.array([
        0 if i != j else dx
        for j in range(0, n)
    ]) + x0
    x = x0
    itercount=0
    while not can_exit(f0):
        if itercount > itermax:
            raise Exception("Convergence not found by newton method.")
        jacob = []
        for i in range(0, n):
            f1_i = func(x1(i, dx))
            df_i = f1_i - f0
            jacob.append(df_i / dx)
        try:
            invjacob = inv(jacob)
        except Exception as e:
            print(jacob)
            raise e
        ndx = numpy.dot(-f0, invjacob)
        x += (ndx / stepdiv)
        f0 = func(x)
        itercount+=1
    return x

if __name__ == "__main__":
    from numpy import sin, cos, exp
    def func(x):
        x1, x2, x3 = x
        return numpy.array([
            3*x1 - cos(x2*x3) - 1.5,
            4*x1**2 - 625*x2**2 + 2*x3 - 1,
            20*x3 + exp(-x1*x2) + 9
        ])
    print(newton_multi(func, [0,0,0]))


