import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import brentq

def find_curve_params(t1, t0=0, v0=5.5, v1=0):
    A = [[t0**3,     t0**2,    t0, t0**0],
         [t1**3,     t1**2,    t1, t1**0],
         [3*t0**2, 2*t0**1, t0**0,     0],
         [3*t1**2, 2*t1**1, t1**0,     0]]
    b = [v0, v1, 0, 0]
    (a3, a2, a1, a0) = np.linalg.solve(A, b)
    return (a3, a2, a1, a0)


def curve(a3, a2, a1, a0, t):
    return a3*t**3 + a2*t**2 + a1*t + a0


def find_curve_params2_aux(X, t0=0, v0=5.5, v1=0, a3=0.171875):
    # Spacchetto il vettore di ingresso
    a2, a1, a0, t1 = X

    # Calcolo le espressioni da azzerare
    c1 = a3 * t0**3 + a2 * t0**2 + a1 * t0 + a0 - v0
    c2 = a3 * t1**3 + a2 * t1**2 + a1 * t1 + a0 - v1
    c3 = 3*a3 * t0**2 + 2 * a2 * t0 + a1
    c4 = 3*a3 * t1**2 + 2 * a2 * t1 + a1

    return np.array([c1, c2, c3, c4])


def find_curve_params2(X0=[-1, 0.2, 6.0, 4.5]):
    (a2_sol, a1_sol, a0_sol, t1_sol) = fsolve(find_curve_params2_aux, X0)
    return a2_sol, a1_sol, a0_sol, t1_sol


def braking_distance(t1):
    a3, a2, a1, a0 = find_curve_params(t1)
    bd = 1/4 * a3 * t1**4 + 1/3 * a2 * t1**3 + 1/2 * a1 * t1**2 + a0 * t1
    return bd


def find_t1_aux(t1):
    return braking_distance(t1) - 10


def find_t1(a=3, b=12):
    t1_sol = brentq(find_t1_aux, a, b)
    return t1_sol
