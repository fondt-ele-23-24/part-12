import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq


def dstate(X, t, G=6.67408e-11, ME=5.972e24, MM=7.34767309e22, MS=800, D=384400e3):
    # "spacchetto" lo stato
    x, v = X

    # Calcolo le forze
    rse = x
    Fse = -G * MS * ME / (rse * abs(rse))
    rsm = x - D
    Fsm = -G * MS * MM / (rsm * abs(rsm))
    # Calcolo le derivate
    dx = v
    dv = 1/MS * (Fse + Fsm)
    # Restituisco il risultato
    return np.array([dx, dv])


def simulate(v0, rE=6371e3, t=np.linspace(0, 3*3600, 3*3600*100)):
    x0 = [rE, v0]
    X = odeint(dstate, x0, t)
    return X, t


def find_balance_point_aux(x, G=6.67408e-11, ME=5.972e24, MM=7.34767309e22, MS=800, D=384400e3):
    # Calcolo le forze
    rse = x
    Fse = -G * MS * ME / (rse * abs(rse))
    rsm = x - D
    Fsm = -G * MS * MM / (rsm * abs(rsm))

    # Restituisco il risultato
    return Fse + Fsm


def find_balance_point(a=1e6, b=384e6):
    x_sol = brentq(find_balance_point_aux, a=a, b=b)
    return x_sol


def distance_in_1h(v0):
    X, t= simulate(v0)
    res = np.interp(3600, t, X[:, 0])
    return res


def distances_for_v0_range(v0_range):
    return [distance_in_1h(v0) for v0 in v0_range]
