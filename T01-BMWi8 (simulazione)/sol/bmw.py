import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq

class Dstate:
    def __init__(self, Cd, rho=1.25, A=2.5*1.2, m=1539, F=10000):
        self.Cd = Cd
        self.rho = rho
        self.A = A
        self.m = m
        self.F = F

    def __call__(self, X, t):
        # "Spacchetto" lo stato
        x, v = X
        # Calcolo le forze
        Ft = -0.5 * self.rho * self.A * self.Cd * v * np.abs(v)
        # Calcolo le derivate
        dx = v
        dv = 1/self.m * (self.F + Ft)
        return np.array([dx, dv])


def simulate(Cd, x0=[0, 0], t=np.linspace(0, 60, 60000)):
    f = Dstate(Cd)
    X = odeint(f, x0, t)
    return X, t


def find_terminal_speed_aux(v, Cd=0.82, rho=1.25, A=2.5*1.2, m=1539, F=10000):
    # Calcolo la forza di trascinamento
    Ft = -0.5 * rho * A * Cd * v * np.abs(v)
    # Restituisco il risultato
    return F + Ft


def find_terminal_speed(a=0, b=100):
    v_sol = brentq(find_terminal_speed_aux, a, b)
    return v_sol


def speed_in_5_seconds(Cd):
    X, t = simulate(Cd)
    res = np.interp(5, t, X[:, 1])
    return res


def find_Cd_aux(Cd):
    return speed_in_5_seconds(Cd) - 31


def find_Cd(a=0.2, b=1.0):
    Cd_sol = brentq(find_Cd_aux, a=a, b=b)
    return Cd_sol
