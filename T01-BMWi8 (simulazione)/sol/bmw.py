import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq

class Dstate:
    def __init__(self, Cd):
        self.Cd = Cd
    
    def __call__(self, X, t):
        # Parametri fissi
        rho = 1.25 # Densita' dell'aria
        A = 2.5 * 1.2 # Superficie della seziojne
        m = 1539 # Massa dell'auto
        F = 10000 # Forza di accelerazione
        # "Spacchetto" lo stato
        x, v = X
        # Calcolo le forze
        Ft = -0.5 * rho * A * self.Cd * v * np.abs(v)
        # Calcolo le derivate
        dx = v
        dv = 1/m * (F + Ft)
        return np.array([dx, dv])


def simulate(Cd):
    x0 = [0, 0]
    t = np.linspace(0, 60, 60000)
    f = Dstate(Cd)
    X = odeint(f, x0, t)
    return X, t


def acceleration(Cd):
    f = Dstate(Cd)
    return 


def find_terminal_speed_aux(v):
    # Parametri fissi
    rho = 1.25 # Densita' dell'aria
    A = 2.5 * 1.2 # Superficie della seziojne
    m = 1539 # Massa dell'auto
    F = 10000 # Forza di accelerazione
    Cd = 0.82
    
    # Calcolo la forza di trascinamento
    Ft = -0.5 * rho * A * Cd * v * np.abs(v)
    # Restituisco il risultato
    return F + Ft
    

def find_terminal_speed():
    a, b = 0, 100
    v_sol = brentq(find_terminal_speed_aux, a, b)
    return v_sol


def speed_in_5_seconds(Cd):
    X, t = simulate(Cd)
    res = np.interp(5, t, X[:, 1])
    return res


def find_Cd_aux(Cd):
    return speed_in_5_seconds(Cd) - 31


def find_Cd():
    a, b = 0.2, 1.0
    Cd_sol = brentq(find_Cd_aux, a=a, b=b)
    return Cd_sol