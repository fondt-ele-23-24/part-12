import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq

def dstate(X, t):
    # Parametri fissi
    rhoA = 0.0020 # Densita' aria
    rhoM = 3.32e3 # Densita' meteorite
    Cd = 0.47 # Coefficiente di trascinamento
    g = 9.7 # Accelerazione di gravita' nella mesosfera
    br = 0.75e-4 # Tasso di riduzione del raggio
    # "Spacchetto" lo stato
    x, v, r = X
    # Calcolo la massa della meteora
    m = rhoM * 4 / 3 * np.pi * r**3
    # Calcolo la superficie delle meteora
    S = np.pi * r**2
    # Calcolo la forza di gravità
    Fg = - m * g
    # Calcolo la forza di trascinamento
    Ft = -0.5 * rhoA * Cd * S * v * abs(v)
    # Calcolo le derivate
    dx = v
    dv = 1/m * (Fg+Ft)
    dr = -br * r * abs(v)
    return np.array([dx, dv, dr])

def simulate():
    v0 = -18e3
    d0 = 51413
    r0 = 0.11
    X0 = [d0, v0, r0]
    t = np.linspace(0, 10, 10 * 10000)
    X = odeint(dstate, X0, t)
    return X, t


def find_terminal_speed_aux(v):
    # Parametri fissi
    rhoA = 0.0020 # Densita' aria
    rhoM = 3.32e3 # Densita' meteorite
    Cd = 0.47 # Coefficiente di trascinamento
    g = 9.7 # Accelerazione di gravita' nella mesosfera
    br = 0.71e-5 # Tasso di riduzione del raggio
    r = 1e-6 # Raggio della meteora (fisso)
    # Calcolo la massa della meteora
    m = rhoM * 4 / 3 * np.pi * r**3
    # Calcolo la superficie delle meteora
    S = np.pi * r**2
    # Calcolo la forza di gravità
    Fg = - m * g
    # Calcolo la forza di trascinamento
    Ft = -0.5 * rhoA * Cd * S * v * abs(v)
    # Calcolo e restituisco la somma delle forze
    return Ft + Fg


def find_terminal_speed():
    a, b = -1000, 1
    vsol = brentq(find_terminal_speed_aux, a, b)
    return vsol


def find_vanishing_time_aux(t):
    X, tt = simulate()
    return np.interp(t, tt, X[:, 2]) - 1e-3


def find_vanishing_time():
    a, b = 0, 60 * 60
    tsol = brentq(find_vanishing_time_aux, a, b)
    return tsol