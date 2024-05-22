import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq

def force_balance(r):
    # Parametri fissi
    rhof = 979 # Densita' del galleggiante
    rhow = 1000 # Densita' dell'acqua
    Cd = 0.62 # Coefficiente di trascinamento
    ms = 8000 # Massa equivalente della cabina in acqua
    mb = 9000 - 1143 # Massa equivalente della zavorra in acqua
    g = 9.81 # Accelerazione di gravita' (costante, per semplicita')
    vf = -0.9 # Velocita' di discesa desiderata
    # Calcolo il volume del galleggiante
    Vf = 4/3 * np.pi * r**3
    # Calcolo la superficie del galleggiante
    Sf = np.pi * r**2
    # Calcolo la massa del galleggiante
    mf = rhof * Vf
    # Calcolo la forza di gravità
    Fg = - g * (ms + mb + mf)
    # Calcolo la forza di galleggiamento
    Fb = g * rhow * Vf
    # Calcolo la forza di trascinamento
    Ft = -0.5 * rhow * Cd * Sf * vf * abs(vf)
    # Calcolo e restituisco la somma delle forze
    return Fg + Fb + Ft


def find_r():
    a, b = 1, 10
    rsol = brentq(force_balance, a, b)
    return rsol


def dstate(X, t):
    # Parametri fissi
    rhof = 979 # Densita' del galleggiante
    rhow = 1000 # Densita' dell'acqua
    Cd = 0.62 # Coefficiente di trascinamento
    ms = 8000 # Massa equivalente della cabina in acqua
    mb = 9000 - 1143 # Massa equivalente della zavorra in acqua
    g = 9.81 # Accelerazione di gravita' (costante, per semplicita')
    r = 5.3 # Raggio del galleggiante
    # "Spacchetto" lo stato
    x, v = X
    # Calcolo il volume del galleggiante
    Vf = 4/3 * np.pi * r**3
    # Calcolo la superficie del galleggiante
    Sf = np.pi * r**2
    # Calcolo la massa del galleggiante
    mf = rhof * Vf
    # Calcolo la forza di gravità
    Fg = - g * (ms + mb + mf)
    # Calcolo la forza di galleggiamento
    Fb = g * rhow * Vf
    # Calcolo la forza di trascinamento
    Ft = -0.5 * rhow * Cd * Sf * v * abs(v)
    # Calcolo le derivate
    dx = v
    dv = 1 / (ms + mb + mf) * (Fg + Fb + Ft)
    return np.array([dx, dv])

    
def simulate():
    v0 = 0
    x0 = 0
    X0 = [x0, v0]
    t = np.linspace(0, 300, 300 * 10)
    X = odeint(dstate, X0, t)
    return X, t



def find_t250_aux(t):
    X, tt = simulate()
    return np.interp(t, tt, X[:, 0]) + 250


def find_t250():
    a, b = 0, 60 * 60
    tsol = brentq(find_t250_aux, a, b)
    return tsol