import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve

def dstate(X, t=0):
    # Parametri fissi
    a = 0.2 # Tasso di crescita della popolazione di conigli
    b = 0.0006 # Tasso di predazione
    d = 0.00004 # Bonus riproduttivo per predazione
    g = 0.1 # Tasso di calo della popolazione di volpi
    # "Spacchetto" lo stato
    r, f = X
    # Calcolo le derivate
    dr = a*r - b*r*f
    df = d*r*f - g*f
    # Restituisco le derivate
    return np.array([dr, df])


def simulate():
    x0 = [2000, 200]
    t = np.linspace(0, 12*10, 12*10*10)
    X = odeint(dstate, x0, t)
    return X, t


def find_eq():
    x0 = [2000, 200]
    Xsol = fsolve(dstate, x0)
    if max(np.abs(dstate(Xsol))) < 1e-6:
        return Xsol
    else:
        return None

    
def find_params_aux(X):
    # Parametri fissi
    a = 0.2 # Tasso di crescita della popolazione di conigli
    g = 0.1 # Tasso di calo della popolazione di volpi
    r = 2500 # Numero di conigli all'equilibrio
    f = 400 # Numero di volpi all'equilibrio
    # "Spacchetto" lo stato
    b, d = X
    # Calcolo le derivate
    dr = a*r - b*r*f
    df = d*r*f - g*f
    # Restituisco le derivate
    return np.array([dr, df])


def find_params():
    x0 = [0.0006, 0.00004]
    Xsol = fsolve(find_params_aux, x0)
    if max(np.abs(find_params_aux(Xsol))) < 1e-6:
        return Xsol
    else:
        return None