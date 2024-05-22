import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve


class Dstate:
    def __init__(self, rho):
        self.rho = rho
    
    def __call__(self, X, t):
        # Parametri del problema
        sigma = 10
        beta = 8/3
        
        # "Spacchetto" lo stato
        x, y, z = X
        # Calcolo le derivate
        dx = sigma * (y - x)
        dy = x * (self.rho - z) - y
        dz = x * y - beta * z
        return np.array([dx, dy, dz])

def simulate(rho):
    X0 = [8, 8, 25]
    t = np.linspace(0, 60, 60 * 10)
    f = Dstate(rho)
    X = odeint(f, X0, t)
    return X, t


def find_eq_aux(X):
    # Parametri del problema
    sigma = 10
    beta = 8/3
    rho = 28

    # "Spacchetto" lo stato
    x, y, z = X
    # Calcolo le derivate
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return np.array([dx, dy, dz])

    
def find_eq():
    X0 = [8, 8, 25]
    Xsol = fsolve(find_eq_aux, X0)
    if np.max(np.abs(find_eq_aux(Xsol))) <= 1e-6:
        return Xsol
    else:
        return None
    

def find_rho_aux(Xt):
    # Spacchetto lo stato
    x, y, rho = Xt
    f = Dstate(rho)
    return f([x, y, 25], None)


def find_rho():
    X0 = [8, 8, 25]
    Xsol = fsolve(find_rho_aux, X0)
    if np.max(np.abs(find_rho_aux(Xsol))) <= 1e-6:
        return Xsol
    else:
        return None

    
