import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve

class Dstate:
    def __init__(self, a, b):
        self.a = a
        self.b = b
    
    def __call__(self, X, t):
        # "Spacchetto" lo stato
        x, y = X
        # Calcolo le derivate
        dx = self.a + x**2 * y - self.b * x - x
        dy = self.b * x - x**2 * y
        return np.array([dx, dy])


def simulate(a, b):
    x0 = [1, 1]
    t = np.linspace(0, 120, 10 * 120)
    f = Dstate(a, b)
    X = odeint(f, x0, t)
    return X, t


def find_eq_aux(X):
    a, b = 0.4, 1.3
    f = Dstate(a=a, b=b)
    return f(X, t=0)


def find_eq():
    X0 = [1, 1]
    Xsol = fsolve(find_eq_aux, X0)
    if np.max(np.abs(find_eq_aux(Xsol))) <= 1e-6:
        return Xsol
    else:
        return None


def control_eq(x, y):
    A = np.array([[1, -x],
                  [0,  x]])
    B = np.array([x - x**2 * y, x**2 * y])
    sol = np.linalg.solve(A, B)
    a, b = sol
    return a, b

