import numpy as np
from base import util
from scipy.integrate import odeint
from scipy.optimize import fsolve

class Dstate:
    def __init__(self, C1, C2, C3, R12, R23, R31):
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.R12 = R12
        self.R23 = R23
        self.R31 = R31

    def __call__(self, X, t):
        # "spacchetto" lo stato
        P1, P2, P3 = X
        # Calcolo i flussi
        q12 = 1/self.R12 * (P1 - P2)
        q23 = 1/self.R23 * (P2 - P3)
        q31 = 1/self.R31 * (P3 - P1)
        # Calcolo le derivate
        dP1 = 1/self.C1 * (q31 - q12)
        dP2 = 1/self.C2 * (q12 - q23)
        dP3 = 1/self.C3 * (q23 - q31)
        # Restituisco il risultato
        return np.array([dP1, dP2, dP3])

    
def simulate(f, X0, t):
    X = odeint(f, X0, t)
    return X, t


def find_equilibrium_aux(X):
    # Dati intermedi
    g = 9.81
    S1 = 1 # Superfici
    S2 = 1
    S3 = 1
    h1 = 3 # Livelli dei serbatoi
    h2 = 2
    h3 = 1
    qmax12 = 0.0002 # Portata per unita' di pressione
    qmax23 = 0.0007
    qmax31 = 0.0008
    # Capacita'
    C1 = S1/g
    C2 = S2/g
    C3 = S3/g
    # Resistenze
    R12 = 1/qmax12
    R23 = 1/qmax23
    R31 = 1/qmax31
    # Pressioni iniziali
    P1_0 = h1 * g
    P2_0 = h2 * g
    P3_0 = h3 * g
    
    # "spacchetto" lo stato
    P1, P2, P3 = X
    # Calcolo i flussi
    q12 = 1/R12 * (P1 - P2)
    q23 = 1/R23 * (P2 - P3)
    q31 = 1/R31 * (P3 - P1)
    # Calcolo le derivate
    dP1 = 1/C1 * (q31 - q12)
    dP2 = 1/C2 * (q12 - q23)
    # Imposto la condizione sulle pressioni
    p_tot = P1+P2+P3-1
    # Restituisco il risultato
    return np.array([dP1, dP2, p_tot])

def find_equilibrium():
    x0 = [5,2,12] # Pressioni iniziali
    x_sol = fsolve(find_equilibrium_aux, x0)
    return x_sol

def find_resistances_aux(x):
    # Pressioni
    P1,P2,P3 = [5,4,6]
    
    # "spacchetto" le resistenze 
    R12,R23,R31 = x
    # Calcolo i flussi
    q12 = 1/R12 * (P1 - P2)
    q23 = 1/R23 * (P2 - P3)
    q31 = 1/R31 * (P3 - P1)
    # Imposto le condizioni sui flussi
    c1 = q12-1
    c2 = q12-1
    c3 = q12-1
    
    return np.array([c1, c2, c3])

def find_resistances():
    x0 = [4,5,7] # Resistenze iniziali
    R12,R23,R31 = fsolve(find_resistances_aux, x0)
    return R12,R23,R31