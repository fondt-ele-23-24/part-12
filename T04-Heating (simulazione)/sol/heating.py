import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve


def dstate(X, t):
    # Parametri fissi
    g = 9.81
    vA = 62 # Volume dell'aria
    vW = 0.25 * 16 * 2.7 # Volume dei muri
    mA = 1.225 * vA / g # Massa dell'aria
    mW = 1050 * vW / g # Massa dei muri
    Ca = 1005 * mA # Capacita' termica dell'aria
    Cw = 1000 * mW # Capacita' termica dei muri
    Rca = 0.35 # Resisitivita' termica convettore-aria
    Raw = 0.5 # Resistivita' termica aria-muro
    Rwo = 3.0 # Resistivita' termica muro-esterno
    Tc = 23 # Temperatura del convettore
    To = 15 # Temperatura esterna

    # "spacchetto" lo stato
    Ta, Tw = X
    # Calcolo le correnti
    wca = 1/Rca * (Tc - Ta)
    waw = 1/Raw * (Ta - Tw)
    wwo = 1/Rwo * (Tw - To)
    # Calcolo le derivate
    dTa = 1/Ca * (wca - waw)
    dTw = 1/Cw * (waw - wwo)
    # Restituisco il risultato
    return np.array([dTa, dTw])

    
def simulate():
    X0 = [19.5, 19.5]
    t = np.linspace(0, 3600 * 2, 3600 * 2)
    X = odeint(dstate, X0, t)
    return X, t


def temp_in_20min(X, t):
    Ta20 = np.interp(60*20, t, X[:, 0])
    Tw20 = np.interp(60*20, t, X[:, 1])
    return Ta20, Tw20

def find_equilibrium_aux(X):
    # Parametri fissi
    g = 9.81
    vA = 62 # Volume dell'aria
    vW = 0.25 * 16 * 2.7 # Volume dei muri
    mA = 1.225 * vA / g # Massa dell'aria
    mW = 1050 * vW / g # Massa dei muri
    Ca = 1005 * mA # Capacita' termica dell'aria
    Cw = 1000 * mW # Capacita' termica dei muri
    Rca = 0.35 # Resisitivita' termica convettore-aria
    Raw = 0.5 # Resistivita' termica aria-muro
    Rwo = 3.0 # Resistivita' termica muro-esterno
    Tc = 23 # Temperatura del convettore
    To = 15 # Temperatura esterna

    # "spacchetto" lo stato
    Ta, Tw = X
    # Calcolo le correnti
    wca = 1/Rca * (Tc - Ta)
    waw = 1/Raw * (Ta - Tw)
    wwo = 1/Rwo * (Tw - To)
    # Calcolo le derivate
    dTa = 1/Ca * (wca - waw)
    dTw = 1/Cw * (waw - wwo)

    return np.array([dTa, dTw])
    

def find_equilibrium():
    x0 = [15,25] # Temperature iniziali
    x_sol = fsolve(find_equilibrium_aux, x0)
    return x_sol
    