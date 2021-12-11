from scipy.integrate import odeint
import numpy as np
from scipy.integrate import solve_ivp


# dgl dMquer/dPquer in Abhängigkeit von Masse m, Druck p, Energiedichte epsilon (abhängig von p) und y = länge ** 2
def dmdp(mquer, pquer, epsilonquer, y):
    return -4 * np.pi * epsilonquer(pquer) * abs(y) ** (3 / 2) * (abs(y) ** (1 / 2) - 2 * mquer) / (
            (epsilonquer(pquer) + pquer) * (mquer + 4 * np.pi * abs(y) ** (3 / 2) * pquer))


# dgl dy/dPquer in Abhängigkeit von y = länge ** 2, Druck p, Energiedichte epsilon (abhängig von p) und Masse m
def dydp(y, pquer, epsilonquer, mquer):
    return -2 * abs(y) * (abs(y) ** (1 / 2) - 2 * mquer) / ((epsilonquer(pquer) + pquer) * (mquer + 4 * np.pi * abs(y) ** (3 / 2) *
                                                                                  pquer))


# dgl als Vektor geschrieben. m in der 0 und y in der 1
def ddp(mquer_y, pquer, epsilonquer):
    return np.array([dmdp(mquer_y[0], pquer, epsilonquer, mquer_y[1]), dydp(mquer_y[1], pquer, epsilonquer, mquer_y[0])])


def ddp_different_order(pquer, mquer_y, epsilonquer):
    return np.array([dmdp(mquer_y[0], pquer, epsilonquer, mquer_y[1]), dydp(mquer_y[1], pquer, epsilonquer, mquer_y[0])])


# löst die DGL in Abhängigkeit vom zentraldruck p_center, energiedichte epsilon und der anzahl an stützstellen n.
def mass_y_p(p_center, epsilon, n):
    p = []
    for i in np.arange(0, n, 1):
        p = np.append(p, p_center/2**i)
    if p[-1] > 1e-16:
        p = np.append(p, 1e-16)
    solved = odeint(ddp, y0=np.array([1e-20, 1e-20]), t=p, args=(epsilon,), tcrit=[0])
    return solved


def mass_y_p_ivp(p_center, epsilon):
    solved = solve_ivp(ddp_different_order, t_span=(p_center, 1e-16), y0=np.array([1e-40, 1e-40]), args=(epsilon,), method='RK45', rtol=1e-7)
    return solved


# Gibt M und R in Sonnenmassen bzw. km aus. Abhängig vom Zentraldruck (km^-2), der Energiedichte und Stützstellen.
def M_R_relation(p_center, epsilon, n):
    if p_center == 0:
        return np.array([0, 0])
    mass_y = mass_y_p(p_center, epsilon, n)
    m_y = mass_y[-1]
    m = m_y[0]/1.47661
    r = np.sqrt(m_y[1])
    return r, m

def M_R_relation_ivp(p_center, epsilon, n):
    if p_center == 0:
        return np.array([0, 0])
    mass_y = mass_y_p_ivp(p_center, epsilon)
    m = mass_y.y[0][-1] / 1.47661
    r = np.sqrt(mass_y.y[1][-1])
    return r, m
