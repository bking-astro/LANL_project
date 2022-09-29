'''
We want a function in script form (real_TOVsolver.py) we can call and get a mass-radius curve by solving the
Tolman-Oppenheimer-Volkoff (TOV) equations. This function will be given baryon density (0.16 fm-3),
pressure (MeV/fm3), and energy density (MeV/fm3) and return the neutron star's mass (Mo), radius (km),
and the tidal deformability (unitless). We solver the TOV equations for a range of central densities/pressures,
so we also pass the function a maximum pressure (MeV/fm3) up to which we integrate.
'''

# some libraries we will need for these functions
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import scipy.interpolate
import pandas as pd

# some constants to make things neater
pi = np.pi
MeV_to_km = 1.323e-6
Mo_to_km = 1.4766
p_c = 4.42e-5 # km^-2


# first we define some auxiliary functions to get the tidal deformability

def F_ode(p, y, m, EOS_e):
    a = np.sqrt(y) / (np.sqrt(y) - 2 * m)
    b = 4 * pi * (EOS_e(p) - p)
    return a * (1 - b)


def yQ_ode(p, y, m, EOS_e):
    a = 4 * pi * y * np.sqrt(y)
    b = np.sqrt(y) - 2 * m

    one = (a / b) * (5 * EOS_e(p) + 9 * p + EOS_e.derivative(nu=1)(p) * (EOS_e(p) + p))
    two = -6 * np.sqrt(y) / b
    three = -4 * (m / b + a * p / b) ** 2

    return one + two + three


# then functions to calculate lambda

def k_2(C, x_R):
    num = ((1 - 2 * C) ** 2) * (2 - x_R + 2 * C * (x_R - 1))

    a = 6 - 3 * x_R + 3 * C * (5 * x_R - 8)
    b = 13 - 11 * x_R + C * (3 * x_R - 2) + (2 * C ** 2) * (1 + x_R)
    c = 2 - x_R + 2 * C * (x_R - 1)

    den = 2 * C * a + 4 * (C ** 3) * b + 3 * ((1 - 2 * C) ** 2) * c * np.log(1 - 2 * C)

    return (8 * (C ** 5) * num) / (5 * den)


def big_lamb(k2, R, M):
    """
    k2: dimensionless tidal Love number
    R: radius of NS in km
    M: mass of NS in Mo (solar masses)
    """
    return (2 / 3) * k2 * ((R / (Mo_to_km * M)) ** 5)


# finally we define the TOV equations in the code

def TOV(y, p, EOS_e):
    """
    y: numpy array containing the vector we are solving (y=r^2 [km2], M [km2], x[unitless])
    p: scalar pressure at which we are evaluating the equations
    EOS_e: interpolation of the energy density function
    """

    dy = np.zeros(y.shape)

    num0 = -2 * y[0] * (np.sqrt(y[0]) - 2 * y[1])
    den0 = (EOS_e(p) + p) * (y[1] + 4 * pi * p * y[0] ** (3 / 2))

    dy[0] = num0 / den0

    num1 = -4 * np.pi * EOS_e(p) * (y[0] ** (3 / 2)) * (np.sqrt(y[0]) - 2 * y[1])
    den1 = (EOS_e(p) + p) * (y[1] + 4 * pi * p * y[0] ** (3 / 2))

    dy[1] = num1 / den1

    num2 = (y[2] ** 2 + y[2] * F_ode(p, y[0], y[1], EOS_e)
            + yQ_ode(p, y[0], y[1], EOS_e)) * (np.sqrt(y[0]) - 2 * y[1])
    den2 = (EOS_e(p) + p) * (y[1] + 4 * pi * p * y[0] ** (3 / 2))

    dy[2] = num2 / den2

    return dy


def solve(EOS_table, max_pressure):
    """
    EOS_table: numpy array containing the baryon density (0.16 fm-3), pressure (MeV/fm3), and energy density (MeV/fm3).
    max_pressure: Scalar value for maximum central pressure (MeV/fm3) used when solving the TOV equations.

    MRL_table: numpy array containing mass (Mo), radius (km), Lambda (unitless)
    """

    # get pressure and energy density in numpy array
    p = MeV_to_km * EOS_table[:,1]
    e = MeV_to_km * EOS_table[:,2]

    # make interpolation for the energy density
    EOS_e = scipy.interpolate.CubicSpline(p, e)

    # define the central pressures over which we will solve TOV eqs
    size = 100
    start = np.log(p_c * 1e-1)  # start/stop for log spacing
    stop = np.log(MeV_to_km * max_pressure)
    central_ps = np.e ** (np.linspace(start, stop, size))  # central pressures in km^-2

    MRL_table = np.zeros((central_ps.shape[0], 3))  # (mass_Mo, radius_km, Lambda)

    i = 0
    for pc in central_ps:
        t0 = pc
        y0 = np.array([1e-7, 1e-10, 2])
        tf = 1e-10 * t0

        leng = int(1e3)
        t = np.linspace(t0, tf, leng)

        sol = scipy.integrate.odeint(TOV, y0, t, args=(EOS_e,))

        MRL_table[i, 0] = max(sol[:, 1] / Mo_to_km)
        MRL_table[i, 1] = max(np.sqrt(sol[:, 0]))

        x_R = sol[-1, 2]
        C = Mo_to_km * MRL_table[i, 0] / MRL_table[i, 1]
        k2 = k_2(C, x_R)

        MRL_table[i, 2] = big_lamb(k2, MRL_table[i, 1], MRL_table[i, 0])

        i += 1

    return MRL_table