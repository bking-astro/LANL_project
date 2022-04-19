# import some libraries

import scipy.integrate
import scipy.interpolate
import numpy as np
import scipy.stats as sts

# define some constants
n0 = 0.16 #MeV/fm^3
pi = np.pi
MeV_to_km = 1.323e-6
max_p = 350 #MeV/fm^3

# define some functions to sample EOS space

def sample_ns(num_points, n_start, max_n):
    """
    Function to create an array of randomly sampled density points in range (n_start, max_n)

    Inputs:
    num_points: int defining number of samples taken
    n_start: first point in returned density array
    max_n: last point in array

    Outputs:
    ns: np array of length = num_points + 2 of density points ordered from least to greatest
    """

    # sample random points in density
    epsilon = 1e-17
    loc_n = n_start + epsilon  # to guarentee we don't randomly pull n_start
    scale_n = max_n - n_start - epsilon
    sample_n = sts.uniform.rvs(loc=loc_n, scale=scale_n, size=num_points - 1)
    # order these points
    n_sort = np.sort(sample_n)

    # construct the arrays
    n_0 = n_start * np.ones(1)
    n = np.append(n_0, n_sort)
    ns = np.append(n, max_n)

    return ns


def sample_cs(num_points, cs_start):
    """
    Function to create an array of randomly sampled speed of sound points in range (0, 1)

    Inputs:
    num_points: int defining number of samples taken
    cs_start: first point in returned array to match EOS we are extending

    Outputs:
    cs: np array of length = num_points + 1 of speed of sound points
    """

    # sample speeds of sound
    sample_cs = sts.uniform.rvs(size=num_points)

    c_0 = cs_start * np.ones(1)
    cs = np.append(c_0, sample_cs)

    return cs


def sample_polytrop(num_points, ns, p_start, max_gamma):
    """
    Function to randomly sample parameters for the polytropic segments extension

    Inputs:
    num_points: int for number of samples taken
    ns: array of randomly sampled and ordered density points from sample_ns
    p_start: int defining the pressure at which the extension starts
    max_gamma: int defining the maximum value the parameter gamma can take

    Outputs:
    gammas: array of length num_points + 1 containing polytropic segment parameters gamma_i
    Ks: array of length num_points + 1 containing polytropic segment parameters K_i
    """

    n_0 = ns[0] * np.ones(1)

    # sample random gammas
    gammas = sts.uniform.rvs(scale=max_gamma, size=num_points + 1)
#     require the first gamma is large enough
    if gammas[0] < 1.5:
        gammas = sts.uniform.rvs(scale=max_gamma, size=num_points + 1)

    # initialize K array
    K0 = p_start * (n_0 ** (-gammas[0]))
    K = np.zeros(num_points)
    Ks = np.append(K0, K)

    for i in range(num_points):
        Ks[i + 1] = Ks[i] * (ns[i + 1] ** (gammas[i] - gammas[i + 1]))

    return gammas, Ks

# functions to extend the EOS based on speed of sound, linear pressure segments, or polytropic segments

def extend_Plin(n_step, ns, starts, cs2):
    """
    Function to extend the EOS with linear segments in pressure by using constant segments in speed of sound

    Inputs:
    n_step: int defining step size in density
    ns: array containing the densities where segments begin/end
    starts: array containing the starting values for the EOS (density, pressure, energy)
    cs2: array ocontaining randomly sampled speed of sound values

    Outputs:
    EOS_ex: array of shape (size, 3) containing the values for EOS extension
    """

    size = int((ns[-1] - starts[0]) / n_step)

    # initialize array
    EOS_ex = np.zeros((size, 3))
    # set starting values at n = 2n0
    EOS_ex[0, 0] = starts[0]
    EOS_ex[0, 1] = starts[1]
    EOS_ex[0, 2] = starts[2]

    i = 0

    for k in range(size - 1):
        # n_i+1
        EOS_ex[k + 1, 0] = EOS_ex[k, 0] + n_step
        # p_i+1
        if ns[i] <= EOS_ex[k, 0] < ns[i + 1]:
            EOS_ex[k + 1, 1] = EOS_ex[k, 1] + n_step * cs2[i] * ((EOS_ex[k, 1] + EOS_ex[k, 2]) / EOS_ex[k, 0])
        else:
            i += 1
            EOS_ex[k + 1, 1] = EOS_ex[k, 1] + n_step * cs2[i] * ((EOS_ex[k, 1] + EOS_ex[k, 2]) / EOS_ex[k, 0])
        # e_i+1
        EOS_ex[k + 1, 2] = EOS_ex[k, 2] + n_step * ((EOS_ex[k, 1] + EOS_ex[k, 2]) / EOS_ex[k, 0])

    return EOS_ex


def extend_cs(n_step, ns, starts, cs2_func):
    """
    Function to extend the EOS by using linear segments in speed of sound

    Inputs:
    n_step: int defining step size in density
    ns: array containing the densities where segments begin/end
    starts: array containing the starting values for the EOS (density, pressure, energy)
    cs2_func: linear interpolation of randomly sampled speed of sound values

    Outputs:
    EOS_ex: array of shape (size, 3) containing the values for EOS extension
    """

    size = int((ns[-1] - starts[0]) / n_step)

    # initialize array
    EOS_ex = np.zeros((size, 3))
    # set starting values at n = 2n0
    EOS_ex[0, 0] = starts[0]
    EOS_ex[0, 1] = starts[1]
    EOS_ex[0, 2] = starts[2]

    for k in range(size - 1):
        # n_i+1
        EOS_ex[k + 1, 0] = EOS_ex[k, 0] + n_step
        # p_i+1
        if cs2_func(EOS_ex[k, 0]) > 1:
            EOS_ex[k + 1, 1] = EOS_ex[k, 1] + n_step * ((EOS_ex[k, 1] + EOS_ex[k, 2]) / EOS_ex[k, 0])
        else:
            EOS_ex[k + 1, 1] = EOS_ex[k, 1] + n_step * (cs2_func(EOS_ex[k, 0])) * (
                        (EOS_ex[k, 1] + EOS_ex[k, 2]) / EOS_ex[k, 0])
        # e_i+1
        EOS_ex[k + 1, 2] = EOS_ex[k, 2] + n_step * ((EOS_ex[k, 1] + EOS_ex[k, 2]) / EOS_ex[k, 0])

    return EOS_ex


def extend_EOS_polytrop(n_step, ns, starts, Ks, gammas):
    """
    Function to extend the EOS by using linear segments in speed of sound

    Inputs:
    n_step: int defining step size in density
    ns: array containing the densities where segments begin/end
    starts: array containing the starting values for the EOS (density, pressure, energy)
    Ks: Array containing K parameters for polytropic segments
    gammas: Array containing gamma parameters for polytropic segments

    Outputs:
    EOS_ex: array of shape (size, 3) containing the values for EOS extension
    causality: boolean value =True if causality was never violated =False if it was at any point in extension
    """

    size = int((ns[-1] - ns[0]) / n_step)

    # initialize array
    EOS_ex = np.zeros((size, 3))
    # set starting values at n = 2n0
    EOS_ex[0, 0] = starts[0]
    EOS_ex[0, 1] = starts[1]
    EOS_ex[0, 2] = starts[2]

    i = 0

    for k in range(size - 1):
        # n_i+1
        EOS_ex[k + 1, 0] = EOS_ex[k, 0] + n_step
        # p_i+1
        if ns[i] < EOS_ex[k + 1, 0] < ns[i + 1]:
            EOS_ex[k + 1, 1] = EOS_ex[k, 1] + n_step * (Ks[i] * gammas[i] * (EOS_ex[k, 0] ** (gammas[i])))
        else:
            i += 1
            EOS_ex[k + 1, 1] = EOS_ex[k, 1] + n_step * (Ks[i] * gammas[i] * (EOS_ex[k, 0] ** (gammas[i])))
        # e_i+1
        EOS_ex[k + 1, 2] = EOS_ex[k, 2] + n_step * ((EOS_ex[k, 1] + EOS_ex[k, 2]) / EOS_ex[k, 0])

        # causality check
        if (EOS_ex[k + 1, 1] - EOS_ex[k, 1]) / (EOS_ex[k + 1, 2] - EOS_ex[k, 2]) > 1:
             break # stop extension when causality is broken

    EOS = EOS_ex[:k+1,:]
    return EOS


def stitch_EOS(small_EOS, EOS_ex):
    """
    Function to stitch EOS extension to original EOS

    Inputs:
    small_EOS: array containing original EOS
    EOS_ex: array containing EOS extension

    Outputs:
    tot_EOS: array containing total EOS
    """
    # get relevant sizes
    size_smol = small_EOS.shape[0] - 1  # -1 becuase we don't want last duplicated entry
    size_ex = EOS_ex.shape[0]

    # initialize array
    tot_EOS = np.zeros((size_smol + size_ex, small_EOS.shape[1]))

    tot_EOS[:size_smol, :] = small_EOS[:size_smol, :]
    tot_EOS[size_smol:, 0] = EOS_ex[:, 0]
    tot_EOS[size_smol:, 1] = EOS_ex[:, 1]
    tot_EOS[size_smol:, 2] = EOS_ex[:, 2]

    return tot_EOS


def extend(small_EOS, nsamp, rho_step=1e-3, rho_max=12, ext_type=None, max_gamma=6):
    """
    Function to simulate and extend an EOS three ways

    inputs
    nsim: int defining number of extensions to produce
    nsamp: int defining number of density points to sample when making extentions
    small_EOS: low density EOS to be extended
    rho_step: int defining the step size in density when extending
    rho_max: the maximum density of extension
    ext_type: the type of extension to do (speed of sound, linear P, and polytropic)
        'cs' (string) or 1 (int)
        'linear' or 2
        'polytrop' or 3
    a_gamma: only needed if doing polytropic extension, but characterized the distribution from
                which gamma's are drawn

    output
    EOS_tot = extended EOS stiched after small_EOS
    """
    n = small_EOS[:, 0]
    p = small_EOS[:, 1]
    e = small_EOS[:, 2]

    starts = [n[-1], p[-1], e[-1]]

    n0 = 0.16  # MeV/fm^3
    max_p = 350  # MeV/fm^3
    n_max = rho_max * n0

    ns = sample_ns(nsamp, n[-1], n_max)

    if ext_type == 'cs' or ext_type == 1:
        # derivative of pressure wrt energy
        dp_de = scipy.interpolate.CubicSpline(p, e).derivative(nu=1)

        # definition of speed of sound
        cs_start = np.sqrt(1 / dp_de(p[-1]))

        # sample cs
        cs = sample_cs(nsamp, cs_start)

        # make speed of sound squared function
        cs2_func = scipy.interpolate.interp1d(ns, cs ** 2)

        # extend EOS
        EOS_ex = extend_cs(rho_step, ns, starts, cs2_func)
        EOS_tot = stitch_EOS(small_EOS, EOS_ex)

        return EOS_tot, ns, cs

    elif ext_type == 'lin' or ext_type == 2:
        # derivative of pressure wrt energy
        dp_de = scipy.interpolate.CubicSpline(p, e).derivative(nu=1)
        # definition of speed of sound
        cs_start = np.sqrt(1 / dp_de(p[-1]))
        # sample cs
        cs = sample_cs(nsamp, cs_start)

        # extend EOS
        EOS_ex = extend_Plin(rho_step, ns, starts, cs)
        EOS_tot = stitch_EOS(small_EOS, EOS_ex)

        return EOS_tot, ns, cs

    elif ext_type == 'poly' or ext_type == 3:
        # sample parameters
        gammas, Ks = sample_polytrop(nsamp, ns, p[-1], max_gamma)

        EOS_ex = extend_EOS_polytrop(rho_step, ns, starts, Ks, gammas)

        EOS_tot = stitch_EOS(small_EOS, EOS_ex)

        return EOS_tot, ns, gammas, Ks

    else:
        print('The ext_type you entered does not match any of the allowed types. Please enter "cs" for a ' +
              'speed of sound extension, "lin" for an extension linear in pressure and "poly" ' +
              'for a polytropic extension')
    return