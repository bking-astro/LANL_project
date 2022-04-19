import numpy as np
from scipy import interpolate
from scipy.integrate import cumtrapz,odeint


CONV_MeV_fm3_to_g_cm3 = 1.78266181e-36 * 1e48
CONV_MeV_fm3_to_dyn_cm2 = 1.78266181e-36*2.99792458e8**2*1e52


G = 6.6743 * 10**(-8)      # Newton's gravitational constant in cgs units
c = 2.99792458 * 10**10    # speed of light in cgs units
M_sun = 1.476 * 10**(5)    # Mass of the sun in geometrized units




def tovrhs(y,t,eos1,eos2,de_dp):
    
    r = y[0]
    m = y[1]
    x = y[2]

    eps = eos1(t)
    p = eos2(t)
    
    dr_dh = -r*(r-2*m)/(m + 4* np.pi * r**3 * p )
    
    dm_dh = 4 * np.pi * r**2 * eps * dr_dh
    
    f_inv = 1/(1 - 2*m/r)
    F = f_inv*(1 - 4*np.pi*(eps-p)*r**2)
    r2Q = f_inv*( 4*np.pi*r**2*(5*eps+9*p+(eps+p)*de_dp(t)) - 6 - f_inv*((4*m**2)/(r**2))*(1 + 4*np.pi*r**3*p/m)**2 )
    dx_dh = dr_dh*(-x**2 - x*F - r2Q)*(1/r)
    
    return [dr_dh, dm_dh, dx_dh]
    

def tov_solve(epsilon,press, size=100):

    epsilon = epsilon  *   CONV_MeV_fm3_to_g_cm3  * G * M_sun**2 / c**2
    press = press     *   CONV_MeV_fm3_to_dyn_cm2     * G * M_sun**2 / c**4
    enthalpy = cumtrapz(1/(epsilon+press), press, initial=0)
    
    enthalpy = np.sort(enthalpy)
    

    spl1 = interpolate.splrep(enthalpy, epsilon,k=3) 
    spl2 = interpolate.splrep(enthalpy, press,k=3) 
    spl3 = interpolate.splrep(press, epsilon,k=3)
    

    def eos1(h):
        # if h < np.amin(enthalpy):
        #     h = np.amin(enthalpy)
        return interpolate.splev(h, spl1, der=0)
    
#     def eos1der(h):
#         # if h < np.amin(enthalpy):
#         #     h = np.amin(enthalpy)
#         return interpolate.splev(h, spl1, der=1)
    
    
    def eos2(h):
        # if h < np.amin(enthalpy):
        #     h = np.amin(enthalpy)
        return interpolate.splev(h, spl2, der=0)
    
#     def eos2der(h):
#         # if h < np.amin(enthalpy):
#         #     h = np.amin(enthalpy)
#         return interpolate.splev(h, spl2, der=1)

    def de_dp(h):
        # if h < np.amin(enthalpy):
        #     h = np.amin(enthalpy)
        p = interpolate.splev(h, spl2, der=0)
        return interpolate.splev(p, spl3, der=1)
    
    start = np.log(0.05)
    stop = np.log(1.8)
    central_enthalpies = np.exp( np.linspace(start,stop,size) )
    
    Radius = np.zeros_like(central_enthalpies)
    Mass = np.zeros_like(central_enthalpies)
    Lambda = np.zeros_like(central_enthalpies)
    
    
    for i,h_c in enumerate(central_enthalpies):
    
        r0 = 1e-5
        m0 = 4/3 * np.pi * r0**3 * eos1(h_c)
        x0 = 2
        
        initial = r0, m0, x0
        
        h = np.linspace(h_c,0,1000)
            
        ans =  odeint(tovrhs,initial,h, args = (eos1,eos2,de_dp),
                      rtol = 1.49012e-7, atol = 1.49012e-7)  
        
        
        Radius[i] = ans[-1,0] * M_sun * 10**(-5)
        Mass[i] =  ans[-1,1]
        C = ans[-1,1]/ans[-1,0]
        xR = ans[-1,2]
        
        #k2
        num = ((1 - 2 * C) ** 2) * (2 - xR + 2 * C * (xR - 1))

        a1 = 6 - 3 * xR + 3 * C * (5 * xR - 8)
        a2 = 13 - 11 * xR + C * (3 * xR - 2) + (2 * C ** 2) * (1 + xR)
        a3 = 2 - xR + 2 * C * (xR - 1)

        den = 2 * C * a1 + 4 * (C ** 3) * a2 + 3 * ((1 - 2 * C) ** 2) * a3 * np.log(1 - 2 * C)
        
        Lambda[i] = (16/ 15) * (num/den)
        
        if Mass[i] < Mass[i-1] and Mass[i] > 1:
            Mass = Mass[:i]
            Radius = Radius[:i]
            Lambda = Lambda[:i]
            break
        
    
    
    return Radius,Mass,Lambda

