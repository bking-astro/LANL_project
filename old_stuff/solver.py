import numpy as np
from scipy import interpolate
from scipy.integrate import cumtrapz,odeint


CONV_MeV_fm3_to_g_cm3 = 1.78266181e-36 * 1e48
CONV_MeV_fm3_to_dyn_cm2 = 1.78266181e-36*2.99792458e8**2*1e52


G = 6.6743 * 10**(-8)      # Newton's gravitational constant in cgs units
c = 2.99792458 * 10**10    # speed of light in cgs units
M_sun = 1.476 * 10**(5)    # Mass of the sun in geometrized units




def tovrhs(y,t,eos1,eos2):
    
    r = y[0]
    m = y[1]

    eps = eos1(t)
    p = eos2(t)
    
    dr_dh = -r*(r-2*m)/(m + 4* np.pi * r**3 * p )
    
    dm_dh = 4 * np.pi * r**2 * eps * dr_dh
    
    return [dr_dh,dm_dh]
    



def tov_solve(epsilon,press, size=100):

    epsilon = epsilon  *   CONV_MeV_fm3_to_g_cm3  * G * M_sun**2 / c**2
    press = press     *   CONV_MeV_fm3_to_dyn_cm2     * G * M_sun**2 / c**4
    enthalpy = cumtrapz(1/(epsilon+press), press, initial=0)
    
    enthalpy = np.sort(enthalpy)
    

    spl1 = interpolate.splrep(enthalpy, epsilon,k=3) 
    spl2 = interpolate.splrep(enthalpy, press,k=3) 
    

    def eos1(h):
        # if h < np.amin(enthalpy):
        #     h = np.amin(enthalpy)
        return interpolate.splev(h, spl1, der=0)
    
    
    def eos2(h):
        # if h < np.amin(enthalpy):
        #     h = np.amin(enthalpy)
        return interpolate.splev(h, spl2, der=0)
    

    central_enthalpies = np.linspace(0.05,1.4,size)
    
    Radius = np.zeros_like(central_enthalpies)
    Mass = np.zeros_like(central_enthalpies)
    
    
    for i,h_c in enumerate(central_enthalpies):
    
        r0 = 0.0001
        m0 = 4/3 * np.pi * r0**3 * eos1(h_c)
        
        initial = r0, m0
        
        h = np.linspace(h_c,0,1000)
            
        ans =  odeint(tovrhs,initial,h, args = (eos1,eos2),
                      rtol = 1.49012e-7, atol = 1.49012e-6)  
        
        
        Radius[i] = ans[-1,0] * M_sun * 10**(-5)
        Mass[i] =  ans[-1,1]
        
        if Mass[i] < Mass[i-1] and Mass[i] > 0.5:
            Mass = Mass[:i]
            Radius = Radius[:i]
            break
        
    
    
    return Radius,Mass

