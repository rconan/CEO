import math
import numpy as np
from scipy import special
import ceo

def variance(r0=None,L0=None,atmosphere=None):
    if atmosphere is not None:
        r0 = atmosphere.r0
        L0 = atmosphere.L0
    L0r0ratio= (L0/r0)**(5./3)
    return (24*math.gamma(6./5)/5.)**(5./6)* \
        (math.gamma(11./6)*math.gamma(5./6)/(2.*math.pi**(8./3)))*L0r0ratio

def covariance(rho,r0=None,L0=None,atmosphere=None):
    if atmosphere is not None:
        r0 = atmosphere.r0
        L0 = atmosphere.L0
    rho = np.array(rho)
    L0r0ratio= (L0/r0)**(5./3)
    cst      = (24.*math.gamma(6./5)/5.)**(5./6)* \
               (math.gamma(11./6)/(2.**(5./6)*math.pi**(8./3)))* \
                L0r0ratio
    out = np.zeros_like(rho)
    for k in range(rho.size):
        if rho[k]==0:
            out[k] = (24.*math.gamma(6./5)/5)**(5./6)* \
                     (math.gamma(11./6)*math.gamma(5./6)/(2.*math.pi**(8./3)))*L0r0ratio
        else:
            u  = 2.*math.pi*rho[k]/L0;
            out[k] = cst*u**(5./6)*special.kv(5./6,u)
    return out

def structure_function(rho,r0=None,L0=None,atmosphere=None):
    return 2*(variance(r0=r0,L0=L0,atmosphere=atmosphere) - covariance(rho,r0=r0,L0=L0,atmosphere=atmosphere))

    
