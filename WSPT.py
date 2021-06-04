import numpy as np
from .Maths import *
from .Masciadri import * 

def Cn2_WSPT(P, T, z, u, v, c = 0.3):
    theta = potential_temp (T, P, P0 = 1000, pow = 0.286)
    S = grad_wind(u, v, z)
    theta_sorted = np.sort(theta)
    delta_theta = theta - theta_sorted
    
    Lw = np.sqrt( (delta_theta)/(np.gradient(theta_sorted, z)) * np.sqrt( (u*v)/(S**2) ) )
    
    ct2 = c*(Lw**(4/3))*(np.gradient(theta_sorted, z)**2)
    
    cn2 = (7.9*1e-5*(P/T**2))**2*ct2
    
    return cn2


def Cn2_HMNSP99(P, T, z, u, v):
    dT = np.gradient(T, z)
    theta = potential_temp(T, P)
    S = grad_wind(u, v, z)
    alt_trop = trop_hght(T,z)
    
    L0_trop = 0.1**(4/3)*10**(0.362+16.728*S-192.347*dT)
    L0_strat = 0.1**(4/3)*10**(0.752+13.819*S-57.784*dT)
    L0_43 = (z<=alt_trop)*L0_trop + (z>alt_trop)*L0_strat
    
    M = M_Mascidari(P,T,theta,z)
    
    return 1.5*(M**2)*L0_43