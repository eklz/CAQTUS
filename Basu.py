import numpy as np
from .Maths import potential_temp


def lt_thrope(theta, alt, fillna = False, size_filter = 11):
    
    sort = np.argsort(theta)
    deltaZ = []
    for i in range(len(theta)):
        delta = alt[i] - alt[sort[i]]
        deltaZ.append(delta)
    deltaZ = np.array(deltaZ)
    if not fillna:
        deltaZ[deltaZ == 0] = np.nan
    
    Lt = abs(moving_average(deltaZ, 11))
    return Lt

def cn2_basu (pres, temp, alt, c1 = 0.02, fillna = False, size_filter = 11):
    theta = potential_temp(temp, pres)
    grad =  np.gradient(np.sort(theta), alt)
    Lt = lt_thrope(theta, alt, fillna=fillna, size_filter=size_filter)
    Ct2 = c1*Lt**(5/3)*moving_average(grad,11)**2
    Cn2 = Cn2 = (7.9*1e-5*(pres/(temp**2)))**2*Ct2
    
    return Cn2