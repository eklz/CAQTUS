import numpy as np
from scipy.ndimage import gaussian_filter1d

def calc_moments(Cn2, alt, windSpeed, lambda_m = 1.55e-6, zenithAngle = 0, 
                 fried = True, seeing = True, isoplanatic = True, coherenceTime = True, rmzeros = True):
    
    def int(x,t):
        x = np.array(x)
        t = np.array(t)
        dt = t[1:]-t[:-1]
        return np.sum(x[:-1]*dt)

    if rmzeros : 
        nzeros = Cn2.nonzero()[0]
        Cn2 = Cn2[nzeros]
        alt = alt[nzeros]
        windSpeed = windSpeed[nzeros]

    r0 = np.power((0.423*np.cos(zenithAngle)**(-1)*
                   (2*np.pi/lambda_m)**2*int(Cn2, alt)), (-3/5))
    
    res = {}

    if fried:
        res['r0'] = r0

    if seeing:
        res['seeing'] = (0.98*lambda_m/r0)

    if isoplanatic:
        res['theta0'] = np.power((2.914*(2*np.pi/lambda_m)**2 *
                                  np.cos(zenithAngle)**(-8/3) *
                                  int(Cn2*np.power(alt, (5/3)), alt)), (-3/5))*1e6

    if coherenceTime:
        res['tau0'] = np.power((2.914*(2*np.pi/lambda_m)**2 *
                                np.cos(zenithAngle)**(-8/3) *
                                int(Cn2*np.power(windSpeed, (5/3)), alt)), (-3/5))

    return res


def grad_wind(u, v, z):
    grad_u = np.gradient(u, z)
    grad_v = np.gradient(v, z)

    return np.sqrt(grad_u**2 + grad_v**2)




def potential_temp (T, P, P0 = 1000, pow = 0.286):
    """Calculate the potential temperature from temperature and pressure

    Args:
        T (aray[float]): Temperature profile in K
        P (array[float]): Pressure profile in hPa
        P0 (int, optional): Standard reference pressure. Defaults to 1000 hPa
        pow (float, optional): R/c_{p}=0.286 for air (meteorology). R is the gas constant of air, 
        and c_{p} is the specific heat capacity at a constant pressure.

    Returns:
        array[float]: Potential temperature in K 
    """
    return T*((P0/P)**pow)


def trop_hght(T,z):
    """Calculate the height of the tropopause
    
    Args:
        T (aray[float]): Temperature profile in K of len n 
        z (aray[float]): Altitue profile in m of len n

    Returns:
        float: Tropopause height in m
    """
    dT = -np.gradient(T,z)
    for i in range(len(dT)):
        if dT[i]<2.0:
            z_2km=(z>z[i]) & (z<=(z[i]+2000))
            if (np.mean(dT[z_2km])<2) and (np.mean(dT[z_2km])<0)  and (z[i]>8000):
                ztrop = z[i]
                return ztrop
    return 0




def moving_average(list, size):
    return np.convolve(list, np.ones(size), 'same')/size

def rm_zeros(x):
    return x[x.nonzeros()[0]]

