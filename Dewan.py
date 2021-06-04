from pylab import *
from scipy.ndimage import gaussian_filter1d
from .Maths import *

def L0_43_Dewan (z, alt_trop, S):
    """Calculate the outer scale of turbulence with Dewan's approach

    Args:
        z (aray[float]): Altitue profile in m of len n
        alt_trop (float): Altitude of the tropopause as given by trop_hght(T,z) in m
        S (aray[float]): Wind shear as given by  Moments.grad_wind(u, v, z)

    Returns:
        aray[float]: List of outer scale length for every altitude in z (len of len(z)) in m
    """

    L0_trop = 0.1**(4/3)*10**(1.64+42*S)
    L0_strat = 0.1**(4/3)*10**(0.5+50*S)
    L0_trop = L0_trop*(S<=0.02)
    L0_strat = L0_strat*(S<=0.045)

    res = (z<=alt_trop)*L0_trop + (z>alt_trop)*L0_strat
    return res


def M_Dewan(P, T, z):
    """DEWAN, Edmond M. A Model for C2n (optical turbulence) profiles using radiosonde data. 
    Directorate of Geophysics, Air Force Materiel Command, 1993.

     Args:
        P (array[float]): Pressure profile in hPa
        T (aray[float]): Temperature profile in K
        z (array[float]): Altitude profile in m

    Returns:
        [int]
    """
    dT = np.gradient(T,z)
    M=-((79e-6*P)/(T**2))*(dT+9.8e-3)
    
    return M


def Cn2_Dewan(P, T, z, S,  gamma = 2.8):
    """DEWAN, Edmond M. A Model for C2n (optical turbulence) profiles using radiosonde data. 
    Directorate of Geophysics, Air Force Materiel Command, 1993.

     Args:
        P (array[float]): Pressure profile in hPa
        T (aray[float]): Temperature profile in K
        z (array[float]): Altitude profile in m
        s (array[float]): Wind Shear
    Returns:
        Cn2 (array[float]): Cn2 m^(-2/3)
    """
    
    alt_trop = trop_hght(T,z)
    L0_43 = L0_43_Dewan (z, alt_trop, S)
    M = M_Dewan(P,T,z)
    
    return gamma*(M**2)*L0_43



def calc_cn2_date_Dewan(df, date, sigma = 0):
    
    
    data_day = df[df.date == date]
    u = data_day.u.values
    v = data_day.v.values
    z = data_day.alt.values
    T = data_day.temp.values
    P = data_day.press.values
    S = grad_wind(u, v, z)
    Cn2 = Cn2_Dewan(P, T, z, S, gamma = 2.8)
    if sigma != 0:
        Cn2 = gaussian_filter1d(Cn2, sigma=sigma)

    Cn2_nn = Cn2[Cn2 != 0]  # non 0 values
    z_nn = z[Cn2 != 0]  # non 0 values
    Ws_nn = data_day.wspeed.values[Cn2 != 0]  # non 0 values

    return Cn2_nn, z_nn, Ws_nn