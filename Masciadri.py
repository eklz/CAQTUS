from pylab import *
from scipy.ndimage import gaussian_filter1d
from .Maths import *
from .Dewan import *


def M_Mascidari(P,T,theta,z):
    """DOI: 10.1093/mnras/stw3111

    Args:
        P (array[float]): Pressure profile in hPa
        T (aray[float]): Temperature profile in K
        theta (aray[float]): Potential temperature profil in K
        z (array[float]): Altitude profile in m

    Returns:
        [int]
    """
    dTheta = np.gradient(theta,z)
    M=-((78e-6*P)/(theta*T))*dTheta

    return M


def Cn2_Mascidari(P, T ,z, S, gamma = 1.5):
    theta = potential_temp (T, P, P0 = 1000, pow = 0.286)
    alt_trop = trop_hght(T,z)

    L0_43 = L0_43_Dewan (z, alt_trop, S)

    M = M_Mascidari(P, T, theta, z)

    return gamma*(M**2)*L0_43
