from scipy import special
import numpy as np


def find_coeffs(r0, theta0, h0, hmax, lambdam=1.55e-6, theta=0):
    H1 = 3*h0**(4/3)*((1/hmax)**(1/3)-(1/h0)**(1/3))
    C1 = 4.05*10**(-13)*(1-np.exp(-hmax/1500))
    k = 2*np.pi/lambdam
    R0 = r0**(-5/3)*np.cos(theta)/(0.423*k**2)
    H2 = 3/4*h0**(4/3)*(hmax**(4/3)-h0**(4/3))
    C2 = 2.7*10**(-16)*1500**(8/3)*special.gamma(8 /
                                                 3)*(1-special.gammaincc(8/3, hmax/1500)) 
    # gamaincc : Regularized upper incomplete gamma function.
    THETA0 = theta0**(-5/3)*np.cos(theta)**(8/3)/(2.914*k**2)
    GAMMA = 10**5*special.gamma(38/3)*(1-special.gammaincc(38/3, hmax/1000))/(
        special.gamma(11)*(1-special.gammaincc(11, hmax/1000)))

    CN20 = (THETA0 - C2 - GAMMA*R0 + GAMMA*C1)/(GAMMA*H1 + H2)

    V = 10**17/(special.gamma(11) *
                (1-special.gammaincc(11, hmax/1000)))*(R0-C1+H1*CN20)

    v = np.sqrt(V/0.00594)*27

    return v, CN20

def gen_HAP(v, cn20, h0, linspace):
    h = np.array(linspace)
    cn2 = 0.00594*(v/27)**2*(1e-5*h)**(10)*np.exp (-h/1000)+1e-16*2.7*np.exp(-h/1500)+cn20*(h0/h)**(4/3)
    
    return cn2