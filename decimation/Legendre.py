# %%
import numpy as np
import numpy.polynomial.polynomial as nppol
import matplotlib.pyplot as plt
import cmath
# %%


def Legendre(n, t):
    Pk0 = np.poly1d([1])
    Pk1 = np.poly1d([1, -1/len(t)*np.sum(t)])

    for _ in range(2, n+1):
        Bk = np.sum([i*Pk1(i)**2 for i in t])
        Bk = Bk/np.sum([Pk1(i)**2 for i in t])

        Gk = np.sum([i*Pk1(i)*Pk0(i) for i in t])
        Gk = Gk/np.sum([Pk0(i)**2 for i in t])

        Pk2 = Pk1*np.poly1d([1, -Bk]) - Gk*Pk0

        Pk0, Pk1 = Pk1, Pk2

    return Pk2


# %%


def Laguerre(m, Pm, sigma=0.001):
    k, Pk = m, Pm
    res = []

    while k > 0:
        tau0 = complex(0)
        eps = 1
        Pk1 = np.polyder(Pk)
        Pk2 = np.polyder(Pk, 2)
        Hk = (k-1)*((k-1)*Pk1**2 - k*Pk*Pk2)
        print(k)
        while eps >= sigma:
            Pk2tau = Pk2(tau0)
            Pk1tau = Pk1(tau0)
            PkTau = Pk(tau0)

            Hktau = Hk(tau0)

            if np.abs((Pk1tau+cmath.sqrt(Hktau))) >= np.abs((Pk1tau-cmath.sqrt(Hktau))):
                tau1 = tau0 - (k*PkTau/(Pk1tau+cmath.sqrt(Hktau)))

            else:
                tau1 = tau0 - (k*PkTau/(Pk1tau-cmath.sqrt(Hktau)))

            eps = np.abs(tau1-tau0)
            tau0 = tau1
            print(k, eps)
        res.append(tau1)
        Pk = (Pk/np.poly1d([1, -tau1]))[0]
        k = k-1

    return np.sort(res)


def roots_leg(n, t):
    Pk0 = np.poly1d([1])  # 1
    Pk1 = np.poly1d([1, -1/len(t)*np.sum(t)])  # x - (1/n)*sum(t)
    Bk = np.sum([i*Pk1(i)**2 for i in t])
    Bk = Bk/np.sum([Pk1(i)**2 for i in t])
    J = np.zeros((n, n))
    for j in range(0, n-1):

        alpha = Bk

        Bk = np.sum([i*Pk1(i)**2 for i in t])
        Bk = Bk/np.sum([Pk1(i)**2 for i in t])
        Gk = np.sum([i*Pk1(i)*Pk0(i) for i in t])
        Gk = Gk/np.sum([Pk0(i)**2 for i in t])

        beta = np.sqrt(Gk)
        Pk2 = Pk1*np.poly1d([1, -Bk]) - Gk*Pk0

        Pk0, Pk1 = Pk1, Pk2

        J[j][j] = alpha
        J[j+1][j] = beta
        J[j][j+1] = beta
    J[n-1][n-1] = Bk
    return np.linalg.eig(J)[0]


def Lagrange(i, roots):
    L = np.poly1d([1])  # 1
    ai = roots[i]
    for a in np.delete(roots, i):
        L = L*1/(ai-a)*np.poly1d([1, -a])
    return L


def interp_gauss_legendre(n, t, y, norm=False):
    res = []
    min = np.min(t)
    max = np.max(t)
    if norm:
        t = (t-min)/(max-min)

    roots = roots_leg(n, t)
    roots = np.sort(roots)
    if len(t) != len(y):
        raise ValueError(
            f't and y must have the same length, t has length {len(t)} and y {len(y)}')

    for i in range(n):
        L = Lagrange(i, roots)
        norm = np.sum([L(j)**2 for j in t])
        sum = np.sum([L(t[k])*y[k] for k in range(len(y))])
        res.append(1/norm * sum)

    if norm:
        roots = roots*(max-min) + min

    return roots, res


def interp_gauss_legendre_pp(n, lt, ly):
    roots, res = [], []
    if len(lt) != len(ly):
        raise AttributeError('lt and ly must have the same length')
    for i in range(len(lt)):
        t = lt[i]
        y = ly[i]
        X, Y = interp_gauss_legendre(n, t, y, norm=True)
        roots = np.concatenate((roots, X))
        res = np.concatenate((res, Y))

    return roots, res


# %%
if __name__ == "__main__":
    Pm = Legendre(10, [i/1000 for i in range(1000)])
    Laguerre(10, Pm)
# %%
