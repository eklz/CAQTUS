# %%
from .Legendre import *
# %%


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


# %%
J = roots_leg(10, [i/1000 for i in range(1000)])
J
# %%
