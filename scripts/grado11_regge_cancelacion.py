"""
Fase 2, paso 29: evaluar los 15 elementos base de grado 11 (parte 28) en el limite de Regge de n=5
(parte 23), y buscar la mejor combinacion lineal posible que minimice el crecimiento.
"""
import itertools
import numpy as np
from regge_limit_n5 import D, eta, build_kinematics_n5, FdotF, Cabc


def Msandwich(ps, eps, i, j, a, b):
    pi, pj, pa, pb = ps[i], ps[j], ps[a], ps[b]
    ei, ej = eps[i], eps[j]
    Vpj = eta(pa, pi) * eta(ei, pj) - eta(pa, ei) * eta(pi, pj)
    Vej = eta(pa, pi) * eta(ei, ej) - eta(pa, ei) * eta(pi, ej)
    return Vpj * eta(ej, pb) - Vej * eta(pj, pb)


def orbit_sum(ps, eps, fn):
    total = 0.0
    for perm in itertools.permutations([1, 2, 3, 4, 5]):
        total += fn(ps, eps, *perm)
    return total


def S1(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, i, j) * FdotF(ps, eps, l, m) * Cabc(ps, eps, k, i, j) * eta(ps[i], ps[j]) * eta(ps[i], ps[l])


def S2(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, i, j) * FdotF(ps, eps, l, m) * Cabc(ps, eps, k, i, j) * eta(ps[i], ps[l]) ** 2


def S3(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, i, j) * FdotF(ps, eps, l, m) * Cabc(ps, eps, k, i, j) * eta(ps[i], ps[l]) * eta(ps[i], ps[m])


def S4(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, m) * Cabc(ps, eps, k, i, j) * eta(ps[i], ps[l])


def S5(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, m) * Cabc(ps, eps, k, i, j) * eta(ps[k], ps[l])


def S6(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, m) * Cabc(ps, eps, k, i, l) * eta(ps[l], ps[m])


def S7(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, m) * Cabc(ps, eps, k, i, l) * eta(ps[i], ps[j])


def S8(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, m) * Cabc(ps, eps, k, i, l) * eta(ps[i], ps[l])


def S9(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, m) * Cabc(ps, eps, k, i, l) * eta(ps[k], ps[l])


def S10(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, k) * Cabc(ps, eps, k, i, j) * eta(ps[i], ps[l])


def S11(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, k) * Cabc(ps, eps, k, i, l) * eta(ps[l], ps[m])


def S12(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, k) * Cabc(ps, eps, k, i, l) * eta(ps[i], ps[l])


def S13(ps, eps, i, j, l, m, k):
    return Msandwich(ps, eps, i, j, l, m) * Msandwich(ps, eps, l, m, i, k) * Cabc(ps, eps, k, i, j)


def S14(ps, eps, i, j, l, m, k):
    return Msandwich(ps, eps, i, j, l, m) * Msandwich(ps, eps, l, m, i, k) * Cabc(ps, eps, k, l, m)


def S15(ps, eps, i, j, l, m, k):
    return Msandwich(ps, eps, i, j, l, m) * Msandwich(ps, eps, l, m, i, k) * Cabc(ps, eps, k, i, l)


BASIS = [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15]


if __name__ == "__main__":
    rng = np.random.RandomState(11)
    t0, m0sq = -1.0, 3.0
    seeds_pol = {i: rng.randn(D) for i in range(1, 6)}
    growing_partner = {1: 2, 2: 1, 3: 2, 4: 1, 5: 1}

    def make_eps(ps):
        eps = {}
        for i in range(1, 6):
            e = seeds_pol[i].copy()
            other = growing_partner[i]
            e = e - (eta(e, ps[i]) / eta(ps[other], ps[i])) * ps[other]
            eps[i] = e
        return eps

    Es = [300, 1000, 3000, 10000, 30000, 100000]
    data = np.zeros((len(Es), len(BASIS)))
    for row, E in enumerate(Es):
        ps = build_kinematics_n5(E, t0, m0sq, rng)
        eps = make_eps(ps)
        for col, fn in enumerate(BASIS):
            data[row, col] = orbit_sum(ps, eps, fn)

    Es_a = np.array(Es)
    print("Pendientes individuales de cada elemento base (S1..S15):")
    for col in range(len(BASIS)):
        vals = np.abs(data[:, col])
        slope = np.polyfit(np.log(Es_a[-3:]), np.log(vals[-3:]), 1)[0]
        print(f"  S{col+1}: slope in E = {slope:.4f}  (s12 = {slope/2:.4f})  value(E=1e5)={data[-1,col]:.4e}")
