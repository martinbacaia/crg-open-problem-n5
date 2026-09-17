"""
Fase 2, paso 27: version de ALTA PRECISION (mpmath) del limite de Regge de n=5 + los 3 invariantes
base de grado 9 (parte 26) -- para explorar con precision suficiente la familia completa de 2
parametros que cancela el termino lider E^8, dado que en float64 solo un punto particular de esa
familia (v1=B1-2*B2) dio un resultado limpio (E^6) mientras otros puntos dieron ruido (leccion ya
aprendida en grado 11: float64 no alcanza para explorar cancelaciones finas).
"""

import itertools
import mpmath as mp

mp.mp.dps = 100

D = 6


def eta(a, b):
    return -a[0] * b[0] + sum(a[k] * b[k] for k in range(1, D))


def boost_from_rest(P, q_rest):
    M = mp.sqrt(-eta(P, P))
    P0 = P[0]
    Pvec = P[1:]
    beta_vec = [x / P0 for x in Pvec]
    beta2 = sum(b * b for b in beta_vec)
    gamma = P0 / M
    q0, qvec = q_rest[0], q_rest[1:]
    bq = sum(beta_vec[k] * qvec[k] for k in range(D - 1))
    q0p = gamma * (q0 + bq)
    if beta2 < mp.mpf('1e-40'):
        qvecp = qvec[:]
    else:
        qvecp = [qvec[k] + (gamma - 1) * (bq / beta2) * beta_vec[k] + gamma * q0 * beta_vec[k]
                  for k in range(D - 1)]
    return [q0p] + qvecp


def build_kinematics_n5(E, t0, m0sq, phi0):
    n1 = [mp.mpf(0)] * (D - 1)
    n1[0] = mp.mpf(1)
    E3 = E - m0sq / (4 * E)
    cos_theta = 1 + t0 / (2 * E * E3)
    sin_theta = mp.sqrt(1 - cos_theta ** 2)
    n3 = [mp.mpf(0)] * (D - 1)
    n3[0] = cos_theta
    n3[1] = sin_theta

    p1 = [-E] + [-E * x for x in n1]
    p2 = [-E] + [-E * (-x) for x in n1]
    p3 = [E3] + [E3 * x for x in n3]

    P45 = [-(p1[c] + p2[c] + p3[c]) for c in range(D)]

    e_orth = [mp.mpf(0)] * (D - 1)
    e_orth[2] = mp.mpf(1)
    n45 = [mp.cos(phi0) * n3[c] + mp.sin(phi0) * e_orth[c] for c in range(D - 1)]
    norm45 = mp.sqrt(sum(x * x for x in n45))
    n45 = [x / norm45 for x in n45]

    m0 = mp.sqrt(m0sq)
    p4_rest = [m0 / 2] + [(m0 / 2) * x for x in n45]
    p5_rest = [m0 / 2] + [(m0 / 2) * (-x) for x in n45]

    p4 = boost_from_rest(P45, p4_rest)
    p5 = boost_from_rest(P45, p5_rest)

    ps = {1: p1, 2: p2, 3: p3, 4: p4, 5: p5}
    return ps


def FdotF(ps, eps, i, j):
    pi, pj, ei, ej = ps[i], ps[j], eps[i], eps[j]
    return 2 * (eta(pi, pj) * eta(ei, ej) - eta(pi, ej) * eta(ei, pj))


def Cabc(ps, eps, k, a, b):
    pk, ek = ps[k], eps[k]
    return eta(pk, ps[a]) * eta(ek, ps[b]) - eta(pk, ps[b]) * eta(ek, ps[a])


def Msandwich(ps, eps, i, j, a, b):
    pi, pj, pa, pb = ps[i], ps[j], ps[a], ps[b]
    ei, ej = eps[i], eps[j]
    Vpj = eta(pa, pi) * eta(ei, pj) - eta(pa, ei) * eta(pi, pj)
    Vej = eta(pa, pi) * eta(ei, ej) - eta(pa, ei) * eta(pi, ej)
    return Vpj * eta(ej, pb) - Vej * eta(pj, pb)


def Basis1(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, i, j) * FdotF(ps, eps, l, m) * Cabc(ps, eps, k, i, j) * eta(ps[i], ps[l])


def Basis2(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, m) * Cabc(ps, eps, k, i, l)


def Basis3(ps, eps, i, j, l, m, k):
    return FdotF(ps, eps, l, m) * Msandwich(ps, eps, i, j, l, k) * Cabc(ps, eps, k, i, l)


def orbit_sum(ps, eps, fn):
    total = mp.mpf(0)
    for perm in itertools.permutations([1, 2, 3, 4, 5]):
        total += fn(ps, eps, *perm)
    return total


def make_eps(ps, seeds_pol, growing_partner):
    eps = {}
    for i in range(1, 6):
        e = seeds_pol[i][:]
        other = growing_partner[i]
        coef = eta(e, ps[i]) / eta(ps[other], ps[i])
        e = [e[c] - coef * ps[other][c] for c in range(D)]
        eps[i] = e
    return eps


if __name__ == "__main__":
    import random
    random.seed(11)
    t0, m0sq, phi0 = mp.mpf(-1), mp.mpf(3), mp.mpf('0.3')
    seeds_pol = {i: [mp.mpf(random.uniform(-3, 3)) for _ in range(D)] for i in range(1, 6)}
    growing_partner = {1: 2, 2: 1, 3: 2, 4: 1, 5: 1}

    Es = [mp.mpf(1000), mp.mpf(10000), mp.mpf(100000), mp.mpf(1000000)]
    data = []
    for E in Es:
        ps = build_kinematics_n5(E, t0, m0sq, phi0)
        eps = make_eps(ps, seeds_pol, growing_partner)
        b1 = orbit_sum(ps, eps, Basis1)
        b2 = orbit_sum(ps, eps, Basis2)
        b3 = orbit_sum(ps, eps, Basis3)
        data.append((E, b1, b2, b3))
        print(f"E={mp.nstr(E,6)}  B1={mp.nstr(b1,15)}  B2={mp.nstr(b2,15)}  B3={mp.nstr(b3,15)}")

    # coeficientes lider A_i = lim B_i/E^8, via richardson de alta precision (2 iteraciones)
    def richardson(vals, Es, power=8):
        A = [v / E ** power for v, E in zip(vals, Es)]
        while len(A) > 1:
            newA = []
            for i in range(len(A) - 1):
                E1, E2 = Es[i], Es[i + 1]
                newA.append((A[i + 1] * E2 - A[i] * E1) / (E2 - E1))
            A = newA
            Es = Es[1:]
        return A[0]

    A1 = richardson([r[1] for r in data], Es)
    A2 = richardson([r[2] for r in data], Es)
    A3 = richardson([r[3] for r in data], Es)
    print(f"\nA1={mp.nstr(A1,20)}\nA2={mp.nstr(A2,20)}\nA3={mp.nstr(A3,20)}")
    print(f"A1/A2 = {mp.nstr(A1/A2,20)}   A2/A3 = {mp.nstr(A2/A3,20)}")
