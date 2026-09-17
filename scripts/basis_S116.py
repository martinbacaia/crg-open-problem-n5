"""
BASIS S1..S16 (S1..S15 = familias A+B+C+D ya conocidas, dimension 15 confirmada en la parte 35;
S16 = la direccion nueva de la topologia "sandwich triple" T, parte 36, que sube el rango a 16).

S16(i,j,l,m,k) := T(i,j,k; l,m) * C(l;i,j) * C(m;i,k)
  (T con field-strengths i,j,k, patas de sandwich a,b=l,m -- la variante mas limpia encontrada en
  grado11_familia_T_triple_sandwich.py). Tsandwich_cached verificada exacta contra la cadena
  generica de referencia en tsandwich_cached_check.py.

Backend rapido (lseries.py, LSeries), reusable para las 3 familias de Regge (solo cambia como se
construyen ps_ls/eps_ls).
"""
from lseries import LSeries, eta_vec

P = list(range(1, 6))


def compute_S116(ps_ls, eps_ls, D):
    def eta(a, b):
        return eta_vec(a, b, D)

    eta_pp, eta_pe, eta_ee = {}, {}, {}
    for a in P:
        for b in P:
            eta_pp[(a, b)] = eta(ps_ls[a], ps_ls[b])
            eta_pe[(a, b)] = eta(ps_ls[a], eps_ls[b])
            eta_ee[(a, b)] = eta(eps_ls[a], eps_ls[b])

    def FdotF(i, j):
        return (eta_pp[(i, j)] * eta_ee[(i, j)] - eta_pe[(j, i)] * eta_pe[(i, j)]) * 2

    def Cabc(k, a, b):
        return eta_pp[(k, a)] * eta_pe[(b, k)] - eta_pp[(k, b)] * eta_pe[(a, k)]

    def Msandwich(i, j, a, b):
        pi_pa = eta_pp[(a, i)]; ei_pj = eta_pe[(j, i)]; pa_ei = eta_pe[(a, i)]; pi_pj = eta_pp[(i, j)]
        Vpj = pi_pa * ei_pj - pa_ei * pi_pj
        ei_ej = eta_ee[(i, j)]; pi_ej = eta_pe[(i, j)]
        Vej = pi_pa * ei_ej - pa_ei * pi_ej
        ej_pb = eta_pe[(b, j)]; pj_pb = eta_pp[(j, b)]
        return Vpj * ej_pb - Vej * pj_pb

    def Tsandwich(i, j, k, a, b):
        vp = {x: eta_pp[(a, x)] for x in P}
        ve = {x: eta_pe[(a, x)] for x in P}
        for m in (i, j, k):
            new_vp = {}
            new_ve = {}
            for x in P:
                new_vp[x] = vp[m] * eta_pe[(x, m)] - ve[m] * eta_pp[(m, x)]
                new_ve[x] = vp[m] * eta_ee[(m, x)] - ve[m] * eta_pe[(m, x)]
            vp, ve = new_vp, new_ve
        return vp[b]

    FF, Cc, Msand = {}, {}, {}
    for i in P:
        for j in P:
            if i < j:
                FF[(i, j)] = FdotF(i, j)
    for k in P:
        for a in P:
            for b in P:
                if len({k, a, b}) == 3:
                    Cc[(k, a, b)] = Cabc(k, a, b)
    for i in P:
        for j in P:
            if i == j:
                continue
            for a in P:
                for b in P:
                    if len({i, j, a, b}) == 4:
                        Msand[(i, j, a, b)] = Msandwich(i, j, a, b)

    def getFF(i, j):
        return FF[(i, j)] if i < j else FF[(j, i)]

    def E(i, j):
        return eta_pp[(i, j)]

    def S1(i, j, l, m, k):  return getFF(i, j) * getFF(l, m) * Cc[(k, i, j)] * E(i, j) * E(i, l)
    def S2(i, j, l, m, k):  return getFF(i, j) * getFF(l, m) * Cc[(k, i, j)] * E(i, l) ** 2
    def S3(i, j, l, m, k):  return getFF(i, j) * getFF(l, m) * Cc[(k, i, j)] * E(i, l) * E(i, m)
    def S4(i, j, l, m, k):  return getFF(l, m) * Msand[(i, j, l, m)] * Cc[(k, i, j)] * E(i, l)
    def S5(i, j, l, m, k):  return getFF(l, m) * Msand[(i, j, l, m)] * Cc[(k, i, j)] * E(k, l)
    def S6(i, j, l, m, k):  return getFF(l, m) * Msand[(i, j, l, m)] * Cc[(k, i, l)] * E(l, m)
    def S7(i, j, l, m, k):  return getFF(l, m) * Msand[(i, j, l, m)] * Cc[(k, i, l)] * E(i, j)
    def S8(i, j, l, m, k):  return getFF(l, m) * Msand[(i, j, l, m)] * Cc[(k, i, l)] * E(i, l)
    def S9(i, j, l, m, k):  return getFF(l, m) * Msand[(i, j, l, m)] * Cc[(k, i, l)] * E(k, l)
    def S10(i, j, l, m, k): return getFF(l, m) * Msand[(i, j, l, k)] * Cc[(k, i, j)] * E(i, l)
    def S11(i, j, l, m, k): return getFF(l, m) * Msand[(i, j, l, k)] * Cc[(k, i, l)] * E(l, m)
    def S12(i, j, l, m, k): return getFF(l, m) * Msand[(i, j, l, k)] * Cc[(k, i, l)] * E(i, l)
    def S13(i, j, l, m, k): return Msand[(i, j, l, m)] * Msand[(l, m, i, k)] * Cc[(k, i, j)]
    def S14(i, j, l, m, k): return Msand[(i, j, l, m)] * Msand[(l, m, i, k)] * Cc[(k, l, m)]
    def S15(i, j, l, m, k): return Msand[(i, j, l, m)] * Msand[(l, m, i, k)] * Cc[(k, i, l)]
    def S16(i, j, l, m, k): return Tsandwich(i, j, k, l, m) * Cc[(l, i, j)] * Cc[(m, i, k)]

    BASIS = [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16]

    import itertools
    totals = [LSeries.zero() for _ in BASIS]
    for perm in itertools.permutations(P):
        i, j, l, m, k = perm
        for idx, fn in enumerate(BASIS):
            totals[idx] = totals[idx] + fn(i, j, l, m, k)
    return totals
