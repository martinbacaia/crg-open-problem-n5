"""
Fase 2, paso 18b: version de ALTA PRECISION (mpmath, 60 digitos) del calculo de rango en grado 11,
para diagnosticar si la inestabilidad observada con float64 (rank oscilando entre 3, 5, 9, 10, 12
segun la muestra -- ver 21-teorema-paridad-grados-8-10-11.md) es un problema de precision numerica
o un fenomeno real del muestreo.
"""

import itertools
import mpmath as mp

mp.mp.dps = 60  # 60 digitos decimales

D = 6


def gen_kinematics(seed):
    # generador congruencial lineal simple para pseudo-aleatoriedad reproducible en alta precision
    state = [seed * 2654435761 % (2**32)]

    def rnd():
        state[0] = (state[0] * 1103515245 + 12345) % (2**31)
        return mp.mpf(state[0]) / mp.mpf(2**31)

    def randn():
        # Box-Muller con alta precision
        u1, u2 = rnd(), rnd()
        if u1 < mp.mpf('1e-30'):
            u1 = mp.mpf('1e-30')
        return mp.sqrt(-2 * mp.log(u1)) * mp.cos(2 * mp.pi * u2)

    def rand_null():
        nvec = [randn() for _ in range(D - 1)]
        norm = mp.sqrt(sum(x * x for x in nvec))
        nvec = [x / norm for x in nvec]
        E = mp.mpf(1) + rnd() * 2
        return E, nvec

    E2, n2 = rand_null()
    E3, n3 = rand_null()
    E4, n4 = rand_null()
    n5 = [randn() for _ in range(D - 1)]
    norm5 = mp.sqrt(sum(x * x for x in n5))
    n5 = [x / norm5 for x in n5]

    def p1_mass2(E5):
        E0 = E2 + E3 + E4 + E5
        spatial = [E2 * n2[i] + E3 * n3[i] + E4 * n4[i] + E5 * n5[i] for i in range(D - 1)]
        return -E0 ** 2 + sum(x * x for x in spatial)

    E5 = mp.findroot(p1_mass2, mp.mpf(5))

    def vec(E, n):
        return [E] + [E * x for x in n]

    p2v, p3v, p4v, p5v = vec(E2, n2), vec(E3, n3), vec(E4, n4), vec(E5, n5)
    p1v = [-(p2v[k] + p3v[k] + p4v[k] + p5v[k]) for k in range(D)]

    def eta(a, b):
        return -a[0] * b[0] + sum(a[k] * b[k] for k in range(1, D))

    ps = {1: p1v, 2: p2v, 3: p3v, 4: p4v, 5: p5v}
    for k in range(1, 6):
        assert abs(eta(ps[k], ps[k])) < mp.mpf('1e-40'), f"p{k}^2 no nulo: {eta(ps[k],ps[k])}"
    total = [sum(ps[k][c] for k in range(1, 6)) for c in range(D)]
    assert max(abs(x) for x in total) < mp.mpf('1e-40')

    eps = {}
    for i in range(1, 6):
        e = [randn() for _ in range(D)]
        other = 2 if i != 2 else 3
        coef = eta(e, ps[i]) / eta(ps[other], ps[i])
        e = [e[c] - coef * ps[other][c] for c in range(D)]
        eps[i] = e

    return ps, eps, eta


def FdotF(eta, ps, eps, i, j):
    pi, pj, epsi, epsj = ps[i], ps[j], eps[i], eps[j]
    return 2 * (eta(pi, pj) * eta(epsi, epsj) - eta(pi, epsj) * eta(epsi, pj))


def Cabc(eta, ps, eps, k, a, b):
    pk, epsk = ps[k], eps[k]
    return eta(pk, ps[a]) * eta(epsk, ps[b]) - eta(pk, ps[b]) * eta(epsk, ps[a])


def orbit_sum(eta, ps, eps, base_fn):
    total = mp.mpf(0)
    for perm in itertools.permutations([1, 2, 3, 4, 5]):
        total += base_fn(eta, ps, eps, *perm)
    return total


SLOT_NAMES = ["i", "j", "l", "m", "k"]
ALL_PAIRS = list(itertools.combinations(SLOT_NAMES, 2))


def make_seed_A(c_contract, mult_two_pairs):
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        slot_val = {"i": i, "j": j, "l": l, "m": m, "k": k}
        if c_contract == "P":
            c = Cabc(eta, ps, eps, k, i, j)
        elif c_contract == "Q":
            c = Cabc(eta, ps, eps, k, l, m)
        else:
            c = Cabc(eta, ps, eps, k, i, l)
        mult = mp.mpf(1)
        for (x, y) in mult_two_pairs:
            a, b = slot_val[x], slot_val[y]
            mult *= eta(ps[a], ps[b])
        return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * c * mult

    return seed_fn


def make_seed_B_v1():
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        return (
            FdotF(eta, ps, eps, i, j)
            * Cabc(eta, ps, eps, l, i, j)
            * Cabc(eta, ps, eps, m, i, j)
            * Cabc(eta, ps, eps, k, i, j)
        )
    return seed_fn


def make_seed_B_v2():
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        return (
            FdotF(eta, ps, eps, i, j)
            * Cabc(eta, ps, eps, l, i, j)
            * Cabc(eta, ps, eps, m, i, k)
            * Cabc(eta, ps, eps, k, j, l)
        )
    return seed_fn


def mp_matrix_rank(rows, tol=mp.mpf('1e-30'), verbose=False):
    """rango numerico via eliminacion gaussiana con pivoteo parcial, en alta precision."""
    A = [row[:] for row in rows]
    nrows = len(A)
    ncols = len(A[0]) if nrows else 0
    rank = 0
    pivots = []
    for col in range(ncols):
        pivot_row = None
        pivot_val = tol
        for r in range(rank, nrows):
            if abs(A[r][col]) > pivot_val:
                pivot_val = abs(A[r][col])
                pivot_row = r
        if pivot_row is None:
            continue
        A[rank], A[pivot_row] = A[pivot_row], A[rank]
        pv = A[rank][col]
        pivots.append(pv)
        for r in range(nrows):
            if r != rank and abs(A[r][col]) > 0:
                factor = A[r][col] / pv
                for c in range(col, ncols):
                    A[r][c] -= factor * A[rank][c]
        rank += 1
        if rank == nrows:
            break
    if verbose:
        print("  pivotes encontrados:", [mp.nstr(p, 5) for p in pivots])
        # despues del ultimo pivote real, mostrar el mayor residual restante (deberia ser ~0)
        max_residual = mp.mpf(0)
        for col in range(ncols):
            for r in range(rank, nrows):
                max_residual = max(max_residual, abs(A[r][col]))
        print("  mayor residual tras el ultimo pivote:", mp.nstr(max_residual, 5))
    return rank


if __name__ == "__main__":
    seeds = {}
    mono4 = list(itertools.combinations_with_replacement(ALL_PAIRS, 2))
    for c in ["P", "Q", "M"]:
        for mp_ in mono4:
            seeds[f"A_{c}__d{mp_[0]}*d{mp_[1]}"] = make_seed_A(c, mp_)
    seeds["B_v1"] = make_seed_B_v1()
    seeds["B_v2"] = make_seed_B_v2()
    seed_names = list(seeds.keys())
    print(f"Total semillas: {len(seed_names)}")

    import sys as _sys
    n_samples = int(_sys.argv[1]) if len(_sys.argv) > 1 else 15
    rows = []
    for kseed in range(1, n_samples + 1):
        ps, eps, eta = gen_kinematics(kseed)
        row = [orbit_sum(eta, ps, eps, seeds[name]) for name in seed_names]
        rows.append(row)
        print(f"  muestra {kseed}/{n_samples} calculada")

    # descartar columnas identicamente nulas
    norms = []
    for c in range(len(seed_names)):
        norms.append(mp.sqrt(sum(rows[r][c] ** 2 for r in range(n_samples))))
    nz_idx = [c for c, nrm in enumerate(norms) if nrm > mp.mpf('1e-25')]
    print(f"\nSemillas no nulas: {len(nz_idx)} de {len(seed_names)}")

    Mn = [[rows[r][c] / norms[c] for c in nz_idx] for r in range(n_samples)]
    rank = mp_matrix_rank(Mn, verbose=True)
    print(f"\nRango (alta precision, {mp.mp.dps} digitos, {n_samples} muestras) = {rank}")
