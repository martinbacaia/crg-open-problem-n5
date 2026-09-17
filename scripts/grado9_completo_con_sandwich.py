"""
Fase 2, paso 26: recalculo RIGUROSO de la dimension real del espacio de invariantes de grado 9,
incluyendo sistematicamente la topologia "sandwich" (M^{ij}_{ab} := p_a.F^i.F^j.p_b) encontrada en
la parte 25 -- que la parte 20 no habia considerado (de ahi la correccion: "dimension 1" paso a
"dimension >= 2").

Familias de semillas de grado 9 (grado = 2 + 4 + 3, cubriendo las 5 patas exactamente una vez):

  (A) "2 parejas F:F + 1 singlete C" x corrector Mandelstam grado 2 (ya usada en la parte 20,
      dimension 1 encontrada AHI, pero puede no ser toda la historia combinada con B).
  (B) "1 pareja F:F + 1 sandwich + 1 singlete C" (grado 2+4+3=9, SIN corrector, estructura
      genuinamente distinta) -- se enumeran sistematicamente las variantes de contraccion:
      el sandwich M^{ij}_{ab} necesita 2 momentos externos (a,b) elegidos de {l,m,k} (6 ordenados,
      ya que M no es simetrico en (a,b) -- verificado en la parte 25); el singlete C^k_{ab} necesita
      2 momentos de {i,j,l,m} (6 no ordenados, C es antisimetrico en (a,b), asi que el signo ya lo
      cubre el orden).

Metodo: exactamente el mismo que ya se uso con exito en partes 20-21 -- simetrizar cada semilla
sobre la orbita completa de S5, evaluar en VARIAS cinematicas numericas independientes, armar la
matriz (muestras x semillas), y calcular su RANGO NUMERICO (no asumido). Dado que grado 9 en
float64 ya demostro ser confiable en la parte 20 (salto de 14 ordenes de magnitud con pocas
muestras) -- a diferencia de grado 11 -- se usa float64 aqui, pero se vuelve a chequear el salto de
valores singulares con cuidado antes de confiar en el resultado (y se re-verifica con mas muestras
si el salto no es limpio).
"""

import itertools
import numpy as np
from scipy.optimize import brentq

D = 6


def gen_kinematics(seed):
    rng = np.random.RandomState(seed)

    def rand_null():
        nvec = rng.randn(D - 1)
        nvec /= np.linalg.norm(nvec)
        E = rng.uniform(1, 3)
        return E, nvec

    E2, n2 = rand_null()
    E3, n3 = rand_null()
    E4, n4 = rand_null()
    n5 = rng.randn(D - 1)
    n5 /= np.linalg.norm(n5)

    def p1_mass2(E5):
        E0 = E2 + E3 + E4 + E5
        spatial = E2 * n2 + E3 * n3 + E4 * n4 + E5 * n5
        return -E0 ** 2 + spatial.dot(spatial)

    xs = np.linspace(-50, 50, 400)
    vals = [p1_mass2(x) for x in xs]
    E5 = None
    for i in range(len(xs) - 1):
        if vals[i] * vals[i + 1] < 0:
            E5 = brentq(p1_mass2, xs[i], xs[i + 1])
            break
    assert E5 is not None

    def vec(E, n):
        return np.concatenate([[E], E * n])

    p2v, p3v, p4v, p5v = vec(E2, n2), vec(E3, n3), vec(E4, n4), vec(E5, n5)
    p1v = -(p2v + p3v + p4v + p5v)

    def eta(a, b):
        return -a[0] * b[0] + a[1:].dot(b[1:])

    ps = {1: p1v, 2: p2v, 3: p3v, 4: p4v, 5: p5v}
    for k in range(1, 6):
        assert abs(eta(ps[k], ps[k])) < 1e-8
    assert np.linalg.norm(sum(ps.values())) < 1e-8

    eps = {}
    for i in range(1, 6):
        e = rng.randn(D)
        other = 2 if i != 2 else 3
        e = e - (eta(e, ps[i]) / eta(ps[other], ps[i])) * ps[other]
        eps[i] = e

    return ps, eps, eta


def FdotF(eta, ps, eps, i, j):
    pi, pj, epsi, epsj = ps[i], ps[j], eps[i], eps[j]
    return 2 * (eta(pi, pj) * eta(epsi, epsj) - eta(pi, epsj) * eta(epsi, pj))


def Cabc(eta, ps, eps, k, a, b):
    pk, epsk = ps[k], eps[k]
    return eta(pk, ps[a]) * eta(epsk, ps[b]) - eta(pk, ps[b]) * eta(epsk, ps[a])


def Msandwich(eta, ps, eps, i, j, a, b):
    pi, pj, pa, pb = ps[i], ps[j], ps[a], ps[b]
    ei, ej = eps[i], eps[j]
    Vpj = eta(pa, pi) * eta(ei, pj) - eta(pa, ei) * eta(pi, pj)
    Vej = eta(pa, pi) * eta(ei, ej) - eta(pa, ei) * eta(pi, ej)
    return Vpj * eta(ej, pb) - Vej * eta(pj, pb)


def orbit_sum(eta, ps, eps, base_fn):
    total = 0.0
    for perm in itertools.permutations([1, 2, 3, 4, 5]):
        total += base_fn(eta, ps, eps, *perm)
    return total


SLOTS = ["i", "j", "l", "m", "k"]
ALL_PAIRS_UNORD = list(itertools.combinations(SLOTS, 2))


# familia A (ya usada en la parte 20): 2 F:F + 1 C_{contraccion} * d(par)
def make_seed_A(c_contract, mult_pair):
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        slot_val = {"i": i, "j": j, "l": l, "m": m, "k": k}
        if c_contract == "P":
            c = Cabc(eta, ps, eps, k, i, j)
        elif c_contract == "Q":
            c = Cabc(eta, ps, eps, k, l, m)
        else:
            c = Cabc(eta, ps, eps, k, i, l)
        a, b = slot_val[mult_pair[0]], slot_val[mult_pair[1]]
        mult = eta(ps[a], ps[b])
        return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * c * mult

    return seed_fn


# familia B (nueva, parte 25): 1 F:F(l,m) + 1 sandwich M(i,j;sand_a,sand_b) + 1 C^k_{c_a,c_b}
def make_seed_B(sand_ab, c_ab):
    """sand_ab, c_ab: tuplas de nombres de slot (ordenados para sand_ab, cualquiera para c_ab)."""
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        slot_val = {"i": i, "j": j, "l": l, "m": m, "k": k}
        sa, sb = slot_val[sand_ab[0]], slot_val[sand_ab[1]]
        ca, cb = slot_val[c_ab[0]], slot_val[c_ab[1]]
        if sa == ca == i or sa == i:
            pass  # no-op, solo legibilidad
        M = Msandwich(eta, ps, eps, i, j, sa, sb)
        C = Cabc(eta, ps, eps, k, ca, cb)
        return FdotF(eta, ps, eps, l, m) * M * C

    return seed_fn


if __name__ == "__main__":
    seeds = {}

    # familia A: 3 contracciones x 10 multiplicadores (igual que parte 20)
    for c in ["P", "Q", "M"]:
        for mp in ALL_PAIRS_UNORD:
            seeds[f"A_{c}__d{mp}"] = make_seed_A(c, mp)

    # familia B: sand_ab de {l,m,k} ordenado (6), c_ab de {i,j,l,m} sin k (6 no ordenados)
    remaining_for_sand = ["l", "m", "k"]
    sand_pairs_ordered = [(a, b) for a in remaining_for_sand for b in remaining_for_sand if a != b]
    remaining_for_c = ["i", "j", "l", "m"]
    c_pairs = list(itertools.combinations(remaining_for_c, 2))
    for sand_ab in sand_pairs_ordered:
        for c_ab in c_pairs:
            seeds[f"B_sand{sand_ab}_C{c_ab}"] = make_seed_B(sand_ab, c_ab)

    print(f"Total de semillas: {len(seeds)}  (familia A: 30, familia B: {len(sand_pairs_ordered)*len(c_pairs)})\n")

    seed_names = list(seeds.keys())
    n_samples = 30
    M_mat = np.zeros((n_samples, len(seed_names)))

    for row, kseed in enumerate(range(1, n_samples + 1)):
        ps, eps, eta = gen_kinematics(kseed)
        for col, name in enumerate(seed_names):
            M_mat[row, col] = orbit_sum(eta, ps, eps, seeds[name])

    norms = np.linalg.norm(M_mat, axis=0)
    nonzero_cols = norms > 1e-6
    print(f"Semillas identicamente nulas: {sum(~nonzero_cols)} de {len(seed_names)}")

    Mn = M_mat[:, nonzero_cols] / norms[nonzero_cols]
    sv = np.linalg.svd(Mn, compute_uv=False)
    rank6 = np.linalg.matrix_rank(Mn, tol=1e-6)
    rank9 = np.linalg.matrix_rank(Mn, tol=1e-9)
    print(f"\nRango (tol=1e-6): {rank6}   Rango (tol=1e-9): {rank9}")
    print(f"Semillas no nulas consideradas: {Mn.shape[1]} de {len(seed_names)}, {n_samples} muestras")
    print("\nPrimeros 15 valores singulares:")
    print("  " + ", ".join(f"{v:.4e}" for v in sv[:15]))
