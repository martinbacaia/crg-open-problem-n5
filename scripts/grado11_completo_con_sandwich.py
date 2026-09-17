"""
Fase 2, paso 27: recalculo RIGUROSO de la dimension de grado 11, incluyendo sistematicamente la
topologia "sandwich" (M^{ij}_{ab}) -- corrige la advertencia de consistencia dejada en
27-grado9-dimension3-y-mejora-a-s3.md: si grado 9 estaba incompleto por no incluir el sandwich,
grado 11 (parte 21, "dimension 3") probablemente tambien lo estaba.

Familias de semillas de grado 11:
  (A) "2 parejas F:F + 1 singlete C" x corrector Mandelstam grado 4 -- ya usada en la parte 21
      (165 semillas: 3 contracciones del singlete x 55 monomios de grado 4).
  (B) "1 pareja F:F + 3 singletes C" -- ya usada en la parte 21 (2 variantes representativas, NO
      exhaustiva).
  (C) NUEVA: "1 pareja F:F + 1 sandwich + 1 singlete C" x corrector Mandelstam grado 2
      (grado 2+4+3+2=11) -- analogo directo de la familia A de grado 9 (parte 26), con un
      corrector extra de grado 2 para llegar a 11.
  (D) NUEVA: "2 sandwiches + 1 singlete C" (grado 4+4+3=11, SIN corrector) -- topologia que no
      tiene analogo en grado 9 (ahi no entraban 2 sandwiches).

Se usa ALTA PRECISION (mpmath, 60 digitos) desde el principio, dada la leccion de la parte 18b/21:
float64 no es confiable para determinar rango en grado 11.
"""

import itertools
from hilbert_series_grado11_altaprecision import (
    mp, D, gen_kinematics, FdotF, Cabc, orbit_sum, mp_matrix_rank,
    ALL_PAIRS, SLOT_NAMES,
)


def Msandwich(eta, ps, eps, i, j, a, b):
    pi, pj, pa, pb = ps[i], ps[j], ps[a], ps[b]
    ei, ej = eps[i], eps[j]
    Vpj = eta(pa, pi) * eta(ei, pj) - eta(pa, ei) * eta(pi, pj)
    Vej = eta(pa, pi) * eta(ei, ej) - eta(pa, ei) * eta(pi, ej)
    return Vpj * eta(ej, pb) - Vej * eta(pj, pb)


# familia A (ya usada en la parte 21)
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


# familia C (NUEVA): F^l:F^m * M^{ij}_{sand_a,sand_b} * C^k_{c_a,c_b} * d(mult_pair)
def make_seed_C(sand_ab, c_ab, mult_pair):
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        slot_val = {"i": i, "j": j, "l": l, "m": m, "k": k}
        sa, sb = slot_val[sand_ab[0]], slot_val[sand_ab[1]]
        ca, cb = slot_val[c_ab[0]], slot_val[c_ab[1]]
        M = Msandwich(eta, ps, eps, i, j, sa, sb)
        C = Cabc(eta, ps, eps, k, ca, cb)
        mult_a, mult_b = slot_val[mult_pair[0]], slot_val[mult_pair[1]]
        mult = eta(ps[mult_a], ps[mult_b])
        return FdotF(eta, ps, eps, l, m) * M * C * mult

    return seed_fn


# familia D (NUEVA): M^{ij}_{a1,b1} * M^{lm}_{a2,b2} * C^k_{c1,c2}  (grado 4+4+3=11, sin corrector)
def make_seed_D(sand1_ab, sand2_ab, c_ab):
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        slot_val = {"i": i, "j": j, "l": l, "m": m, "k": k}
        a1, b1 = slot_val[sand1_ab[0]], slot_val[sand1_ab[1]]
        a2, b2 = slot_val[sand2_ab[0]], slot_val[sand2_ab[1]]
        ca, cb = slot_val[c_ab[0]], slot_val[c_ab[1]]
        M1 = Msandwich(eta, ps, eps, i, j, a1, b1)
        M2 = Msandwich(eta, ps, eps, l, m, a2, b2)
        C = Cabc(eta, ps, eps, k, ca, cb)
        return M1 * M2 * C

    return seed_fn


if __name__ == "__main__":
    seeds = {}

    # familia A: 55 monomios de grado 4 x 3 contracciones
    mono4 = list(itertools.combinations_with_replacement(ALL_PAIRS, 2))
    for c in ["P", "Q", "M"]:
        for mp_ in mono4:
            seeds[f"A_{c}__d{mp_[0]}*d{mp_[1]}"] = make_seed_A(c, mp_)

    seeds["B_v1"] = make_seed_B_v1()
    seeds["B_v2"] = make_seed_B_v2()

    # familia C: sandwich contraido con {l,m,k} (6 ordenados), C^k contraido con {i,j,l,m} (6),
    # corrector de grado 2 (10 pares) -- se restringe a un subconjunto razonable para no explotar
    # combinatoriamente: se usan los 2 tipos de contraccion del sandwich ya encontrados utiles en
    # grado 9 (con l/m y con l/k), los 6 posibles C^k, y los 10 correctores.
    remaining_for_sand = ["l", "m", "k"]
    sand_pairs_ordered = [(a, b) for a in remaining_for_sand for b in remaining_for_sand if a != b]
    remaining_for_c = ["i", "j", "l", "m"]
    c_pairs = list(itertools.combinations(remaining_for_c, 2))
    # subconjunto razonable de correctores (no los 10, para mantener el computo tratable en alta
    # precision) -- se eligen 4 que cubren los distintos "tipos" (dentro del par F:F, dentro del
    # par sandwich, mixto, y con el singlete)
    mult_pairs_subset = [("l", "m"), ("i", "j"), ("i", "l"), ("k", "l")]
    for sand_ab in sand_pairs_ordered:
        for c_ab in c_pairs:
            for mult_pair in mult_pairs_subset:
                seeds[f"C_sand{sand_ab}_C{c_ab}_d{mult_pair}"] = make_seed_C(sand_ab, c_ab, mult_pair)

    print(f"Total de semillas (A+B+C, D se agrega aparte por volumen): {len(seeds)}")

    seed_names = list(seeds.keys())
    n_samples = 15
    rows = []
    for kseed in range(1, n_samples + 1):
        ps, eps, eta = gen_kinematics(kseed)
        row = [orbit_sum(eta, ps, eps, seeds[name]) for name in seed_names]
        rows.append(row)
        print(f"  muestra {kseed}/{n_samples} calculada")

    norms = []
    for c in range(len(seed_names)):
        norms.append(mp.sqrt(sum(rows[r][c] ** 2 for r in range(n_samples))))
    nz_idx = [c for c, nrm in enumerate(norms) if nrm > mp.mpf('1e-25')]
    print(f"\nSemillas no nulas: {len(nz_idx)} de {len(seed_names)}")

    Mn = [[rows[r][c] / norms[c] for c in nz_idx] for r in range(n_samples)]
    rank = mp_matrix_rank(Mn, verbose=True)
    print(f"\nRango (A+B+C, alta precision, {n_samples} muestras) = {rank}")
