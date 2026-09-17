"""
Re-verificacion independiente (sesion de revision adversarial, 2026-09-17) de la dimension de
grado 9 con ALTA PRECISION (mpmath 60 digitos), en vez de float64 como en
grado9_completo_con_sandwich.py -- por las dudas, dado que grado 9 esta cerca del umbral (grado
~10) donde el proyecto ya encontro que float64 deja de ser confiable (ver
21-teorema-paridad-grados-8-10-11.md). Reusa el generador de cinematica ya validado de
hilbert_series_grado11_altaprecision.py y las mismas 66 semillas (familias A+B) de
grado9_completo_con_sandwich.py.
"""
import itertools
import mpmath as mp
from hilbert_series_grado11_altaprecision import gen_kinematics, mp_matrix_rank, FdotF, Cabc, orbit_sum

mp.mp.dps = 60

SLOTS = ["i", "j", "l", "m", "k"]
ALL_PAIRS_UNORD = list(itertools.combinations(SLOTS, 2))


def Msandwich(eta, ps, eps, i, j, a, b):
    pi, pj, pa, pb = ps[i], ps[j], ps[a], ps[b]
    ei, ej = eps[i], eps[j]
    Vpj = eta(pa, pi) * eta(ei, pj) - eta(pa, ei) * eta(pi, pj)
    Vej = eta(pa, pi) * eta(ei, ej) - eta(pa, ei) * eta(pi, ej)
    return Vpj * eta(ej, pb) - Vej * eta(pj, pb)


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


def make_seed_B(sand_ab, c_ab):
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        slot_val = {"i": i, "j": j, "l": l, "m": m, "k": k}
        sa, sb = slot_val[sand_ab[0]], slot_val[sand_ab[1]]
        ca, cb = slot_val[c_ab[0]], slot_val[c_ab[1]]
        M = Msandwich(eta, ps, eps, i, j, sa, sb)
        C = Cabc(eta, ps, eps, k, ca, cb)
        return FdotF(eta, ps, eps, l, m) * M * C
    return seed_fn


def build_seeds():
    seeds = {}
    for c in ["P", "Q", "M"]:
        for mp_ in ALL_PAIRS_UNORD:
            seeds[f"A_{c}__d{mp_}"] = make_seed_A(c, mp_)
    remaining_for_sand = ["l", "m", "k"]
    sand_pairs_ordered = [(a, b) for a in remaining_for_sand for b in remaining_for_sand if a != b]
    remaining_for_c = ["i", "j", "l", "m"]
    c_pairs = list(itertools.combinations(remaining_for_c, 2))
    for sand_ab in sand_pairs_ordered:
        for c_ab in c_pairs:
            seeds[f"B_sand{sand_ab}_C{c_ab}"] = make_seed_B(sand_ab, c_ab)
    return seeds


if __name__ == "__main__":
    seeds = build_seeds()
    seed_names = list(seeds.keys())
    print(f"Total de semillas: {len(seed_names)} (familia A: 30, familia B: {len(seed_names)-30})")

    for batch_name, seed_range in [("lote1 (1-20)", range(1, 21)), ("lote2 (301-320)", range(301, 321))]:
        rows = []
        for kseed in seed_range:
            ps, eps, eta = gen_kinematics(kseed)
            rows.append([orbit_sum(eta, ps, eps, seeds[name]) for name in seed_names])
        n_samples = len(rows)
        norms = [mp.sqrt(sum(rows[r][c] ** 2 for r in range(n_samples))) for c in range(len(seed_names))]
        nz_idx = [c for c, nn in enumerate(norms) if nn > mp.mpf('1e-25')]
        print(f"\n{batch_name}: semillas no nulas: {len(nz_idx)} de {len(seed_names)}")
        Mn = [[rows[r][c] / norms[c] for c in nz_idx] for r in range(n_samples)]
        rank = mp_matrix_rank(Mn, verbose=True)
        print(f"RANGO ({batch_name}, {mp.mp.dps} digitos, {n_samples} muestras) = {rank}")
