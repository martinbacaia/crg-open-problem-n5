"""Extension: familias A+B+C (ya en grado11_completo_con_sandwich.py) + familia D (2 sandwiches)."""
import itertools
from grado11_completo_con_sandwich import (
    mp, gen_kinematics, orbit_sum, mp_matrix_rank, ALL_PAIRS,
    make_seed_A, make_seed_B_v1, make_seed_B_v2, make_seed_C, make_seed_D,
)

if __name__ == "__main__":
    seeds = {}
    mono4 = list(itertools.combinations_with_replacement(ALL_PAIRS, 2))
    for c in ["P", "Q", "M"]:
        for mp_ in mono4:
            seeds[f"A_{c}__d{mp_[0]}*d{mp_[1]}"] = make_seed_A(c, mp_)
    seeds["B_v1"] = make_seed_B_v1()
    seeds["B_v2"] = make_seed_B_v2()

    remaining_for_sand = ["l", "m", "k"]
    sand_pairs_ordered = [(a, b) for a in remaining_for_sand for b in remaining_for_sand if a != b]
    remaining_for_c = ["i", "j", "l", "m"]
    c_pairs = list(itertools.combinations(remaining_for_c, 2))
    mult_pairs_subset = [("l", "m"), ("i", "j"), ("i", "l"), ("k", "l")]
    for sand_ab in sand_pairs_ordered:
        for c_ab in c_pairs:
            for mult_pair in mult_pairs_subset:
                seeds[f"C_sand{sand_ab}_C{c_ab}_d{mult_pair}"] = make_seed_C(sand_ab, c_ab, mult_pair)

    sand1_opts = [("l", "m"), ("l", "k"), ("m", "k")]
    sand2_opts = [("i", "j"), ("i", "k"), ("j", "k")]
    c_opts = [("i", "j"), ("l", "m"), ("i", "l")]
    for s1 in sand1_opts:
        for s2 in sand2_opts:
            for c in c_opts:
                seeds[f"D_{s1}_{s2}_{c}"] = make_seed_D(s1, s2, c)

    print(f"Total de semillas (A+B+C+D): {len(seeds)}")
    seed_names = list(seeds.keys())
    n_samples = 20
    rows = []
    for kseed in range(1, n_samples + 1):
        ps, eps, eta = gen_kinematics(kseed)
        row = [orbit_sum(eta, ps, eps, seeds[name]) for name in seed_names]
        rows.append(row)
        print(f"  muestra {kseed}/{n_samples} calculada")

    norms = [mp.sqrt(sum(rows[r][c] ** 2 for r in range(n_samples))) for c in range(len(seed_names))]
    nz_idx = [c for c, n in enumerate(norms) if n > mp.mpf('1e-25')]
    print(f"Semillas no nulas: {len(nz_idx)} de {len(seed_names)}")
    Mn = [[rows[r][c] / norms[c] for c in nz_idx] for r in range(n_samples)]
    rank = mp_matrix_rank(Mn, verbose=True)
    print(f"\nRango (A+B+C+D, {n_samples} muestras) = {rank}")
