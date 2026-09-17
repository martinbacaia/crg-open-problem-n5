"""
Fase 2, ampliacion de la parte 28: en vez de los subconjuntos representativos usados ahi para las
familias C y D (por miedo al costo computacional), se prueba la COMBINATORIA COMPLETA:

- Familia C: sandwich contraido con las 6 combinaciones ordenadas de {l,m,k} (ya era completo) x
  C^k contraido con las 6 combinaciones de {i,j,l,m} (ya era completo) x TODOS los 10 correctores
  posibles de grado 2 (antes solo 4 representativos) = 360 semillas (antes 144).
- Familia D: sandwich 1 con TODAS las 6 combinaciones ordenadas de {l,m,k} (antes 3) x sandwich 2
  con TODAS las 6 de {i,j,k} (antes 3) x C^k con TODOS los 10 pares posibles de {i,j,l,m,k} (antes
  3) = 360 semillas (antes 27).

Total: A(165)+B(2)+C(360)+D(360) = 887 semillas (antes 338). Estimado por timing: ~330s para 20
muestras -- se corre en 2 lotes independientes (semillas 1-20 y 201-220) para verificar
estabilidad del rango, como en el resto del proyecto.
"""
import itertools, time, sys
from grado11_con_familia_D import (
    mp, gen_kinematics, orbit_sum, mp_matrix_rank, ALL_PAIRS,
    make_seed_A, make_seed_B_v1, make_seed_B_v2, make_seed_C, make_seed_D,
)

def build_seeds():
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
    ALL5 = ["i", "j", "l", "m", "k"]
    mult_pairs_full = list(itertools.combinations(ALL5, 2))  # 10, antes 4
    for sand_ab in sand_pairs_ordered:
        for c_ab in c_pairs:
            for mult_pair in mult_pairs_full:
                seeds[f"C_sand{sand_ab}_C{c_ab}_d{mult_pair}"] = make_seed_C(sand_ab, c_ab, mult_pair)

    sand1_full = [(a, b) for a in remaining_for_sand for b in remaining_for_sand if a != b]  # 6, antes 3
    remaining_for_sand2 = ["i", "j", "k"]
    sand2_full = [(a, b) for a in remaining_for_sand2 for b in remaining_for_sand2 if a != b]  # 6, antes 3
    c_full = list(itertools.combinations(ALL5, 2))  # 10, antes 3
    for s1 in sand1_full:
        for s2 in sand2_full:
            for c in c_full:
                seeds[f"D_{s1}_{s2}_{c}"] = make_seed_D(s1, s2, c)

    return seeds


if __name__ == "__main__":
    seeds = build_seeds()
    print(f"Total de semillas (A+B+C+D, combinatoria COMPLETA): {len(seeds)}")
    seed_names = list(seeds.keys())

    batch_start = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n_samples = int(sys.argv[2]) if len(sys.argv) > 2 else 20

    t0 = time.time()
    rows = []
    for kseed in range(batch_start, batch_start + n_samples):
        ts = time.time()
        ps, eps, eta = gen_kinematics(kseed)
        row = [orbit_sum(eta, ps, eps, seeds[name]) for name in seed_names]
        rows.append(row)
        print(f"  muestra {kseed} calculada ({time.time()-ts:.1f}s, acumulado {time.time()-t0:.1f}s)", flush=True)
    print(f"tiempo total: {time.time()-t0:.1f}s", flush=True)

    norms = [mp.sqrt(sum(rows[r][c] ** 2 for r in range(n_samples))) for c in range(len(seed_names))]
    nz_idx = [c for c, n in enumerate(norms) if n > mp.mpf('1e-25')]
    print(f"Semillas no nulas: {len(nz_idx)} de {len(seed_names)}")
    Mn = [[rows[r][c] / norms[c] for c in nz_idx] for r in range(n_samples)]
    rank = mp_matrix_rank(Mn, verbose=True)
    print(f"\nRango (combinatoria completa, semillas {batch_start}-{batch_start+n_samples-1}) = {rank}")

    import pickle
    with open(f"grado11_combinatoria_completa_batch{batch_start}.pkl", "wb") as f:
        pickle.dump({"seed_names": seed_names, "nz_idx": nz_idx, "rank": rank}, f)
