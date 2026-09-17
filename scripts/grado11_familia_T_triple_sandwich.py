"""
Familia T (NUEVA, ver topologia_nueva_triple_sandwich.py): T(i,j,k;a,b) * C(l;c,d) * C(m;e,f),
grado exacto 5+3+3=11, sin corrector.

Se combina con las filas ya calculadas de A+B+C+D (combinatoria_completa_rows.pkl, mismas 20
muestras cinematicas de gen_kinematics(1..20)) para ver si el rango total sube mas alla de 15.
"""
import itertools, pickle, time
import mpmath as mp
from hilbert_series_grado11_altaprecision import gen_kinematics, Cabc, orbit_sum, mp_matrix_rank
from topologia_nueva_triple_sandwich import Tsandwich

SLOTS = ['i', 'j', 'k', 'l', 'm']
ALL_PAIRS5 = list(itertools.combinations(SLOTS, 2))


def make_seed_T(a_b, c1_ab, c2_ab):
    def seed_fn(eta, ps, eps, i, j, k, l, m):
        slot = {'i': i, 'j': j, 'k': k, 'l': l, 'm': m}
        a, b = slot[a_b[0]], slot[a_b[1]]
        T = Tsandwich(eta, ps, eps, i, j, k, a, b)
        c1a, c1b = slot[c1_ab[0]], slot[c1_ab[1]]
        c2a, c2b = slot[c2_ab[0]], slot[c2_ab[1]]
        C1 = Cabc(eta, ps, eps, l, c1a, c1b)
        C2 = Cabc(eta, ps, eps, m, c2a, c2b)
        return T * C1 * C2
    return seed_fn


def build_seeds_T():
    seeds = {}
    remaining_for_c = ['i', 'j', 'k']
    c_pairs = list(itertools.combinations(remaining_for_c, 2))  # 3 pairs: ij, ik, jk
    for a_b in ALL_PAIRS5:  # 10
        for c1_ab in c_pairs:  # 3
            for c2_ab in c_pairs:  # 3
                seeds[f"T_{a_b}_{c1_ab}_{c2_ab}"] = make_seed_T(a_b, c1_ab, c2_ab)
    return seeds


if __name__ == "__main__":
    with open("combinatoria_completa_rows.pkl", "rb") as f:
        old_data = pickle.load(f)
    old_seed_names = old_data["seed_names"]
    old_rows = old_data["rows"]
    n_samples = len(old_rows)
    kseeds = sorted(old_rows.keys())
    print(f"reusando {n_samples} muestras ya calculadas (A+B+C+D): {kseeds}")

    seeds_T = build_seeds_T()
    T_names = list(seeds_T.keys())
    print(f"familia T: {len(T_names)} semillas nuevas")

    t0 = time.time()
    new_rows = {}
    for kseed in kseeds:
        ps, eps, eta = gen_kinematics(kseed)
        new_rows[kseed] = [orbit_sum(eta, ps, eps, seeds_T[n]) for n in T_names]
        print(f"  muestra {kseed} (familia T) calculada, acumulado {time.time()-t0:.1f}s", flush=True)

    with open("familia_T_rows.pkl", "wb") as f:
        pickle.dump({"seed_names": T_names, "rows": new_rows}, f)

    all_names = old_seed_names + T_names
    rows_combined = []
    for k in kseeds:
        rows_combined.append(old_rows[k] + new_rows[k])

    norms = [mp.sqrt(sum(rows_combined[r][c] ** 2 for r in range(n_samples))) for c in range(len(all_names))]
    nz_idx = [c for c, n in enumerate(norms) if n > mp.mpf('1e-25')]
    print(f"semillas no nulas: {len(nz_idx)} de {len(all_names)}")
    Mn = [[rows_combined[r][c] / norms[c] for c in nz_idx] for r in range(n_samples)]
    rank = mp_matrix_rank(Mn, verbose=True)
    print(f"\nRANGO COMBINADO (A+B+C+D + T, {n_samples} muestras) = {rank}")
    print(f"(comparar contra 15, que era el rango sin la familia T)")
