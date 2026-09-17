import itertools, sys, time, pickle, os
import mpmath as mp
from hilbert_series_grado11_altaprecision import gen_kinematics, FdotF, orbit_sum
from grado11_completo_con_sandwich import Msandwich
from topologia_nueva_triple_sandwich import Tsandwich
from grado11_familia_T_triple_sandwich import SLOTS

PKL = "familia_TFFTM_rows.pkl"
all_pairs = list(itertools.combinations(SLOTS, 2))
ijk_pairs = list(itertools.combinations(['i', 'j', 'k'], 2))


def seed_TFF(a_b, mult4):
    def fn(eta, ps, eps, i, j, k, l, m):
        slot = {'i': i, 'j': j, 'k': k, 'l': l, 'm': m}
        a, b = slot[a_b[0]], slot[a_b[1]]
        T = Tsandwich(eta, ps, eps, i, j, k, a, b)
        FF = FdotF(eta, ps, eps, l, m)
        mult = mp.mpf(1)
        for (x, y) in mult4:
            mult *= eta(ps[slot[x]], ps[slot[y]])
        return T * FF * mult
    return fn


def seed_TM(a_b, cd, mult2):
    def fn(eta, ps, eps, i, j, k, l, m):
        slot = {'i': i, 'j': j, 'k': k, 'l': l, 'm': m}
        a, b = slot[a_b[0]], slot[a_b[1]]
        c, d = slot[cd[0]], slot[cd[1]]
        T = Tsandwich(eta, ps, eps, i, j, k, a, b)
        M = Msandwich(eta, ps, eps, l, m, c, d)
        x, y = slot[mult2[0]], slot[mult2[1]]
        mult = eta(ps[x], ps[y])
        return T * M * mult
    return fn


def build_seeds():
    seeds = {}
    mult4_opts = list(itertools.combinations_with_replacement(ijk_pairs, 2))
    for a_b in all_pairs:
        for mult4 in mult4_opts:
            seeds[f"TFF_{a_b}_{mult4}"] = seed_TFF(a_b, mult4)
    cd_opts = ijk_pairs + [('l', 'm')]
    mult2_opts = ijk_pairs + [('l', 'm')]
    for a_b in all_pairs:
        for cd in cd_opts:
            for mult2 in mult2_opts:
                seeds[f"TM_{a_b}_{cd}_{mult2}"] = seed_TM(a_b, cd, mult2)
    return seeds


def load():
    if os.path.exists(PKL):
        with open(PKL, "rb") as f:
            return pickle.load(f)
    return {"seed_names": None, "rows": {}}


def save(data):
    with open(PKL, "wb") as f:
        pickle.dump(data, f)


if __name__ == "__main__":
    kseeds_to_add = [int(x) for x in sys.argv[1:]]
    seeds = build_seeds()
    names = list(seeds.keys())
    data = load()
    if data["seed_names"] is None:
        data["seed_names"] = names
    for kseed in kseeds_to_add:
        if kseed in data["rows"]:
            print(f"kseed {kseed} ya calculado, salteando")
            continue
        t0 = time.time()
        ps, eps, eta = gen_kinematics(kseed)
        row = [orbit_sum(eta, ps, eps, seeds[n]) for n in names]
        data["rows"][kseed] = row
        save(data)
        print(f"kseed {kseed} listo ({time.time()-t0:.1f}s)", flush=True)
