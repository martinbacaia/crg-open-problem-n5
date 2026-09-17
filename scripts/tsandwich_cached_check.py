"""
Verifica que una version de Tsandwich basada en CACHES eta_pp/eta_pe/eta_ee (analoga al patron
usado para Msandwich en regge_n5_simbolico_grado11.py / fast_grado11_familia*.py) coincide
exactamente con la version de referencia (cadena generica sobre vectores explicitos,
topologia_nueva_triple_sandwich.py, ya validada en la parte 36).

Necesario antes de meter T en el pipeline de Regge (que usa eta_pp/pe/ee cacheados, no vectores).
"""
import itertools
import mpmath as mp
from hilbert_series_grado11_altaprecision import gen_kinematics
from topologia_nueva_triple_sandwich import Tsandwich as Tsandwich_ref

P = list(range(1, 6))


def Tsandwich_cached(eta_pp, eta_pe, eta_ee, i, j, k, a, b):
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


if __name__ == "__main__":
    ok = True
    for seed in (1, 2, 3):
        ps, eps, eta = gen_kinematics(seed)
        eta_pp, eta_pe, eta_ee = {}, {}, {}
        for a in P:
            for b in P:
                eta_pp[(a, b)] = eta(ps[a], ps[b])
                eta_pe[(a, b)] = eta(ps[a], eps[b])
                eta_ee[(a, b)] = eta(eps[a], eps[b])

        n_checked = 0
        for i, j, k, a, b in itertools.product(P, repeat=5):
            if len({i, j, k, a, b}) < 3:
                continue
            v_ref = Tsandwich_ref(eta, ps, eps, i, j, k, a, b)
            v_cached = Tsandwich_cached(eta_pp, eta_pe, eta_ee, i, j, k, a, b)
            if abs(v_ref - v_cached) > mp.mpf('1e-40'):
                print("MISMATCH", seed, i, j, k, a, b, v_ref, v_cached)
                ok = False
            n_checked += 1
        print(f"seed {seed}: {n_checked} tuplas chequeadas, todas coinciden: {ok}")
    print("RESULTADO FINAL:", "OK" if ok else "FALLO")
