"""
Fase 2, paso 22: generalizacion a n=5 del limite de Regge validado para n=4
(23-regge-limit-n4-validado.md).

Construccion (exacta, sin bisection -- se puede resolver todo analiticamente):

  - Cluster "duro" creciente: p1, p2 (como en n=4), s_12 = 4E^2 -> infinito.
  - p3: tercera particula "dura", con energia E3(E) y angulo theta(E) respecto de n1 ELEGIDOS
    para mantener FIJOS simultaneamente 2 invariantes:
      t13 := -(p1+p3)^2 = t0                 (analogo directo del "t" de n=4)
      M45^2 := -(p1+p2+p3)^2 = m0^2          (masa invariante FIJA del sistema recular {4,5})
    Se deriva analiticamente (ver mas abajo) que esto se logra con:
      E3 = E + m0^2/(4E)
      cos(theta) = 1 - t0/(2 E E3)   ->  1  (limite hacia adelante), igual que en n=4
  - P45 := -(p1+p2+p3) es un vector TIPO TIEMPO de masa invariante FIJA m0 (no depende de E) que
    se aleja del resto del sistema con energia creciente -- se bombea (boost estandar) una
    configuracion de decaimiento de 2 cuerpos FIJA en el frame de reposo de P45 (angulo y reparto
    de energia fijos, eligidos una sola vez) al frame de laboratorio, dando p4(E), p5(E).

Esta construccion es exacta (sin necesidad de bisection numerica): p1,p2,p3 son automaticamente
nulos por construccion (se escriben como E_i*(1,n_i) con |n_i|=1), y p4,p5 son nulos por
construccion del boost de un decaimiento de 2 cuerpos sin masa (p4_rest, p5_rest cada uno con
energia m0/2 en el frame de reposo de P45, cuya norma se preserva bajo boost).

Antes de evaluar cualquier invariante de polarizacion, se CLASIFICAN numericamente todos los pares
(i,j) en "canal creciente" (crece con E) vs "canal fijo" (se mantiene acotado) -- leccion de
23-regge-limit-n4-validado.md: usar un par de canal fijo como referencia de gauge para las
polarizaciones produce un artefacto de crecimiento espurio.
"""

import itertools
import numpy as np

D = 6


def eta(a, b):
    return -a[0] * b[0] + a[1:].dot(b[1:])


def boost_from_rest(P, q_rest):
    """Bombea un 4-vector q_rest (definido en el frame de reposo de P, P=(M,0,...,0)) al frame de
    laboratorio donde el vector de reposo (M,0,...,0) se mapea a P. Formula estandar de boost de
    Lorentz (independiente de la convencion de signo de la metrica)."""
    M = np.sqrt(-eta(P, P))
    P0 = P[0]
    Pvec = P[1:]
    beta_vec = Pvec / P0          # beta * n_hat
    beta2 = beta_vec.dot(beta_vec)
    gamma = P0 / M

    q0, qvec = q_rest[0], q_rest[1:]
    q0p = gamma * (q0 + beta_vec.dot(qvec))
    if beta2 < 1e-30:
        qvecp = qvec.copy()
    else:
        qvecp = qvec + (gamma - 1) * (beta_vec.dot(qvec) / beta2) * beta_vec + gamma * q0 * beta_vec
    return np.concatenate([[q0p], qvecp])


def build_kinematics_n5(E, t0, m0sq, rng):
    n1 = np.zeros(D - 1)
    n1[0] = 1.0

    E3 = E - m0sq / (4 * E)
    cos_theta = 1 + t0 / (2 * E * E3)
    assert -1 <= cos_theta <= 1, f"t0/m0sq no realizables a esta energia (cos_theta={cos_theta})"
    sin_theta = np.sqrt(1 - cos_theta ** 2)

    n3 = np.zeros(D - 1)
    n3[0] = cos_theta
    n3[1] = sin_theta

    p1 = -np.concatenate([[E], E * n1])
    p2 = -np.concatenate([[E], E * (-n1)])
    p3 = np.concatenate([[E3], E3 * n3])

    P45 = -(p1 + p2 + p3)
    M45sq = -eta(P45, P45)
    assert abs(M45sq - m0sq) < 1e-6 * max(1, m0sq, E ** 2 * 1e-9), \
        f"M45^2={M45sq} != m0sq={m0sq}"
    M45 = np.sqrt(m0sq)  # usar el valor exacto conocido -- evita cancelacion catastrofica en floats

    # configuracion de decaimiento de 2 cuerpos FIJA (no depende de E), elegida una sola vez.
    # NOTA (encontrado empiricamente): una direccion n45 completamente ALEATORIA da crecimiento
    # NO limpio (oscilante, sin ley de potencias definida) en los pares cruzados {1,2,3}x{4,5} --
    # diagnosticado como una interferencia entre la escala de boost y la orientacion generica del
    # decaimiento. Alinear n45 EXACTAMENTE con el eje del boost (n3) da crecimiento limpio pero
    # degenera (algunos pares se anulan identicamente -- configuracion de medida cero). La solucion
    # que da crecimiento limpio en TODOS los pares SIN degenerar es una inclinacion FIJA pequena
    # (phi0) respecto del eje del boost.
    phi0 = 0.3
    e_orth = np.zeros(D - 1)
    e_orth[2] = 1.0
    n45 = np.cos(phi0) * n3 + np.sin(phi0) * e_orth
    n45 /= np.linalg.norm(n45)
    p4_rest = np.concatenate([[M45 / 2], (M45 / 2) * n45])
    p5_rest = np.concatenate([[M45 / 2], (M45 / 2) * (-n45)])

    p4 = boost_from_rest(P45, p4_rest)
    p5 = boost_from_rest(P45, p5_rest)

    ps = {1: p1, 2: p2, 3: p3, 4: p4, 5: p5}
    for k in range(1, 6):
        assert abs(eta(ps[k], ps[k])) < 1e-6 * E ** 2, f"p{k}^2 = {eta(ps[k], ps[k])}"
    assert np.linalg.norm(sum(ps.values())) < 1e-6 * E

    return ps


def FdotF(ps, eps, i, j):
    pi, pj, epsi, epsj = ps[i], ps[j], eps[i], eps[j]
    return 2 * (eta(pi, pj) * eta(epsi, epsj) - eta(pi, epsj) * eta(epsi, pj))


def Cabc(ps, eps, k, a, b):
    pk, epsk = ps[k], eps[k]
    return eta(pk, ps[a]) * eta(epsk, ps[b]) - eta(pk, ps[b]) * eta(epsk, ps[a])


if __name__ == "__main__":
    rng = np.random.RandomState(11)
    t0, m0sq = -1.0, 3.0

    print("=" * 100)
    print("PASO 1: clasificar los 10 pares (i,j) en 'canal creciente' vs 'canal fijo'")
    print("=" * 100)
    Es = [50, 500, 5000, 50000]
    pair_vals = {p: [] for p in itertools.combinations(range(1, 6), 2)}
    for E in Es:
        ps = build_kinematics_n5(E, t0, m0sq, rng)
        for (i, j) in pair_vals:
            pair_vals[(i, j)].append(eta(ps[i], ps[j]))

    growing_pairs = []
    fixed_pairs = []
    for (i, j), vals in pair_vals.items():
        ratio = abs(vals[-1] / vals[0])
        slope = np.polyfit(np.log(Es), np.log(np.abs(vals)), 1)[0]
        tag = "CRECIENTE" if slope > 0.5 else "FIJO"
        print(f"  ({i},{j}): eta ~ E^{slope:.3f}   valores={['%.3e' % v for v in vals]}   [{tag}]")
        (growing_pairs if slope > 0.5 else fixed_pairs).append((i, j))

    print(f"\nPares crecientes: {growing_pairs}")
    print(f"Pares fijos: {fixed_pairs}")

    print("\n" + "=" * 100)
    print("PASO 2: construir epsilons usando SOLO pares crecientes como referencia de gauge")
    print("=" * 100)
    # para cada particula, elegir como referencia un companero de un par CRECIENTE
    growing_partner = {}
    for i in range(1, 6):
        for (a, b) in growing_pairs:
            if a == i:
                growing_partner[i] = b
                break
            if b == i:
                growing_partner[i] = a
                break
    print("Referencia de gauge elegida por particula:", growing_partner)
    assert len(growing_partner) == 5, "no se encontro referencia creciente para alguna particula"

    seeds = {i: rng.randn(D) for i in range(1, 6)}

    def make_eps(ps):
        eps = {}
        for i in range(1, 6):
            e = seeds[i].copy()
            other = growing_partner[i]
            e = e - (eta(e, ps[i]) / eta(ps[other], ps[i])) * ps[other]
            eps[i] = e
        return eps

    print("\n" + "=" * 100)
    print("PASO 3: crecimiento del candidato de grado 9 ya confirmado (parte 20)")
    print("=" * 100)

    def T_grado9(ps, eps, i, j, l, m, k):
        d_ik = eta(ps[i], ps[k])
        d_jk = eta(ps[j], ps[k])
        return FdotF(ps, eps, i, j) * FdotF(ps, eps, l, m) * Cabc(ps, eps, k, i, j) * (d_ik - d_jk)

    def orbit_sum(ps, eps, base_fn):
        total = 0.0
        for perm in itertools.permutations([1, 2, 3, 4, 5]):
            total += base_fn(ps, eps, *perm)
        return total

    Es2 = [50, 100, 300, 1000, 3000, 10000, 30000, 100000]
    rows = []
    for E in Es2:
        ps = build_kinematics_n5(E, t0, m0sq, rng)
        eps = make_eps(ps)
        s12 = -eta(ps[1] + ps[2], ps[1] + ps[2])
        val = orbit_sum(ps, eps, T_grado9)
        rows.append((E, s12, val))
        print(f"  E={E:>8}  s12={s12: .6e}  S_grado9={val: .6e}")

    Es_arr = np.array([r[0] for r in rows])
    s12_arr = np.array([r[1] for r in rows])
    val_arr = np.array([abs(r[2]) for r in rows])

    def slope(x, y, k=4):
        lx, ly = np.log(x[-k:]), np.log(y[-k:])
        return np.polyfit(lx, ly, 1)[0]

    print(f"\nS_grado9 ~ E^{slope(Es_arr, val_arr):.4f}   ~ s12^{slope(s12_arr, val_arr):.4f}")
    print("(grado 9 total en momentos -> se esperaria, por conteo dimensional ingenuo sin "
          "cancelaciones, S ~ E^9 si TODOS los factores crecieran como E; pero como algunos pares "
          "son de 'canal fijo', la potencia real puede ser menor -- ver clasificacion del paso 1)")
