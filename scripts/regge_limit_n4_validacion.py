"""
Fase 2, paso 20: construccion EXPLICITA (a nivel de vectores, no solo simbolos s,t,u) del limite de
Regge para n=4 -- "s -> infinito a t fijo" -- como paso de VALIDACION antes de intentar generalizar
a n=5 (ver 22-hacia-regge-growth-n5.md, donde se encontro que una generalizacion ingenua del limite
de cono de luz no funciona por conservacion de momento + imposibilidad de un momento nulo generico
con overlap cero en ambas direcciones de referencia).

Construccion (estandar en la literatura de scattering elastico de alta energia -- 2->2 en el CM
frame, angulo de scattering theta -> 0 correlacionado con la energia E -> infinito de forma que el
momentum transfer t se mantenga EXACTAMENTE fijo):

  p1 = -E(1, n1)      (incoming, convencion all-outgoing: negado)
  p2 = -E(1, n2),  n2 = -n1
  p3 =  E(1, n3)      (outgoing)
  p4 =  E(1, n4),  n4 = -n3

con n1,n3 vectores unitarios espaciales en un plano de scattering de 2 dimensiones (embebido en las
D-1 dimensiones espaciales), angulo entre ellos theta. Formulas estandar (metrica mostly-plus):
  s = 4 E^2                              (crece como Lambda^2 si E=Lambda)
  t = -2 E^2 (1 - cos(theta))            (se mantiene FIJO si cos(theta) = 1 + t0/(2E^2))
  u = -s - t                             (crece como -4E^2, junto con s)

Esta es la realizacion fisica estandar del limite de Regge/alta energia a momentum transfer fijo --
la diferencia clave con el intento ingenuo de la parte 22 es que ACA LAS 4 PARTICULAS PARTICIPAN del
boost de alta energia (ninguna se queda "quieta"), y t se mantiene fijo correlacionando el angulo con
la energia, no dejando 2 particulas sin boostear.

Validacion: se mide numericamente (ajuste log-log de varias escalas de Lambda=E) la tasa de
crecimiento de 2 invariantes construidos con el formalismo F^i:F^j ya usado toda la sesion, y se
compara contra el conteo dimensional ingenuo (grado en momentos == potencia de Lambda, y por lo
tanto potencia de s^(grado/2)):
  - F^1:F^2 (grado 2 en momentos)              -> se espera que crezca como Lambda^2 ~ s^1
  - (F^1:F^2)*(F^3:F^4) (grado 4 en momentos)  -> se espera que crezca como Lambda^4 ~ s^2

(coincide con la tabla del paper, seccion 2.10: estructuras genericas de 2 derivadas crecen como s,
estructuras genericas de 4 derivadas crecen como s^2 -- el paper solo encuentra combinaciones
ESPECIALES, no genericas, que crecen mas lento a mayor orden derivativo; una combinacion generica
sigue el conteo dimensional ingenuo).
"""

import numpy as np

D = 6  # dimension total (1 tiempo + 5 espacio), consistente con el resto del proyecto


def build_kinematics(Lambda, t0, phi_extra=None, rng=None):
    """Cinematica de scattering elastico 2->2 (todo-outgoing) en D=6, con s=4*Lambda^2 y t=t0 FIJO."""
    E = Lambda
    cos_theta = 1 + t0 / (2 * E ** 2)
    assert -1 <= cos_theta <= 1, f"t0 no realizable a esta energia (cos_theta={cos_theta})"
    sin_theta = np.sqrt(1 - cos_theta ** 2)

    n1 = np.zeros(D - 1)
    n1[0] = 1.0
    n2 = -n1

    n3 = np.zeros(D - 1)
    n3[0] = cos_theta
    n3[1] = sin_theta
    n4 = -n3

    def vec(sign, E_, n):
        return sign * np.concatenate([[E_], E_ * n])

    p1 = -vec(1, E, n1)
    p2 = -vec(1, E, n2)
    p3 = vec(1, E, n3)
    p4 = vec(1, E, n4)

    def eta(a, b):
        return -a[0] * b[0] + a[1:].dot(b[1:])

    ps = {1: p1, 2: p2, 3: p3, 4: p4}
    for k in range(1, 5):
        assert abs(eta(ps[k], ps[k])) < 1e-6, f"p{k}^2 = {eta(ps[k],ps[k])}"
    assert np.linalg.norm(sum(ps.values())) < 1e-6

    s = -eta(p1 + p2, p1 + p2)
    t = -eta(p1 + p3, p1 + p3)
    u = -eta(p1 + p4, p1 + p4)

    if rng is None:
        rng = np.random.RandomState(0)
    eps = {}
    for i in range(1, 5):
        e = rng.randn(D)
        other = 2 if i != 2 else 3
        e = e - (eta(e, ps[i]) / eta(ps[other], ps[i])) * ps[other]
        eps[i] = e

    return ps, eps, eta, s, t, u


def FdotF(eta, ps, eps, i, j):
    pi, pj, epsi, epsj = ps[i], ps[j], eps[i], eps[j]
    return 2 * (eta(pi, pj) * eta(epsi, epsj) - eta(pi, epsj) * eta(epsi, pj))


if __name__ == "__main__":
    t0 = -1.0
    rng = np.random.RandomState(7)
    # generar epsilons UNA VEZ como "semillas" (direcciones aleatorias fijas), reproyectadas a cada
    # Lambda -- para que la comparacion entre escalas use "la misma" polarizacion generica
    seeds = {i: rng.randn(D) for i in range(1, 5)}

    # IMPORTANTE (bug encontrado y corregido): la referencia usada para fijar el gauge de cada
    # eps_i NO puede ser un socio del "canal t" (pares (1,3) y (2,4) en esta cinematica), porque
    # eta(p1,p3)=eta(p2,p4)=-t/2 se mantiene FIJO por construccion (es literalmente el momentum
    # transfer que estamos manteniendo constante) -- usarlo como denominador de la proyeccion de
    # gauge hace que la parte transversal de esa polarizacion crezca espuriamente sin limite fisico
    # (se detecto porque (F1:F2)(F3:F4) crecia como Lambda^5 en vez de Lambda^4 esperado). La
    # referencia correcta es el socio del "canal s" (pares (1,2) y (3,4), que SI crecen como
    # Lambda^2), evitando cualquier denominador que se mantenga chico en el limite.
    other_ref = {1: 2, 2: 1, 3: 4, 4: 3}

    Lambdas = [10, 30, 100, 300, 1000, 3000, 10000, 30000]
    rows = []
    for Lam in Lambdas:
        ps, eps_dummy, eta, s, t, u = build_kinematics(Lam, t0)
        eps = {}
        for i in range(1, 5):
            e = seeds[i].copy()
            other = other_ref[i]
            e = e - (eta(e, ps[i]) / eta(ps[other], ps[i])) * ps[other]
            eps[i] = e
        F12 = FdotF(eta, ps, eps, 1, 2)
        F34 = FdotF(eta, ps, eps, 3, 4)
        prod = F12 * F34
        rows.append((Lam, s, t, u, F12, prod))
        print(f"Lambda={Lam:>8.1f}  s={s: .6e}  t={t: .6e}  u={u: .6e}  "
              f"F1:F2={F12: .6e}  (F1:F2)(F3:F4)={prod: .6e}")

    print("\nVerificacion: t se mantiene FIJO en todas las escalas (por construccion) -- OK si "
          "todos los valores de t de arriba son iguales a t0 =", t0)

    Lams = np.array([r[0] for r in rows])
    F12s = np.array([abs(r[4]) for r in rows])
    prods = np.array([abs(r[5]) for r in rows])
    ss = np.array([r[1] for r in rows])

    # ajuste log-log entre las 2 escalas mas grandes (regimen ya asintotico)
    def slope(x, y):
        lx, ly = np.log(x), np.log(y)
        return np.polyfit(lx, ly, 1)[0]

    slope_F12_vs_Lambda = slope(Lams[-4:], F12s[-4:])
    slope_prod_vs_Lambda = slope(Lams[-4:], prods[-4:])
    slope_F12_vs_s = slope(ss[-4:], F12s[-4:])
    slope_prod_vs_s = slope(ss[-4:], prods[-4:])

    print(f"\nF1:F2           ~ Lambda^{slope_F12_vs_Lambda:.4f}   ~ s^{slope_F12_vs_s:.4f}"
          f"   (esperado: Lambda^2 ~ s^1)")
    print(f"(F1:F2)(F3:F4)  ~ Lambda^{slope_prod_vs_Lambda:.4f}   ~ s^{slope_prod_vs_s:.4f}"
          f"   (esperado: Lambda^4 ~ s^2)")
