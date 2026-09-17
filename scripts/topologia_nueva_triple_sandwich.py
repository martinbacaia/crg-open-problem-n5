"""
Fase 2 -- topologia GENUINAMENTE NUEVA para grado 11: "triple sandwich"
T(i,j,k;a,b) := p_a . F^i . F^j . F^k . p_b

Razonamiento previo (ver 36-topologia-nueva-triple-sandwich.md): los 3 atomos conocidos
(F^i:F^j grado-momento 2 / 2 slots de polarizacion; C^k_{ab} grado 3 / 1 slot; M^{ij}_{ab}
"sandwich doble" grado 4 / 2 slots) agotan TODAS las formas de partir los 5 slots de polarizacion
en bloques de tamano {1,2} sin exceder grado 11 (demostrado por conteo, no solo por busqueda:
particiones 2+2+1 dan familias A/C/D, particion 2+1+1+1 da familia B; la particion 1+1+1+1+1 ya
excede grado 11). Family completa confirmada (parte 35): dimension 15 exacta, sin margen.

Un atomo NUEVO con 3 slots de polarizacion (grado momento 5) abre nuevas particiones: 3+2 (T+F:F o
T+M, con corrector) y **3+1+1 (T+C+C, SIN corrector, grado exacto 5+3+3=11)** -- la mas prometedora
por analogia con B y D (las familias "exactas", sin corrector, resultaron mas ricas).

Implementacion: se propaga un vector v=p_a a traves de F^i, F^j, F^k aplicando la formula
v.F^m = (v.p_m) eps_m - (v.eps_m) p_m, representando v por sus productos punto con TODOS los p_x,
eps_x (nunca componentes explicitas) -- generaliza directamente el patron ya usado y verificado
para Msandwich (sandwich doble).
"""
import itertools
import mpmath as mp
from hilbert_series_grado11_altaprecision import (
    mp as mpmod, D, gen_kinematics, FdotF, Cabc, orbit_sum, mp_matrix_rank, ALL_PAIRS,
)
from grado11_completo_con_sandwich import Msandwich


def apply_F(vp, ve, m, eta, ps, eps):
    """v -> v.F^m, representado por sus productos punto con todos los p_x, eps_x."""
    new_vp = {}
    new_ve = {}
    for x in range(1, 6):
        new_vp[x] = vp[m] * eta(eps[m], ps[x]) - ve[m] * eta(ps[m], ps[x])
        new_ve[x] = vp[m] * eta(eps[m], eps[x]) - ve[m] * eta(ps[m], eps[x])
    return new_vp, new_ve


def start_vec(a, eta, ps, eps):
    vp = {x: eta(ps[a], ps[x]) for x in range(1, 6)}
    ve = {x: eta(ps[a], eps[x]) for x in range(1, 6)}
    return vp, ve


def Msandwich_via_chain(eta, ps, eps, i, j, a, b):
    vp, ve = start_vec(a, eta, ps, eps)
    vp, ve = apply_F(vp, ve, i, eta, ps, eps)
    vp, ve = apply_F(vp, ve, j, eta, ps, eps)
    return vp[b]


def Tsandwich(eta, ps, eps, i, j, k, a, b):
    vp, ve = start_vec(a, eta, ps, eps)
    vp, ve = apply_F(vp, ve, i, eta, ps, eps)
    vp, ve = apply_F(vp, ve, j, eta, ps, eps)
    vp, ve = apply_F(vp, ve, k, eta, ps, eps)
    return vp[b]


if __name__ == "__main__":
    ps, eps, eta = gen_kinematics(1)

    print("=" * 80)
    print("PASO 0: verificar que la cadena generica reproduce Msandwich ya validado")
    print("=" * 80)
    ok = True
    for i, j, a, b in itertools.product(range(1, 6), repeat=4):
        if len({i, j, a, b}) < 2:
            continue
        v1 = Msandwich(eta, ps, eps, i, j, a, b)
        v2 = Msandwich_via_chain(eta, ps, eps, i, j, a, b)
        if abs(v1 - v2) > mp.mpf('1e-40'):
            print("MISMATCH", i, j, a, b, v1, v2)
            ok = False
    print("Msandwich via cadena generica == Msandwich original:", ok)

    print()
    print("=" * 80)
    print("PASO 1: Tsandwich no es identicamente nulo (chequeo puntual)")
    print("=" * 80)
    val = Tsandwich(eta, ps, eps, 1, 2, 3, 4, 5)
    print("T(1,2,3;4,5) =", val)
    print("no nulo:", abs(val) > mp.mpf('1e-30'))
