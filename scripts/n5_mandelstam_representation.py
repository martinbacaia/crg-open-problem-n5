"""
Fase 2, paso 2: derivar (no solo citar) como actua el grupo de permutaciones sobre los invariantes
de Mandelstam independientes para SCATTERING DE 5 PARTICULAS, que es la pieza que falta para
extender la clasificacion de S-matrices de n=4 (arXiv:1910.14392) a n=5 (el hueco de la pregunta #44
de Strings 2024 -- ver fase1-mapa/04-profundizacion-crg.md).

Motivacion del calculo
-----------------------
Para n=4 particulas sin masa, hay C(4,2)=6 invariantes s_ij=(p_i+p_j)^2, pero la conservacion de
momento (suma p_i = 0) impone 4 relaciones lineales, dejando solo 2 invariantes independientes
(s,t, con u=-s-t) -- el paper usa esto y encuentra que S4 actua sobre ese espacio de 2 invariantes
via el cociente S4/(Z2xZ2) = S3 (su seccion 2.3).

Para n=5, hay C(5,2)=10 invariantes s_ij, y la conservacion de momento impone 5 relaciones (una por
particula), dejando 10-5=5 invariantes independientes. La pregunta que el paper NO responde (porque
solo trata n=4) es: como actua S5 sobre ese espacio de 5 invariantes independientes?

Enfoque (derivado desde cero, verificado por 2 metodos):
  1. El espacio "crudo" de los 10 invariantes s_ij es la representacion de permutacion de S5 sobre
     los 2-subconjuntos de {1,...,5} (los "bordes" del grafo completo K5).
  2. La conservacion de momento, para cada particula k, dice sum_{j != k} s_kj = 0. Esto es
     exactamente el mapa de incidencia (borde -> vertice) del grafo K5, que es S5-equivariante.
  3. El espacio de invariantes independientes = NUCLEO de ese mapa de incidencia. Por ser el nucleo
     de un mapa equivariante, es un sub-espacio invariante bajo S5 -- es decir, es UNA
     REPRESENTACION genuina de S5, no solo un espacio vectorial.
  4. Se calcula esa representacion explicitamente (matrices 5x5 para cada permutacion), se verifica
     que es irreducible, y se identifica su caracter.

Resultado nuevo de esta sesion (no esta en el paper original, que solo trata n=4): la representacion
de S5 sobre los invariantes de Mandelstam de 5 puntos es IRREDUCIBLE de dimension 5 -- a diferencia
de n=4, donde el espacio de invariantes (s,t) es una representacion de dimension 2 de S4 (que resulta
ser tambien irreducible, la representacion [2,2] de S4, cuya restriccion al subgrupo S3 da la "2M"
usada en el paper). Esto es estructuralmente MUY distinto de "S3 actuando por simple permutacion de
(s,t,u) entre si" -- para n=5 no hay una base natural donde S5 simplemente permute los invariantes.
"""

from fractions import Fraction
from itertools import combinations, permutations

import sympy as sp


def edges(n):
    return list(combinations(range(1, n + 1), 2))


def incidence_matrix(n):
    """Matriz de incidencia vertices x aristas de K_n (entradas 0/1)."""
    edge_list = edges(n)
    M = sp.zeros(n, len(edge_list))
    for col, (i, j) in enumerate(edge_list):
        M[i - 1, col] = 1
        M[j - 1, col] = 1
    return M, edge_list


def permuted_edge_index(edge_list, sigma, edge):
    i, j = edge
    new_edge = tuple(sorted((sigma[i - 1], sigma[j - 1])))
    return edge_list.index(new_edge)


def edge_permutation_matrix(n, edge_list, sigma):
    """Matriz de permutacion (dim = num aristas) inducida por sigma sobre las aristas de K_n."""
    m = len(edge_list)
    P = sp.zeros(m, m)
    for col, e in enumerate(edge_list):
        row = permuted_edge_index(edge_list, sigma, e)
        P[row, col] = 1
    return P


def kernel_basis_matrix(n):
    """Base (como columnas) del nucleo del mapa de incidencia -- el espacio de invariantes
    de Mandelstam independientes, como combinaciones lineales de los s_ij crudos."""
    M, edge_list = incidence_matrix(n)
    ns = M.nullspace()
    dim_expected = n * (n - 3) // 2
    assert len(ns) == dim_expected, f"dimension del nucleo = {len(ns)}, esperada {dim_expected}"
    basis = sp.Matrix.hstack(*ns)
    return basis, edge_list


def action_on_kernel(n, basis, edge_list, sigma):
    """
    Matriz (dim_ker x dim_ker) que representa la accion de sigma sobre el espacio de invariantes
    independientes, expresada en la base `basis` del nucleo.
    """
    P = edge_permutation_matrix(n, edge_list, sigma)
    dim_ker = basis.shape[1]
    cols = []
    for c in range(dim_ker):
        v = basis[:, c]
        pv = P * v
        # expresar pv como combinacion lineal de las columnas de basis: resolver basis @ x = pv
        sol, params = basis.gauss_jordan_solve(pv)
        assert not params, "el nucleo no es invariante bajo esta permutacion (no deberia pasar)"
        cols.append(sol)
    return sp.Matrix.hstack(*cols)


def cycle_type(sigma):
    n = len(sigma)
    seen = [False] * n
    cycles = []
    for i in range(n):
        if seen[i]:
            continue
        length = 0
        j = i
        while not seen[j]:
            seen[j] = True
            j = sigma[j] - 1
            length += 1
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


def sign_of(sigma):
    n = len(sigma)
    ct = cycle_type(sigma)
    return (-1) ** sum(c - 1 for c in ct)


def class_representatives(n):
    """Un representante de cada clase de conjugacion de S_n, mas el tamano de la clase."""
    from sympy.utilities.iterables import partitions as sympy_partitions
    from math import factorial

    reps = {}
    class_sizes = {}
    for sigma in permutations(range(1, n + 1)):
        ct = cycle_type(sigma)
        if ct not in reps:
            reps[ct] = sigma
            class_sizes[ct] = 0
        class_sizes[ct] += 1
    return reps, class_sizes


def hand_derived_character(n, ct):
    """Formula derivada a mano: chi_kernel(g) = C(fix,2) + c2(g) - fix, donde fix = numero de
    puntos fijos y c2 = numero de 2-ciclos en el tipo de ciclo ct."""
    fix = ct.count(1)
    c2 = ct.count(2)
    return sp.Rational(fix * (fix - 1), 2) + c2 - fix


if __name__ == "__main__":
    for n in (4, 5):
        print("=" * 78)
        print(f"n = {n} particulas: espacio de invariantes de Mandelstam independientes")
        print("=" * 78)

        basis, edge_list = kernel_basis_matrix(n)
        dim_ker = basis.shape[1]
        print(f"Numero de invariantes 'crudos' s_ij: {len(edge_list)}  "
              f"(esperado C(n,2)={n*(n-1)//2})")
        print(f"Dimension del nucleo (invariantes independientes): {dim_ker}  "
              f"(esperado n(n-3)/2={n*(n-3)//2})")

        reps, class_sizes = class_representatives(n)
        print(f"\nClases de conjugacion de S{n}: {len(reps)}")

        print(f"\n{'tipo de ciclo':>18} | {'tamano clase':>12} | {'traza (matriz)':>14} | "
              f"{'formula a mano':>14} | match")
        traces = {}
        all_match = True
        for ct, sigma in sorted(reps.items(), key=lambda kv: (-len(kv[0]), kv[0])):
            A = action_on_kernel(n, basis, edge_list, sigma)
            tr_matrix = A.trace()
            tr_formula = hand_derived_character(n, ct)
            ok = sp.simplify(tr_matrix - tr_formula) == 0
            all_match = all_match and ok
            traces[ct] = tr_matrix
            print(f"{str(ct):>18} | {class_sizes[ct]:>12} | {str(tr_matrix):>14} | "
                  f"{str(tr_formula):>14} | {'OK' if ok else 'MISMATCH'}")

        # norma del caracter: <chi,chi> = (1/n!) sum_g |chi(g)|^2 = (1/n!) sum_clases size*chi^2
        from math import factorial
        norm = sum(class_sizes[ct] * traces[ct] ** 2 for ct in reps) / factorial(n)
        print(f"\n||chi||^2 = {norm}  (debe ser 1 si la representacion es irreducible)")
        print(f"Traza en la identidad (= dimension) = {traces[tuple([1]*n)]}")
        print(f"Todas las trazas coinciden con la formula a mano: {all_match}")

        if n == 4:
            print("\nNota: para n=4 esta representacion de dimension 2 de S4 es la que, restringida")
            print("al subgrupo S3 (fijando la particula 4), debe coincidir con la '2M' ya validada")
            print("en s3_partition_functions.py / verificacion_rigurosa.py.")
        print()
