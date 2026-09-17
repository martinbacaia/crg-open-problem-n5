"""
Test de correctitud de lseries.py contra sympy (referencia ya verificada en el resto del proyecto).
Ops MODERADAS a proposito (nada de encadenar pow3 sobre pow3 -- eso hace explotar el TAMANO de los
enteros exactos en AMBOS lados por igual, no es un problema de lseries, pero vuelve el test
innecesariamente lento sin agregar poder de deteccion de bugs). Lo que importa es cubrir bien
+,-,*,**2,**3 y escalares, con ordenes negativos Y positivos mezclados (el caso real: eta_pp, F:F,
etc. tienen ordenes negativos, y algo asi rompio la v1 de lseries -- ver comentario en el modulo).
"""
import random, sympy as sp
from fractions import Fraction
from lseries import LSeries

x = sp.Symbol('x')

def rand_sympy_poly(seed, lo=-5, hi=5):
    rng = random.Random(seed)
    return sum(sp.Rational(rng.randint(-9,9), rng.randint(1,5)) * x**k for k in range(lo,hi))

def sympy_coeffs(expr, lo, hi):
    expr = sp.expand(expr)
    d = {}
    for o in range(lo, hi):
        c = expr.coeff(x, o)
        if c != 0:
            num, den = sp.together(c).as_numer_denom()
            d[o] = Fraction(int(num), int(den))
    return d

def main():
    n_trials = 100
    ok = True
    for trial in range(n_trials):
        rng = random.Random(2000+trial)
        lo, hi = -5, 5
        a_sp = rand_sympy_poly(2*trial, lo, hi)
        b_sp = rand_sympy_poly(2*trial+1, lo, hi)
        c_sp = rand_sympy_poly(2*trial+7919, lo, hi)
        a_ls = LSeries.from_sympy(a_sp, x, lo-1, hi+1)
        b_ls = LSeries.from_sympy(b_sp, x, lo-1, hi+1)
        c_ls = LSeries.from_sympy(c_sp, x, lo-1, hi+1)

        # una secuencia FIJA mezclada de operaciones (representativa de las formulas reales:
        # productos de pocos factores, sumas, una potencia chica, escalares), no encadenamiento
        # patologico de potencias
        op = rng.choice(['A','B','C','D','E'])
        if op == 'A':   # tipo FdotF: 2*(a*b - c*a)
            r_sp = sp.expand(2*(a_sp*b_sp - c_sp*a_sp)); r_ls = (a_ls*b_ls - c_ls*a_ls)*2
        elif op == 'B':  # tipo Cabc: a*b - a*c
            r_sp = sp.expand(a_sp*b_sp - a_sp*c_sp); r_ls = a_ls*b_ls - a_ls*c_ls
        elif op == 'C':  # producto de 4 factores (tipo S_i)
            r_sp = sp.expand(a_sp*b_sp*c_sp*a_sp); r_ls = a_ls*b_ls*c_ls*a_ls
        elif op == 'D':  # con **2 y resta, mas escalar
            s = rng.randint(-4,4)
            r_sp = sp.expand(s*(a_sp**2 - b_sp*c_sp)); r_ls = (a_ls**2 - b_ls*c_ls)*s
        elif op == 'E':  # suma de varios productos (tipo orbit sum de pocos terminos)
            r_sp = sp.expand(a_sp*b_sp + b_sp*c_sp - c_sp*a_sp)
            r_ls = a_ls*b_ls + b_ls*c_ls - c_ls*a_ls

        a = sympy_coeffs(r_sp, -60, 60)
        b = {o: v for o, v in r_ls.c.items() if -60 <= o < 60}
        if a != b:
            print(f"MISMATCH trial{trial} op={op}")
            print("  sympy-only/diff:", {k: (a.get(k), b.get(k)) for k in set(a)|set(b) if a.get(k)!=b.get(k)})
            ok = False

    print("ALL OK" if ok else "FAILURES FOUND")
    return ok

if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
