import pickle, sympy as sp
from fractions import Fraction

with open("S116_familia1.pkl", "rb") as f:
    totals = pickle.load(f)  # LSeries objects, S1..S16

orders = list(range(-10, -4))  # -10..-5, mismo umbral que grado11_extraer_base9.py (bound: >= -4 = s^2)
rows = []
for o in orders:
    rows.append([sp.Rational(t.coeff(o).numerator, t.coeff(o).denominator) for t in totals])
M = sp.Matrix(rows)
print("shape", M.shape, "rank", M.rank())
ns = M.nullspace()
print("nullspace dim (familia 1, con S16):", len(ns))
for i, v in enumerate(ns):
    print(f"U{i+1}:", list(v.T))

with open("grado11_S16_nullspace_familia1.pkl", "wb") as f:
    pickle.dump(ns, f)
