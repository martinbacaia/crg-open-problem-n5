import pickle, sympy as sp
import regge_n5_simbolico_grado11 as G

x = G.x

with open("S115.pkl","rb") as f:
    totals = pickle.load(f)

# orders x^-10 .. x^-5 (6 rows) -> null space of dim 9 (per 30-derivacion-simbolica-cascada.md)
orders = list(range(-10,-4))  # -10,-9,-8,-7,-6,-5
rows = []
for o in orders:
    rows.append([t.coeff(x,o) for t in totals])
M = sp.Matrix(rows)
print("shape", M.shape, "rank", M.rank())

ns = M.nullspace()
print("nullspace dim", len(ns))

# sanity: also print rank including x^-4 (should be 7, i.e. nullity 8) to match table
rows2 = rows + [[t.coeff(x,-4) for t in totals]]
M2 = sp.Matrix(rows2)
print("rank incl x^-4:", M2.rank())

with open("grado11_base9.pkl","wb") as f:
    pickle.dump(ns, f)

for i, v in enumerate(ns):
    v2 = sp.simplify(v.T)
    print(f"V{i+1} (S1..S15 coeffs):", list(v2))
