import pickle, sympy as sp

with open("S116_familia1.pkl", "rb") as f:
    S1_ls = pickle.load(f)
with open("S116_familia2.pkl", "rb") as f:
    S2_ls = pickle.load(f)


def frac_to_rat(v):
    return sp.Rational(v.numerator, v.denominator)


rows = []
for o in range(-10, -4):
    rows.append([frac_to_rat(t.coeff(o)) for t in S1_ls])
for o in range(-11, -4):
    rows.append([frac_to_rat(t.coeff(o)) for t in S2_ls])

M = sp.Matrix(rows)
print("combined matrix shape", M.shape, "rank", M.rank())
ns = M.nullspace()
print("combined nullspace dim (familia1+familia2, con S16):", len(ns))

for i, v in enumerate(ns):
    combo1 = [sum(v[j] * S1_ls[j].coeff(o) for j in range(16)) for o in range(-15, 6)]
    combo2 = [sum(v[j] * S2_ls[j].coeff(o) for j in range(16)) for o in range(-15, 6)]
    lo1 = next((o for o, c in zip(range(-15, 6), combo1) if c != 0), None)
    lo2 = next((o for o, c in zip(range(-15, 6), combo2) if c != 0), None)
    print(f"U{i+1}: familia1 lowest={lo1}  familia2 lowest={lo2}  coefS16={sp.nsimplify(v[15])}")

with open("grado11_S16_nullspace_combinado.pkl", "wb") as f:
    pickle.dump(ns, f)
