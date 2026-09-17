import pickle, sympy as sp

with open("grado11_S16_nullspace_combinado.pkl", "rb") as f:
    ns = pickle.load(f)
with open("S116_familia3_N30.pkl", "rb") as f:
    S3_ls = pickle.load(f)

for i, v in enumerate(ns):
    combo = [sum(v[j] * S3_ls[j].coeff(o) for j in range(16)) for o in range(-24, 6)]
    lo = next((o for o, c in zip(range(-24, 6), combo) if c != 0), None)
    print(f"U{i+1} en familia 3 (N=30): lowest order = {lo}")
