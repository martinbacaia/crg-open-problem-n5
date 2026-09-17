import pickle, sympy as sp
import regge_n5_familia2_grado11 as G2

x = G2.x

with open("S115_familia2.pkl","rb") as f:
    totals2 = pickle.load(f)
with open("grado11_base9.pkl","rb") as f:
    ns = pickle.load(f)

for i, v in enumerate(ns):
    combo = sum(v[j]*totals2[j] for j in range(15))
    combo = sp.expand(combo)
    lowest = None
    for o in range(-24, 24):
        c = combo.coeff(x, o)
        if c != 0:
            lowest = o
            break
    print(f"V{i+1}: lowest nonzero order in familia2 = x^{lowest}  ({'CUMPLE <=s^2 (order>=-4)' if lowest is not None and lowest>=-4 else 'VIOLA s^2' if lowest is not None else 'identicamente 0'})")
