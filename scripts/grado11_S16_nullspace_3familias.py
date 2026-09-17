"""
Analisis DEFINITIVO: nucleo combinado de las 3 familias de Regge a la vez (no solo evaluar U1,U2,U3
por separado), con S16 incluido -- la pregunta pendiente marcada como "mas urgente" al cierre de la
parte 36. Si este nucleo combinado es {0}, confirma que NINGUNA combinacion de S1..S16 satisface
<=s^2 en las 3 familias simultaneamente (no solo que las 3 direcciones ya encontradas fallan).
"""
import pickle, sympy as sp

with open("S116_familia1.pkl", "rb") as f:
    S1_ls = pickle.load(f)
with open("S116_familia2.pkl", "rb") as f:
    S2_ls = pickle.load(f)
with open("S116_familia3_N45.pkl", "rb") as f:
    S3_ls = pickle.load(f)


def frac_to_rat(v):
    return sp.Rational(v.numerator, v.denominator)


rows = []
for o in range(-10, -4):
    rows.append([frac_to_rat(t.coeff(o)) for t in S1_ls])
for o in range(-11, -4):
    rows.append([frac_to_rat(t.coeff(o)) for t in S2_ls])
for o in range(-22, -4):
    rows.append([frac_to_rat(t.coeff(o)) for t in S3_ls])

M = sp.Matrix(rows)
print("matriz combinada (3 familias, con S16) shape", M.shape, "rank", M.rank())
ns = M.nullspace()
print("nucleo combinado (3 familias, con S16) dimension:", len(ns))
for i, v in enumerate(ns):
    print(f"V{i+1}:", list(v.T))

with open("grado11_S16_nullspace_3familias.pkl", "wb") as f:
    pickle.dump(ns, f)
