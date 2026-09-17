"""
Familia 1, S1..S16 (agrega S16 = topologia T nueva, ver basis_S116.py). Cross-check obligatorio:
S1..S15 deben coincidir EXACTO con S115_fast.pkl (ya validado contra el pipeline sympy original)
antes de confiar en S16.
"""
import time, pickle, sympy as sp
import regge_n5_simbolico_kinematics as K
from regge_n5_simbolico_eps import make_eps_qext
from lseries import LSeries
from basis_S116 import compute_S116

N = 16
LO = -6
x = K.x


def build_ls_base():
    Sser = sp.series(sp.sqrt(K.R), x, 0, N).removeO()
    ps = K.build_kinematics_n5_qext()
    eps = make_eps_qext(ps)
    def to_ls(q):
        expr = sp.cancel(q.A + q.B * Sser)
        expr = sp.series(expr, x, 0, N).removeO()
        return LSeries.from_sympy(expr, x, LO, N)
    ps_ls = {i: [to_ls(c) for c in ps[i]] for i in range(1, 6)}
    eps_ls = {i: [to_ls(c) for c in eps[i]] for i in range(1, 6)}
    return ps_ls, eps_ls


if __name__ == "__main__":
    t0 = time.time()
    ps_ls, eps_ls = build_ls_base()
    print("base build", time.time() - t0, flush=True)

    t0 = time.time()
    totals = compute_S116(ps_ls, eps_ls, K.D)
    print("orbit sum (16 basis)", time.time() - t0, flush=True)

    with open("S115_fast.pkl", "rb") as f:
        ref15 = pickle.load(f)

    ok = True
    for idx in range(15):
        match = (totals[idx].c == ref15[idx].c)
        if not match:
            ok = False
            print(f"S{idx+1}: MISMATCH vs S115_fast.pkl")
    print("S1..S15 coinciden con S115_fast.pkl:", ok)

    lo16 = totals[15].lowest_order()
    print(f"S16: lowest order = {lo16}, is_zero = {totals[15].is_zero()}")

    with open("S116_familia1.pkl", "wb") as f:
        pickle.dump(totals, f)
