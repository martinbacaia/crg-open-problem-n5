import time, pickle, sympy as sp
import regge_n5_familia3_kinematics as K
from regge_n5_familia3_eps import make_eps_qext
from lseries import LSeries
from basis_S116 import compute_S116

x = K.x


def build_ls_base(N, LO=-8):
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
    for N in (30, 45):
        t0 = time.time()
        ps_ls, eps_ls = build_ls_base(N, LO=-8)
        totals = compute_S116(ps_ls, eps_ls, K.D)
        print(f"N={N}: orbit sum time", time.time() - t0, flush=True)
        for idx, t in enumerate(totals):
            print(f"  S{idx+1}: lowest order {t.lowest_order()}")
        with open(f"S116_familia3_N{N}.pkl", "wb") as f:
            pickle.dump(totals, f)
