import sympy as sp
from qext import QExt, x
from regge_n5_familia3_kinematics import D, build_kinematics_n5_qext, eta, Q

growing_partner = {1:2,2:1,3:2,4:1,5:1}

SEEDS = {
    1:[2,-3,1,5,-2,4],
    2:[-1,4,-5,2,3,-4],
    3:[3,1,2,-4,5,-1],
    4:[-2,5,3,1,-4,2],
    5:[4,-2,-1,3,1,-5],
}

def make_eps_qext(ps):
    eps = {}
    for i in range(1,6):
        e = [QExt(sp.Rational(v)) for v in SEEDS[i]]
        other = growing_partner[i]
        coef = eta(e, ps[i]) / eta(ps[other], ps[i])
        e2 = [e[c] - coef*ps[other][c] for c in range(D)]
        eps[i] = e2
    return eps

if __name__ == "__main__":
    ps = build_kinematics_n5_qext()
    eps = make_eps_qext(ps)
    print("ok")
