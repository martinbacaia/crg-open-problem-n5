import time, itertools, sympy as sp
import regge_n5_familia3_kinematics as K

N = 24
x = K.x

def build_series(N=N):
    Sser = sp.series(sp.sqrt(K.R), x, 0, N).removeO()
    ps = K.build_kinematics_n5_qext()
    def to_series(q):
        expr = sp.cancel(q.A + q.B*Sser)
        return sp.series(expr, x, 0, N).removeO()
    ps_ser = {i:[to_series(c) for c in ps[i]] for i in range(1,6)}
    return ps_ser

def trunc(expr, N=N):
    return sp.series(sp.expand(expr), x, 0, N).removeO()

def eta(a,b):
    tot = -a[0]*b[0] + sum(a[k]*b[k] for k in range(1,K.D))
    return trunc(tot)

if __name__ == "__main__":
    t0=time.time()
    ps = build_series()
    print("build time", time.time()-t0)

    print("\nmasslessness (region de interes, orden <=10, debe ser 0):")
    for i in range(1,6):
        e = eta(ps[i], ps[i])
        e = sp.expand(e)
        nz = [(o,e.coeff(x,o)) for o in range(-12,11) if e.coeff(x,o)!=0]
        print(f"  p{i}^2 nonzero terms (orden<=10):", nz)

    print("\nmomentum conservation:")
    tot = [sum(ps[i][c] for i in range(1,6)) for c in range(K.D)]
    print(" ", [sp.expand(t) for t in tot])

    print("\npairwise eta(pi,pj), leading order:")
    for i,j in itertools.combinations(range(1,6),2):
        e = eta(ps[i],ps[j])
        e = sp.expand(e)
        lowest = None
        for o in range(-8, N):
            c = e.coeff(x,o)
            if c != 0:
                lowest = o
                break
        print(f"  ({i},{j}): leading x^{lowest}  coeff={e.coeff(x,lowest) if lowest is not None else 0}")
