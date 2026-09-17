import time, itertools, sympy as sp
import regge_n5_familia2_kinematics as K

N = 12
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

    print("\nmasslessness check (should be 0 up to O(x^N)):")
    for i in range(1,6):
        e = eta(ps[i], ps[i])
        print(f"  p{i}^2:", e)

    print("\nmomentum conservation check (sum should be 0):")
    tot = [sum(ps[i][c] for i in range(1,6)) for c in range(K.D)]
    tot = [sp.expand(t) for t in tot]
    print(" ", tot)

    print("\npairwise eta(pi,pj), leading order in x:")
    for i,j in itertools.combinations(range(1,6),2):
        e = eta(ps[i],ps[j])
        e = sp.expand(e)
        # lowest order term
        lowest = None
        for o in range(-6,N):
            c = e.coeff(x,o)
            if c != 0:
                lowest = o
                break
        print(f"  ({i},{j}): leading x^{lowest}  coeff={e.coeff(x,lowest) if lowest is not None else 0}")
