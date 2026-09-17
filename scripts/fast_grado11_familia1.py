"""
Version RAPIDA (lseries.py) de regge_n5_simbolico_grado11.py (S1..S15, familia 1). Ver
fast_grado9_familia1.py para la explicacion del metodo. Compara contra S115.pkl (ya calculado con
el pipeline sympy original) coeficiente a coeficiente antes de confiar en el resultado.
"""
import time, itertools, pickle, sympy as sp
from fractions import Fraction
import regge_n5_simbolico_kinematics as K
from regge_n5_simbolico_eps import make_eps_qext
from lseries import LSeries, eta_vec

N = 16
LO = -6
x = K.x

def build_ls_base():
    Sser = sp.series(sp.sqrt(K.R), x, 0, N).removeO()
    ps = K.build_kinematics_n5_qext()
    eps = make_eps_qext(ps)
    def to_ls(q):
        expr = sp.cancel(q.A + q.B*Sser)
        expr = sp.series(expr, x, 0, N).removeO()
        return LSeries.from_sympy(expr, x, LO, N)
    ps_ls = {i:[to_ls(c) for c in ps[i]] for i in range(1,6)}
    eps_ls = {i:[to_ls(c) for c in eps[i]] for i in range(1,6)}
    return ps_ls, eps_ls

def eta(a,b):
    return eta_vec(a,b,K.D)

if __name__ == "__main__":
    t0=time.time()
    ps_ls, eps_ls = build_ls_base()
    print('base build', time.time()-t0, flush=True)

    P = list(range(1,6))
    t0=time.time()
    eta_pp, eta_pe, eta_ee = {}, {}, {}
    for a in P:
        for b in P:
            eta_pp[(a,b)] = eta(ps_ls[a], ps_ls[b])
            eta_pe[(a,b)] = eta(ps_ls[a], eps_ls[b])
            eta_ee[(a,b)] = eta(eps_ls[a], eps_ls[b])
    print('eta cache', time.time()-t0, flush=True)

    def FdotF(i,j):
        return (eta_pp[(i,j)]*eta_ee[(i,j)] - eta_pe[(j,i)]*eta_pe[(i,j)]) * 2

    def Cabc(k,a,b):
        return eta_pp[(k,a)]*eta_pe[(b,k)] - eta_pp[(k,b)]*eta_pe[(a,k)]

    def Msandwich(i,j,a,b):
        pi_pa = eta_pp[(a,i)]; ei_pj = eta_pe[(j,i)]; pa_ei = eta_pe[(a,i)]; pi_pj = eta_pp[(i,j)]
        Vpj = pi_pa*ei_pj - pa_ei*pi_pj
        ei_ej = eta_ee[(i,j)]; pi_ej = eta_pe[(i,j)]
        Vej = pi_pa*ei_ej - pa_ei*pi_ej
        ej_pb = eta_pe[(b,j)]; pj_pb = eta_pp[(j,b)]
        return Vpj*ej_pb - Vej*pj_pb

    t0=time.time()
    FF, Cc, Msand = {}, {}, {}
    for i in P:
        for j in P:
            if i<j: FF[(i,j)] = FdotF(i,j)
    for k in P:
        for a in P:
            for b in P:
                if len({k,a,b})==3: Cc[(k,a,b)] = Cabc(k,a,b)
    for i in P:
        for j in P:
            if i==j: continue
            for a in P:
                for b in P:
                    if len({i,j,a,b})==4: Msand[(i,j,a,b)] = Msandwich(i,j,a,b)
    print('FF/Cabc/Msandwich cache', time.time()-t0, flush=True)

    def getFF(i,j):
        return FF[(i,j)] if i<j else FF[(j,i)]
    def E(i,j):
        return eta_pp[(i,j)]

    def S1(i,j,l,m,k):  return getFF(i,j)*getFF(l,m)*Cc[(k,i,j)]*E(i,j)*E(i,l)
    def S2(i,j,l,m,k):  return getFF(i,j)*getFF(l,m)*Cc[(k,i,j)]*E(i,l)**2
    def S3(i,j,l,m,k):  return getFF(i,j)*getFF(l,m)*Cc[(k,i,j)]*E(i,l)*E(i,m)
    def S4(i,j,l,m,k):  return getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,j)]*E(i,l)
    def S5(i,j,l,m,k):  return getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,j)]*E(k,l)
    def S6(i,j,l,m,k):  return getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]*E(l,m)
    def S7(i,j,l,m,k):  return getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]*E(i,j)
    def S8(i,j,l,m,k):  return getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]*E(i,l)
    def S9(i,j,l,m,k):  return getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]*E(k,l)
    def S10(i,j,l,m,k): return getFF(l,m)*Msand[(i,j,l,k)]*Cc[(k,i,j)]*E(i,l)
    def S11(i,j,l,m,k): return getFF(l,m)*Msand[(i,j,l,k)]*Cc[(k,i,l)]*E(l,m)
    def S12(i,j,l,m,k): return getFF(l,m)*Msand[(i,j,l,k)]*Cc[(k,i,l)]*E(i,l)
    def S13(i,j,l,m,k): return Msand[(i,j,l,m)]*Msand[(l,m,i,k)]*Cc[(k,i,j)]
    def S14(i,j,l,m,k): return Msand[(i,j,l,m)]*Msand[(l,m,i,k)]*Cc[(k,l,m)]
    def S15(i,j,l,m,k): return Msand[(i,j,l,m)]*Msand[(l,m,i,k)]*Cc[(k,i,l)]

    BASIS = [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15]

    t0=time.time()
    totals = [LSeries.zero() for _ in BASIS]
    for perm in itertools.permutations(P):
        i,j,l,m,k = perm
        for idx, fn in enumerate(BASIS):
            totals[idx] = totals[idx] + fn(i,j,l,m,k)
    print('orbit sum (15 basis)', time.time()-t0, flush=True)

    with open("S115.pkl","rb") as f:
        totals_ref = pickle.load(f)

    def sympy_dict(expr, lo, hi):
        expr = sp.expand(expr)
        d = {}
        for o in range(lo,hi):
            c = expr.coeff(x,o)
            if c != 0:
                num,den = sp.together(c).as_numer_denom()
                d[o] = Fraction(int(num), int(den))
        return d

    ok = True
    lo_check, hi_check = -10, 6
    for idx in range(15):
        d_fast = {o:v for o,v in totals[idx].c.items() if lo_check<=o<hi_check}
        d_ref = sympy_dict(totals_ref[idx], lo_check, hi_check)
        match = (d_fast == d_ref)
        print(f"S{idx+1}: match S115.pkl exact = {match}")
        if not match:
            ok = False
            print("  fast:", d_fast)
            print("  ref :", d_ref)
    print("ALL MATCH" if ok else "MISMATCH FOUND")

    with open("S115_fast.pkl","wb") as f:
        pickle.dump(totals, f)
