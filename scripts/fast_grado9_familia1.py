"""
Version RAPIDA (backend lseries.py, sin sympy en el loop caliente) de
regge_n5_simbolico_grado9.py (B1,B2,B3, familia 1). Reconstruye ps_ser/eps_ser con sympy (rapido,
~1-2s, sin cambios) y despues hace TODA la aritmetica (eta, F:F, C, orbit sum sobre 120
permutaciones) con LSeries (Fraction puro). Al final compara, coeficiente a coeficiente, contra
B123.pkl (ya calculado con el pipeline sympy original, de confianza) -- si no coincide EXACTO, esto
no se usa para nada mas.
"""
import time, itertools, pickle, sympy as sp
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

    def Basis1(i,j,l,m,k): return getFF(i,j)*getFF(l,m)*Cc[(k,i,j)]*eta_pp[(i,l)]
    def Basis2(i,j,l,m,k): return getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]
    def Basis3(i,j,l,m,k): return getFF(l,m)*Msand[(i,j,l,k)]*Cc[(k,i,l)]

    t0=time.time()
    B1 = LSeries.zero(); B2 = LSeries.zero(); B3 = LSeries.zero()
    for perm in itertools.permutations(P):
        i,j,l,m,k = perm
        B1 = B1 + Basis1(i,j,l,m,k)
        B2 = B2 + Basis2(i,j,l,m,k)
        B3 = B3 + Basis3(i,j,l,m,k)
    print('orbit sum', time.time()-t0, flush=True)

    with open("B123.pkl","rb") as f:
        B1_ref, B2_ref, B3_ref = pickle.load(f)

    def sympy_dict(expr, lo, hi):
        expr = sp.expand(expr)
        d = {}
        for o in range(lo,hi):
            c = expr.coeff(x,o)
            if c != 0:
                num,den = sp.together(c).as_numer_denom()
                from fractions import Fraction
                d[o] = Fraction(int(num), int(den))
        return d

    ok = True
    for name, fast, ref in [("B1",B1,B1_ref),("B2",B2,B2_ref),("B3",B3,B3_ref)]:
        lo_check, hi_check = -9, 6
        d_fast = {o:v for o,v in fast.c.items() if lo_check<=o<hi_check}
        d_ref = sympy_dict(ref, lo_check, hi_check)
        match = (d_fast == d_ref)
        print(f"{name}: match B123.pkl exact = {match}")
        if not match:
            ok = False
            print("  fast:", d_fast)
            print("  ref :", d_ref)
    print("ALL MATCH" if ok else "MISMATCH FOUND")
