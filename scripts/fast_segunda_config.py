"""
Chequeo de robustez pendiente (marcado como el mas urgente en 32-...md): repetir TODO el analisis
(grado 11, familias 1 y 2, subespacio universal) con una SEGUNDA configuracion de parametros
distinta en cada familia. Ahora es barato gracias al backend rapido (lseries.py, ya cross-validado
exacto contra sympy en la configuracion de referencia).

Familia 1, config 2: t0=-4, m0sq=9, cos_phi0=5/13, sin_phi0=12/13 (la misma segunda config ya usada
en la parte 30 para el piso s^3 de grado 9 -- reutilizada aca para grado 11 tambien).
Familia 2, config 2: t0=-12, m0sq0=1, cos_phi0=5/13, sin_phi0=12/13 (mismo m0sq0=1 que la config de
referencia -- garantiza sqrt(R) racional, ver busqueda en el chat -- pero t0 y phi0 distintos).
"""
import time, itertools, pickle, sympy as sp
from fractions import Fraction
from lseries import LSeries, eta_vec

def run_family1_cfg2():
    import regge_n5_simbolico_kinematics as K
    import regge_n5_simbolico_eps as EPS
    x = K.x
    # monkeypatch params
    K.t0 = sp.Rational(-4)
    K.m0sq = sp.Rational(9)
    K.cos_phi0 = sp.Rational(5,13)
    K.sin_phi0 = sp.Rational(12,13)
    K.E3 = K.E - K.m0sq/(4*K.E)
    K.cos_theta = sp.cancel(1 + K.t0/(2*K.E*K.E3))
    K.R = sp.cancel(1 - K.cos_theta**2)
    K.set_R(K.R)
    Sser_check = sp.series(sp.sqrt(K.R), x, 0, 8)
    assert not (Sser_check.has(sp.sqrt) or Sser_check.has(sp.I)), "familia1 cfg2: sqrt(R) no racional"

    N = 16; LO = -6
    Sser = sp.series(sp.sqrt(K.R), x, 0, N).removeO()
    ps = K.build_kinematics_n5_qext()
    eps = EPS.make_eps_qext(ps)
    def to_ls(q):
        expr = sp.cancel(q.A + q.B*Sser)
        expr = sp.series(expr, x, 0, N).removeO()
        return LSeries.from_sympy(expr, x, LO, N)
    ps_ls = {i:[to_ls(c) for c in ps[i]] for i in range(1,6)}
    eps_ls = {i:[to_ls(c) for c in eps[i]] for i in range(1,6)}
    return compute_S115(ps_ls, eps_ls, K.D)

def run_family2_cfg2():
    import regge_n5_familia2_kinematics as K
    import regge_n5_familia2_eps as EPS
    x = K.x
    K.t0 = sp.Rational(-12)
    K.m0sq0 = sp.Rational(1)
    K.m0sq = K.m0sq0**2 / x**2
    K.cos_phi0 = sp.Rational(5,13)
    K.sin_phi0 = sp.Rational(12,13)
    K.E3 = K.E - K.m0sq/(4*K.E)
    K.cos_theta = sp.cancel(1 + K.t0/(2*K.E*K.E3))
    K.R = sp.cancel(1 - K.cos_theta**2)
    K.set_R(K.R)
    Sser_check = sp.series(sp.sqrt(K.R), x, 0, 8)
    assert not (Sser_check.has(sp.sqrt) or Sser_check.has(sp.I)), "familia2 cfg2: sqrt(R) no racional"

    N = 20; LO = -6
    Sser = sp.series(sp.sqrt(K.R), x, 0, N).removeO()
    ps = K.build_kinematics_n5_qext()
    eps = EPS.make_eps_qext(ps)
    def to_ls(q):
        expr = sp.cancel(q.A + q.B*Sser)
        expr = sp.series(expr, x, 0, N).removeO()
        return LSeries.from_sympy(expr, x, LO, N)
    ps_ls = {i:[to_ls(c) for c in ps[i]] for i in range(1,6)}
    eps_ls = {i:[to_ls(c) for c in eps[i]] for i in range(1,6)}
    return compute_S115(ps_ls, eps_ls, K.D)

def compute_S115(ps_ls, eps_ls, D):
    P = list(range(1,6))
    def eta(a,b): return eta_vec(a,b,D)
    eta_pp, eta_pe, eta_ee = {}, {}, {}
    for a in P:
        for b in P:
            eta_pp[(a,b)] = eta(ps_ls[a], ps_ls[b])
            eta_pe[(a,b)] = eta(ps_ls[a], eps_ls[b])
            eta_ee[(a,b)] = eta(eps_ls[a], eps_ls[b])

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

    def getFF(i,j): return FF[(i,j)] if i<j else FF[(j,i)]
    def E(i,j): return eta_pp[(i,j)]

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
    totals = [LSeries.zero() for _ in BASIS]
    for perm in itertools.permutations(P):
        i,j,l,m,k = perm
        for idx, fn in enumerate(BASIS):
            totals[idx] = totals[idx] + fn(i,j,l,m,k)
    return totals

def lowest_order(ls):
    if ls.is_zero(): return None
    return ls.lowest_order()

if __name__ == "__main__":
    t0=time.time()
    S1_ls = run_family1_cfg2()
    print("familia1 cfg2 orbit sum done", time.time()-t0, flush=True)
    for i,t in enumerate(S1_ls):
        print(f"  S{i+1} lowest order: {lowest_order(t)}")

    t0=time.time()
    S2_ls = run_family2_cfg2()
    print("familia2 cfg2 orbit sum done", time.time()-t0, flush=True)
    for i,t in enumerate(S2_ls):
        print(f"  S{i+1} lowest order: {lowest_order(t)}")

    # cross-check A(S8)=-1/2 A(S1) y A(S14)=A(S15) (relaciones exactas ya encontradas en la parte 30)
    o1 = lowest_order(S1_ls[0])
    print("familia1 cfg2: A(S8)/A(S1) at leading order", S1_ls[7].coeff(o1)/S1_ls[0].coeff(o1))
    o14 = lowest_order(S1_ls[13])
    print("familia1 cfg2: A(S14)/A(S15) at leading order", S1_ls[13].coeff(o14)/S1_ls[14].coeff(o14))

    # matriz combinada (misma logica que antes) para el subespacio universal
    rows = []
    for o in range(-10,-4):
        rows.append([t.coeff(o) for t in S1_ls])
    for o in range(-11,-4):
        rows.append([t.coeff(o) for t in S2_ls])
    M = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in row] for row in rows])
    print("combined matrix shape", M.shape, "rank", M.rank())
    ns = M.nullspace()
    print("combined nullspace dim (cfg2)", len(ns))

    for i,v in enumerate(ns):
        combo1 = [sum(v[j]*S1_ls[j].coeff(o) for j in range(15)) for o in range(-15,6)]
        combo2 = [sum(v[j]*S2_ls[j].coeff(o) for j in range(15)) for o in range(-15,6)]
        lo1 = next((o for o,c in zip(range(-15,6),combo1) if c!=0), None)
        lo2 = next((o for o,c in zip(range(-15,6),combo2) if c!=0), None)
        print(f"U{i+1} (cfg2): familia1 lowest={lo1}  familia2 lowest={lo2}")

    with open("S115_cfg2_check.pkl","wb") as f:
        pickle.dump((S1_ls, S2_ls, ns), f)
