"""
Version RAPIDA (lseries.py) de S1..S15 para la familia 3 (M45^2~E^4). A diferencia de familias 1/2
(donde los pares "crecientes" son todos ~E^2), aca hay pares que crecen mas rapido (~E^4) -- el
grado (en x) de los S1..S15 puede ser mucho mas negativo que en las familias anteriores, asi que
hace falta un N (orden de truncacion) bastante mas generoso. Se corre con 2 valores de N distintos
y se exige que el resultado en la ventana de interes (ordenes -4 a -30, digamos) sea IDENTICO entre
ambos -- si difiere, N es insuficiente y hay que subirlo (chequeo de estabilidad, no aceptar N
"a ojo").
"""
import time, itertools, sys, sympy as sp
import regge_n5_familia3_kinematics as K
from regge_n5_familia3_eps import make_eps_qext
from lseries import LSeries, eta_vec

x = K.x

def build_ls_base(N, LO=-8):
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

def compute_S115(N, LO=-8):
    ps_ls, eps_ls = build_ls_base(N, LO)
    P = list(range(1,6))
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

if __name__ == "__main__":
    for N in (30, 45):
        t0=time.time()
        totals = compute_S115(N, LO=-8)
        print(f"N={N}: orbit sum time", time.time()-t0, flush=True)
        for idx,t in enumerate(totals):
            lo = t.lowest_order()
            print(f"  S{idx+1}: lowest order {lo}")
        import pickle
        with open(f"S115_familia3_N{N}.pkl","wb") as f:
            pickle.dump(totals, f)
