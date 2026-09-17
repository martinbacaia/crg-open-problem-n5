import time, itertools, sympy as sp
import regge_n5_familia2_kinematics as K
from regge_n5_familia2_eps import make_eps_qext

N = 20
x = K.x

def build_series(N=N):
    Sser = sp.series(sp.sqrt(K.R), x, 0, N).removeO()
    ps = K.build_kinematics_n5_qext()
    eps = make_eps_qext(ps)
    def to_series(q):
        expr = sp.cancel(q.A + q.B*Sser)
        return sp.series(expr, x, 0, N).removeO()
    ps_ser = {i:[to_series(c) for c in ps[i]] for i in range(1,6)}
    eps_ser = {i:[to_series(c) for c in eps[i]] for i in range(1,6)}
    return ps_ser, eps_ser

def trunc(expr, N=N):
    return sp.series(sp.expand(expr), x, 0, N).removeO()

def eta(a,b):
    tot = -a[0]*b[0] + sum(a[k]*b[k] for k in range(1,K.D))
    return trunc(tot)

if __name__ == "__main__":
    t0=time.time()
    ps_ser, eps_ser = build_series()
    print('base series build', time.time()-t0, flush=True)

    P = list(range(1,6))
    t0=time.time()
    eta_pp = {}
    eta_pe = {}
    eta_ee = {}
    for a in P:
        for b in P:
            eta_pp[(a,b)] = eta(ps_ser[a], ps_ser[b])
            eta_pe[(a,b)] = eta(ps_ser[a], eps_ser[b])
            eta_ee[(a,b)] = eta(eps_ser[a], eps_ser[b])
    print('eta cache time', time.time()-t0, flush=True)

    def FdotF(i,j):
        return trunc(2*(eta_pp[(i,j)]*eta_ee[(i,j)] - eta_pe[(j,i)]*eta_pe[(i,j)]))

    def Cabc(k,a,b):
        return trunc(eta_pp[(k,a)]*eta_pe[(b,k)] - eta_pp[(k,b)]*eta_pe[(a,k)])

    def Msandwich(i,j,a,b):
        pi_pa = eta_pp[(a,i)]
        ei_pj = eta_pe[(j,i)]
        pa_ei = eta_pe[(a,i)]
        pi_pj = eta_pp[(i,j)]
        Vpj = trunc(pi_pa*ei_pj - pa_ei*pi_pj)
        ei_ej = eta_ee[(i,j)]
        pi_ej = eta_pe[(i,j)]
        Vej = trunc(pi_pa*ei_ej - pa_ei*pi_ej)
        ej_pb = eta_pe[(b,j)]
        pj_pb = eta_pp[(j,b)]
        return trunc(Vpj*ej_pb - Vej*pj_pb)

    t0=time.time()
    FF = {}
    for i in P:
        for j in P:
            if i<j:
                FF[(i,j)] = FdotF(i,j)
    Cc = {}
    for k in P:
        for a in P:
            for b in P:
                if len({k,a,b})==3:
                    Cc[(k,a,b)] = Cabc(k,a,b)
    Msand = {}
    for i in P:
        for j in P:
            if i==j: continue
            for a in P:
                for b in P:
                    if len({i,j,a,b})==4:
                        Msand[(i,j,a,b)] = Msandwich(i,j,a,b)
    print('FF/Cabc/Msandwich cache time', time.time()-t0, flush=True)

    def getFF(i,j):
        return FF[(i,j)] if i<j else FF[(j,i)]
    def E(i,j):
        return eta_pp[(i,j)]

    def S1(i,j,l,m,k):  return trunc(getFF(i,j)*getFF(l,m)*Cc[(k,i,j)]*E(i,j)*E(i,l))
    def S2(i,j,l,m,k):  return trunc(getFF(i,j)*getFF(l,m)*Cc[(k,i,j)]*E(i,l)**2)
    def S3(i,j,l,m,k):  return trunc(getFF(i,j)*getFF(l,m)*Cc[(k,i,j)]*E(i,l)*E(i,m))
    def S4(i,j,l,m,k):  return trunc(getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,j)]*E(i,l))
    def S5(i,j,l,m,k):  return trunc(getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,j)]*E(k,l))
    def S6(i,j,l,m,k):  return trunc(getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]*E(l,m))
    def S7(i,j,l,m,k):  return trunc(getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]*E(i,j))
    def S8(i,j,l,m,k):  return trunc(getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]*E(i,l))
    def S9(i,j,l,m,k):  return trunc(getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)]*E(k,l))
    def S10(i,j,l,m,k): return trunc(getFF(l,m)*Msand[(i,j,l,k)]*Cc[(k,i,j)]*E(i,l))
    def S11(i,j,l,m,k): return trunc(getFF(l,m)*Msand[(i,j,l,k)]*Cc[(k,i,l)]*E(l,m))
    def S12(i,j,l,m,k): return trunc(getFF(l,m)*Msand[(i,j,l,k)]*Cc[(k,i,l)]*E(i,l))
    def S13(i,j,l,m,k): return trunc(Msand[(i,j,l,m)]*Msand[(l,m,i,k)]*Cc[(k,i,j)])
    def S14(i,j,l,m,k): return trunc(Msand[(i,j,l,m)]*Msand[(l,m,i,k)]*Cc[(k,l,m)])
    def S15(i,j,l,m,k): return trunc(Msand[(i,j,l,m)]*Msand[(l,m,i,k)]*Cc[(k,i,l)])

    BASIS = [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15]

    t0=time.time()
    totals = [sp.Integer(0) for _ in BASIS]
    for perm in itertools.permutations(P):
        i,j,l,m,k = perm
        for idx, fn in enumerate(BASIS):
            totals[idx] = totals[idx] + fn(i,j,l,m,k)
    totals = [sp.expand(t) for t in totals]
    print('orbit sum time (15 basis)', time.time()-t0, flush=True)

    import pickle
    with open("S115_familia2.pkl","wb") as f:
        pickle.dump(totals, f)

    for idx,t in enumerate(totals):
        print(f"S{idx+1}: lowest nonzero order:", end=" ")
        lo = None
        for o in range(-24,24):
            c = t.coeff(x,o)
            if c != 0:
                lo = o
                break
        print(lo)
