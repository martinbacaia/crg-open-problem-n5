import time, itertools, sympy as sp
import regge_n5_simbolico_kinematics as K
from regge_n5_simbolico_eps import make_eps_qext

N = 16
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

    # precompute caches
    P = list(range(1,6))
    t0=time.time()
    eta_pp = {}
    eta_pe = {}   # eta(p_a, eps_b)  (order matters)
    eta_ee = {}
    for a in P:
        for b in P:
            eta_pp[(a,b)] = eta(ps_ser[a], ps_ser[b])
            eta_pe[(a,b)] = eta(ps_ser[a], eps_ser[b])
            eta_ee[(a,b)] = eta(eps_ser[a], eps_ser[b])
    print('eta cache time', time.time()-t0, flush=True)

    def FdotF(i,j):
        return trunc(2*(eta_pp[(i,j)]*eta_ee[(i,j)] - eta_pe[(j,i)]*eta_pe[(i,j)]))
        # note: eta(pi,ej)=eta_pe[(i,j)], eta(ei,pj)=eta_pe[(j,i)]

    def Cabc(k,a,b):
        return trunc(eta_pp[(k,a)]*eta_pe[(b,k)] - eta_pp[(k,b)]*eta_pe[(a,k)])
        # eta(pk,pa)=eta_pp[(k,a)], eta(ek,pb)=eta_pe[(b,k)]... wait eta_pe defined as eta(p_a,eps_b);
        # eta(ek,pb) = eta(pb, ek) by symmetry = eta_pe[(b,k)]. OK.

    def Msandwich(i,j,a,b):
        pi_pa = eta_pp[(a,i)]
        ei_pj = eta_pe[(j,i)]   # eta(pj,ei) = eta_pe[(j,i)]
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

    def Basis1(i,j,l,m,k):
        return trunc(getFF(i,j)*getFF(l,m)*Cc[(k,i,j)]*eta_pp[(i,l)])
    def Basis2(i,j,l,m,k):
        return trunc(getFF(l,m)*Msand[(i,j,l,m)]*Cc[(k,i,l)])
    def Basis3(i,j,l,m,k):
        return trunc(getFF(l,m)*Msand[(i,j,l,k)]*Cc[(k,i,l)])

    t0=time.time()
    B1 = sp.Integer(0); B2 = sp.Integer(0); B3 = sp.Integer(0)
    for perm in itertools.permutations(P):
        i,j,l,m,k = perm
        B1 = B1 + Basis1(i,j,l,m,k)
        B2 = B2 + Basis2(i,j,l,m,k)
        B3 = B3 + Basis3(i,j,l,m,k)
    B1 = sp.expand(B1); B2 = sp.expand(B2); B3 = sp.expand(B3)
    print('orbit sum time', time.time()-t0, flush=True)

    import pickle
    with open("B123.pkl","wb") as f:
        pickle.dump((B1,B2,B3), f)
    print("B1 leading terms around x^-8..x^-4:")
    for pw in range(-8,-3):
        print(pw, B1.coeff(x,pw))
    print("B2:")
    for pw in range(-8,-3):
        print(pw, B2.coeff(x,pw))
    print("B3:")
    for pw in range(-8,-3):
        print(pw, B3.coeff(x,pw))
