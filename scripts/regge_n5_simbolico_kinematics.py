"""
Fase 2, paso 30: reconstruccion SIMBOLICA EXACTA (sympy) del limite de Regge de n=5
(scripts/regge_limit_n5.py, parte 22/23) via x=1/E, x->0.

Idea clave para que sea tratable: en vez de invocar sqrt()/series() sobre expresiones anidadas (muy
lento, ver 30-derivacion-simbolica-cascada.md), se explota que sin_theta=S es la UNICA fuente de
irracionalidad algebraica no trivial en toda la construccion (aparte de sqrt(m0sq), evitada eligiendo
m0sq perfecto cuadrado), y aparece siempre LINEALMENTE (nunca al cuadrado sin reducir) en cada
componente de cada p_i, eps_i. Por eso toda la cinematica vive exactamente en la extension cuadratica
Q(x)(S), S^2=R(x) (ver qext.py) -- aritmetica exacta, sin sqrt() literal, hasta el final, donde se
sustituye S por su serie de potencias real (sp.series(sqrt(R),x,0,N), rapido porque R es una funcion
racional simple, sin nada anidado).

m0sq, t0 se toman como sp.Rational elegidos PERFECTOS CUADRADOS/simples para mantener todo en Q
(ningun sqrt(entero) suelto en el resultado final, mas facil de comparar entre configuraciones).
phi0 (angulo fijo del "leftover" decay, ver regge_limit_n5.py) se toma como un punto racional del
circulo unitario (cos_phi0,sin_phi0) en vez de un angulo generico, por la misma razon.

Para la verificacion cruzada de independencia de parametros (30-derivacion-simbolica-cascada.md):
cambiar t0, m0sq, cos_phi0/sin_phi0 (y las SEEDS de regge_n5_simbolico_eps.py) a otro punto racional
generico y re-correr regge_n5_simbolico_grado9.py / _grado11.py -- deben reproducirse los MISMOS
cocientes exactos (A1/A2=2, A(S8)/A(S1)=-1/2, A(S14)/A(S15)=1) y el mismo orden de "piso" de la
cascada de cancelacion, si estos son identidades cinematicas genuinas y no accidentes numericos.
"""
import sympy as sp
from qext import QExt, x, set_R

D = 6
E = 1/x
t0 = sp.Rational(-1)
m0sq = sp.Rational(4)
cos_phi0 = sp.Rational(3,5)
sin_phi0 = sp.Rational(4,5)

E3 = E - m0sq/(4*E)
cos_theta = sp.cancel(1 + t0/(2*E*E3))
R = sp.cancel(1 - cos_theta**2)   # S^2 = R(x), S = sin_theta
set_R(R)

def Q(v):
    return v if isinstance(v, QExt) else QExt(v)

def eta(a,b):
    tot = -Q(a[0])*Q(b[0])
    for k in range(1,D):
        tot = tot + Q(a[k])*Q(b[k])
    return tot

def boost_from_rest(P, q_rest, m0_val):
    M = m0_val
    P0 = P[0]
    Pvec = P[1:]
    beta_vec = [Q(c)/Q(P0) for c in Pvec]
    beta2 = QExt(0)
    for b in beta_vec:
        beta2 = beta2 + b*b
    gamma = Q(P0)/Q(M)
    q0, qvec = q_rest[0], q_rest[1:]
    bq = QExt(0)
    for k in range(D-1):
        bq = bq + beta_vec[k]*Q(qvec[k])
    q0p = gamma*(Q(q0)+bq)
    qvecp = []
    for k in range(D-1):
        term = Q(qvec[k]) + (gamma-QExt(1))*(bq/beta2)*beta_vec[k] + gamma*Q(q0)*beta_vec[k]
        qvecp.append(term)
    return [q0p]+qvecp

def build_kinematics_n5_qext():
    n1 = [sp.Integer(0)]*(D-1); n1[0]=sp.Integer(1)
    S = QExt(0,1)  # sin_theta
    n3 = [cos_theta, S, 0,0,0]

    p1 = [QExt(-E)] + [QExt(-E*c) for c in n1]
    p2 = [QExt(-E)] + [QExt(-E*(-c)) for c in n1]
    p3 = [QExt(E3)] + [Q(E3)*Q(c) for c in n3]

    P45 = [-(p1[c]+p2[c]+p3[c]) for c in range(D)]

    e_orth = [0,0,1,0,0]  # index2 = x3 direction
    n45_raw = [cos_phi0*Q(n3[c]) + sin_phi0*Q(e_orth[c]) for c in range(D-1)]
    # norm should be exactly 1 (cos_phi0^2+sin_phi0^2=1, cos_theta^2+S^2=1) -- verify later
    n45 = n45_raw

    m0 = sp.sqrt(m0sq)  # =2 rational
    p4_rest = [QExt(m0/2)] + [Q(m0/2)*c for c in n45]
    p5_rest = [QExt(m0/2)] + [-Q(m0/2)*c for c in n45]

    p4 = boost_from_rest(P45, p4_rest, m0)
    p5 = boost_from_rest(P45, p5_rest, m0)

    return {1:p1,2:p2,3:p3,4:p4,5:p5}

if __name__ == "__main__":
    ps = build_kinematics_n5_qext()
    for k in (4,5):
        print(k, [ (c.A, c.B) for c in ps[k] ])
