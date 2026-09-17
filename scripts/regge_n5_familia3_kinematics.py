"""
TERCERA familia de limites de Regge para n=5 -- patron cualitativo distinto de las 2 anteriores.

Intento previo (documentado en el chat, NO usado): escalar t13 directamente (t13~E^1, intermedio
entre fijo y ~E^2) rompe la paridad de R(x) (da sqrt(R) con potencias SEMI-enteras de x, x^(1/2),
x^(3/2), ...) -- un chequeo de racionalidad (`sqrt(R)` no debe tener `sqrt()` ni potencias
fraccionarias sueltas) lo detecto ANTES de construir nada mas. Diagnostico: el mecanismo de "angulo
tunado" que fija t13 solo admite escalas de M45^2 con exponente PAR en E (por que sqrt(R) tenga
potencias enteras) -- consistente con el "teorema de paridad" ya probado en la parte 21 del proyecto
para invariantes multilineales. Se abandono esa via y se uso la SIGUIENTE, valida:

Familia 3 = misma construccion de la familia 1/2 (cluster duro {1,2,3}, t13 FIJO via el angulo ya
derivado) pero con M45^2 escalando MAS RAPIDO que s12 (en vez de mas lento -- familia 1 tiene M45^2
fijo (~E^0), familia 2 tiene M45^2~E^2 (mismo orden que s12), familia 3 tiene M45^2~E^4 -- CRECE MAS
RAPIDO que el propio canal duro s12). Parametros (t0=1, m0sq0=1) encontrados por busqueda chica para
que sqrt(R) sea EXACTAMENTE racional (mismo criterio que en las familias 1 y 2).
"""
import sympy as sp
from qext import QExt, x, set_R

D = 6
E = 1/x
t0 = sp.Rational(1)
m0sq0 = sp.Rational(1)
m0sq = m0sq0**4 / x**4     # M45^2 ~ E^4 -- crece MAS RAPIDO que s12~E^2 (patron nuevo)
cos_phi0 = sp.Rational(3,5)
sin_phi0 = sp.Rational(4,5)

E3 = E - m0sq/(4*E)
cos_theta = sp.cancel(1 + t0/(2*E*E3))
R = sp.cancel(1 - cos_theta**2)
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
    S = QExt(0,1)
    n3 = [cos_theta, S, 0,0,0]

    p1 = [QExt(-E)] + [QExt(-E*c) for c in n1]
    p2 = [QExt(-E)] + [QExt(-E*(-c)) for c in n1]
    p3 = [QExt(E3)] + [Q(E3)*Q(c) for c in n3]

    P45 = [-(p1[c]+p2[c]+p3[c]) for c in range(D)]

    e_orth = [0,0,1,0,0]
    n45_raw = [cos_phi0*Q(n3[c]) + sin_phi0*Q(e_orth[c]) for c in range(D-1)]
    n45 = n45_raw

    m0 = m0sq0**2/x**2   # = sqrt(m0sq), exacto
    p4_rest = [QExt(m0/2)] + [Q(m0/2)*c for c in n45]
    p5_rest = [QExt(m0/2)] + [-Q(m0/2)*c for c in n45]

    p4 = boost_from_rest(P45, p4_rest, m0)
    p5 = boost_from_rest(P45, p5_rest, m0)

    return {1:p1,2:p2,3:p3,4:p4,5:p5}

if __name__ == "__main__":
    ps = build_kinematics_n5_qext()
    for k in (4,5):
        print(k, [ (c.A, c.B) for c in ps[k] ])
