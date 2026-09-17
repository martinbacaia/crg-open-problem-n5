"""
Fase 2, paso 32: SEGUNDA familia de limites de Regge para n=5 -- limite JERARQUICO de 2 escalas.

Reutiliza EXACTAMENTE la construccion de la familia 1 (regge_n5_simbolico_kinematics.py, partes
22/23/30): mismo cluster duro {1,2,3} (s12 creciente, t13 fijo via el angulo theta(E) tunado) +
mismo cluster de retroceso {4,5} (boost de un decaimiento de 2 cuerpos fijo en el rest-frame de
P45). La UNICA diferencia: en vez de mantener m0sq:=M45^2 FIJO (orden x^0), se lo hace crecer con
E a una tasa distinta (mas lenta) que s12~E^2 -- un limite de Regge "iterado"/multi-escala genuino,
NO alcanzable por reetiquetado de la familia 1 (ver 31-base-explicita-grado11-y-problema-segunda-
familia.md para por que el reetiquetado simple no sirve).

Se elige el exponente MAS SIMPLE que preserva potencias enteras de x=1/E en toda la aritmetica
QExt (evitando raices fraccionarias de x): m0sq = (m0sq0)^2 / x^2, es decir m0 = m0sq0/x ~ E ->
M45^2 ~ E^2, EL MISMO orden que s12 (pero segunda escala genuinamente independiente: t13 sigue
EXACTO en t0 por construccion -- la derivacion de cos_theta,E3 que fija t13=t0 y M45^2=m0sq es
valida para CUALQUIER m0sq, incluso dependiente de x, ver comentario en regge_limit_n5.py).

La formula exacta de cos_theta,E3 (derivada analiticamente para CUALQUIER m0sq, ver
regge_limit_n5.py) se mantiene sin cambios -- solo se sustituye el VALOR de m0sq.
"""
import sympy as sp
from qext import QExt, x, set_R

D = 6
E = 1/x
# t0, m0sq0 elegidos (busqueda exhaustiva chica) para que R(x) de una serie sqrt(R) RACIONAL --
# igual que en la familia 1 ("elegidos perfectos cuadrados/simples para mantener todo en Q"). Un
# primer intento con t0=-1,m0sq0=1 daba R con sqrt(3) sin cancelar hasta la suma final sobre S5
# (bug real encontrado por el cross-check con lseries.py, ver 32-...md) -- corregido aca.
t0 = sp.Rational(-3)
m0sq0 = sp.Rational(1)          # escala del segundo parametro jerarquico
m0sq = m0sq0**2 / x**2          # M45^2 ~ E^2 (segunda escala, independiente de t0, misma
                                 # potencia que s12 pero NO ligada a ella por construccion)
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

    m0 = m0sq0/x   # = sqrt(m0sq), exacto (m0sq0>0, x>0)
    p4_rest = [QExt(m0/2)] + [Q(m0/2)*c for c in n45]
    p5_rest = [QExt(m0/2)] + [-Q(m0/2)*c for c in n45]

    p4 = boost_from_rest(P45, p4_rest, m0)
    p5 = boost_from_rest(P45, p5_rest, m0)

    return {1:p1,2:p2,3:p3,4:p4,5:p5}

if __name__ == "__main__":
    ps = build_kinematics_n5_qext()
    for k in (4,5):
        print(k, [ (c.A, c.B) for c in ps[k] ])
