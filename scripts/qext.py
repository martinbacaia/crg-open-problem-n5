"""
Aritmetica exacta en la extension cuadratica Q(x)(S), S^2 = R(x) (funcion racional conocida).
cancel() se aplica tras cada operacion para evitar el crecimiento descontrolado de fracciones
racionales anidadas (imprescindible para que series() al final sea rapido).
"""
import sympy as sp

x = sp.Symbol('x', positive=True)

class QExt:
    __slots__ = ('A','B')
    def __init__(self, A, B=0):
        self.A = sp.sympify(A)
        self.B = sp.sympify(B)
    def __add__(self, o):
        o = o if isinstance(o, QExt) else QExt(o)
        return QExt(sp.cancel(self.A+o.A), sp.cancel(self.B+o.B))
    __radd__ = __add__
    def __sub__(self, o):
        o = o if isinstance(o, QExt) else QExt(o)
        return QExt(sp.cancel(self.A-o.A), sp.cancel(self.B-o.B))
    def __rsub__(self, o):
        return QExt(o) - self
    def __neg__(self):
        return QExt(-self.A, -self.B)
    def __mul__(self, o):
        if not isinstance(o, QExt):
            return QExt(sp.cancel(self.A*o), sp.cancel(self.B*o))
        A = sp.cancel(self.A*o.A + self.B*o.B*R_EXPR)
        B = sp.cancel(self.A*o.B + self.B*o.A)
        return QExt(A,B)
    __rmul__ = __mul__
    def __truediv__(self, o):
        if not isinstance(o, QExt):
            return QExt(sp.cancel(self.A/o), sp.cancel(self.B/o))
        denom = sp.cancel(o.A**2 - o.B**2*R_EXPR)
        num = self * QExt(o.A, -o.B)
        return QExt(sp.cancel(num.A/denom), sp.cancel(num.B/denom))
    def __rtruediv__(self, o):
        return QExt(o)/self
    def __repr__(self):
        return f"QExt({self.A} + ({self.B})*S)"

R_EXPR = None

def set_R(Rexpr):
    global R_EXPR
    R_EXPR = Rexpr
