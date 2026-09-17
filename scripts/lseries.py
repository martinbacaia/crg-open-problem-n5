"""
Serie de Laurent truncada, EXACTA, sobre Q, respaldada por `fractions.Fraction` puro (sin sympy).

Motivacion (medida, no solo esperada): un benchmark directo (multiplicar+truncar 2000 veces series
de ~20 terminos) dio sympy=106s vs. este backend=1.6s (~65x). El cuello de botella de
regge_n5_simbolico_grado9/11.py y sus analogos de familia 2 NO es algebra simbolica genuina (no hay
simbolos libres mas alla de x, ya truncado a una ventana finita antes de este paso) -- es aritmetica
racional de series ya truncadas, para la que sympy paga un overhead de arbol de expresion enorme.

Solo se necesitan +,-,*,potencia entera pequena (NINGUNA division: eta, F:F, C, Msandwich y
S1..S15/B1..B3 se construyen unicamente con sumas/restas/productos, ver los scripts originales) --
por eso no se implementa division aqui (mas simple, menos superficie de bug).

Representacion: diccionario {orden_entero: Fraction}, solo entradas no nulas.

IMPORTANTE (bug real encontrado y corregido por `test_lseries.py`): NO se trunca el orden superior
en cada multiplicacion. Un primer intento truncaba a un orden fijo `TRUNC_HI` despues de CADA `*`
(para acotar el tamano) -- pero eso es INCORRECTO en general para series de Laurent con ordenes
negativos: un termino de orden alto descartado en un paso puede, al multiplicarse mas adelante por
un factor de orden NEGATIVO, "volver a caer" dentro de la ventana de ordenes que sí importan,
dandolos MAL. El test de 200 casos aleatorios con secuencias de operaciones encontro exactamente
este tipo de discrepancia (varios ordenes intermedios, no solo cerca del corte). Como en este
proyecto todas las series base (`ps_ser`, `eps_ser`) ya vienen truncadas UNA VEZ por sympy (con un
`N` generoso) antes de convertirse a `LSeries`, todo lo que sigue (eta, F:F, C, Msandwich, S1..S15,
B1..B3) es aritmetica EXACTA de polinomios de Laurent FINITOS -- no hace falta truncar nunca mas; el
numero de ordenes distintos en cualquier producto esta acotado por la suma de los rangos de los
factores (nunca explota combinatoriamente, solo se acota por el ancho total de ordenes involucrados).
"""
from fractions import Fraction

def _to_fraction(v):
    if isinstance(v, Fraction):
        return v
    if isinstance(v, int):
        return Fraction(v)
    # sympy Rational/Integer: tiene .p, .q
    if hasattr(v, 'p') and hasattr(v, 'q'):
        return Fraction(int(v.p), int(v.q))
    # expresion sympy general (p.ej. Mul de Rational*Integer sin auto-combinar): forzar a
    # numerador/denominador exactos via together() antes de convertir.
    # IMPORTANTE (bug real encontrado: int() sobre un sympy irracional, p.ej. int(sqrt(3))==1,
    # trunca SILENCIOSAMENTE en vez de fallar -- eso corrompio sin avisar un cross-check anterior
    # cuando R(x) resulto tener sqrt(3) en la serie de S antes de cancelar en la suma final sobre
    # S5). Por eso se exige explicitamente que num/den sean enteros sympy antes de convertir.
    import sympy as sp
    num, den = sp.together(sp.sympify(v)).as_numer_denom()
    if not (num.is_Integer and den.is_Integer):
        raise ValueError(
            f"_to_fraction: coeficiente no es racional exacto (contiene un irracional sin "
            f"cancelar, p.ej. sqrt): {v!r} -> num={num!r}, den={den!r}. Esto NO se trunca "
            f"silenciosamente -- si esto aparece, la serie base (antes de la suma sobre la orbita "
            f"de S5) todavia tiene una raiz sin cancelar; hay que revisar los parametros de la "
            f"cinematica (t0, m0sq, etc.) para que R(x) de una serie sqrt(R) racional, como ya se "
            f"hace en el resto del proyecto."
        )
    return Fraction(int(num), int(den))


class LSeries:
    __slots__ = ('c',)

    def __init__(self, c=None):
        self.c = c if c is not None else {}

    @staticmethod
    def zero():
        return LSeries({})

    @staticmethod
    def constant(v):
        v = _to_fraction(v)
        return LSeries({0: v}) if v != 0 else LSeries({})

    @staticmethod
    def from_sympy(expr, x, lo, hi):
        """Extrae coeficientes de una expresion sympy (ya truncada/expandida) para ordenes en
        [lo, hi) y los convierte a Fraction exactos."""
        import sympy as sp
        expr = sp.expand(expr)
        d = {}
        for o in range(lo, hi):
            co = expr.coeff(x, o)
            if co != 0:
                fr = _to_fraction(co)
                if fr != 0:
                    d[o] = fr
        return LSeries(d)

    def coeff(self, o):
        return self.c.get(o, Fraction(0))

    def is_zero(self):
        return len(self.c) == 0

    def lowest_order(self):
        if not self.c:
            return None
        return min(self.c.keys())

    def __eq__(self, other):
        if isinstance(other, (int, Fraction)):
            other = LSeries.constant(other)
        return self.c == other.c

    def __add__(self, other):
        if isinstance(other, (int, Fraction)):
            other = LSeries.constant(other)
        d = dict(self.c)
        for o, v in other.c.items():
            nv = d.get(o, Fraction(0)) + v
            if nv == 0:
                d.pop(o, None)
            else:
                d[o] = nv
        return LSeries(d)

    __radd__ = __add__

    def __neg__(self):
        return LSeries({o: -v for o, v in self.c.items()})

    def __sub__(self, other):
        if isinstance(other, (int, Fraction)):
            other = LSeries.constant(other)
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            fo = _to_fraction(other)
            return LSeries({o: v * fo for o, v in self.c.items() if v * fo != 0})
        d = {}
        for oa, va in self.c.items():
            for ob, vb in other.c.items():
                o = oa + ob
                nv = d.get(o, Fraction(0)) + va * vb
                if nv == 0:
                    d.pop(o, None)
                else:
                    d[o] = nv
        return LSeries(d)

    __rmul__ = __mul__

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        result = LSeries.constant(1)
        base = self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def __repr__(self):
        terms = sorted(self.c.items())
        return "LSeries(" + " + ".join(f"({v})*x^{o}" for o, v in terms) + ")"


def eta_vec(a, b, D):
    """a,b: listas de D componentes LSeries (o convertibles). Metrica mostly-plus: -a0 b0 + sum a_k b_k."""
    tot = -(a[0] * b[0])
    for k in range(1, D):
        tot = tot + a[k] * b[k]
    return tot
