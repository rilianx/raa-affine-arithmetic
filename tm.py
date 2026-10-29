"""First-order Taylor models (Makino and Berz, 2003), used as a second baseline.

A model is  T(x) = c0 + sum_i c_i (x_i - m_i) + I,  where m_i is the midpoint of the
i-th interval of the box and I = [lo, hi] is the interval remainder.

The product of two models keeps the part of the polynomial product up to order 1 and
encloses the rest by B(P_e) + B(P_1) I_2 + B(P_2) I_1 + I_1 I_2 [Makino and Berz, 2003,
Def. 2], where P_e is the second-order part and B(P) a bound of P over the box. The
second-order monomials are bounded one by one, so that the squares (x_i - m_i)^2 are
enclosed in [0, r_i^2] (the "special treatment of squares" of Neumaier, 2003). The
remainders are not re-centered. Univariate functions use the same Chebyshev
linearization as ZCAA and AF2 (zcaa/unary.py). Floating point, no outward rounding.
"""
from .box import get_box


def _imul(a, b, c, d):
    p = (a * c, a * d, b * c, b * d)
    return min(p), max(p)


class TM1:
    __slots__ = ("c0", "coef", "lo", "hi")

    def __init__(self, c0, coef, lo=0.0, hi=0.0):
        self.c0 = float(c0)
        self.coef = coef if isinstance(coef, list) else list(coef)
        self.lo, self.hi = float(lo), float(hi)

    # ---- ranges -------------------------------------------------------------
    def poly_bound(self):
        """Bound of the polynomial part over the box."""
        dev = sum(abs(c) * 0.5 * (b - a) for c, (a, b) in zip(self.coef, get_box()))
        return self.c0 - dev, self.c0 + dev

    def get_range(self):
        a, b = self.poly_bound()
        return a + self.lo, b + self.hi

    def residual_width(self):
        """Vertical distance between the two parallel envelopes (width of the remainder)."""
        return self.hi - self.lo

    def __call__(self, x):
        v = self.c0 + sum(c * (xi - 0.5 * (a + b)) for c, xi, (a, b) in zip(self.coef, x, get_box()))
        return v + self.lo, v + self.hi

    # ---- affine operations ----------------------------------------------------
    def __add__(self, o):
        if isinstance(o, TM1):
            return TM1(self.c0 + o.c0, [a + b for a, b in zip(self.coef, o.coef)],
                       self.lo + o.lo, self.hi + o.hi)
        if isinstance(o, (int, float)):
            return TM1(self.c0 + o, self.coef, self.lo, self.hi)
        return NotImplemented

    def __radd__(self, o):
        return self + o

    def scale(self, c):
        c = float(c)
        lo, hi = (c * self.lo, c * self.hi) if c >= 0 else (c * self.hi, c * self.lo)
        return TM1(self.c0 * c, [a * c for a in self.coef], lo, hi)

    def __neg__(self):
        return self.scale(-1.0)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def add_interval(self, lo, hi):
        return TM1(self.c0, self.coef, self.lo + lo, self.hi + hi)

    # ---- multiplication ---------------------------------------------------------
    def __mul__(self, o):
        if isinstance(o, (int, float)):
            return self.scale(o)
        if not isinstance(o, TM1):
            return NotImplemented
        box = get_box()
        rad = [0.5 * (b - a) for a, b in box]
        a0, b0, a, b = self.c0, o.c0, self.coef, o.coef
        coef = [a0 * bi + b0 * ai for ai, bi in zip(a, b)]
        nz = [i for i in range(len(a)) if a[i] != 0 or b[i] != 0]
        # second-order part, monomial by monomial
        elo = ehi = 0.0
        for k, i in enumerate(nz):
            c = a[i] * b[i]
            sq = rad[i] * rad[i]
            if c > 0:
                ehi += c * sq
            else:
                elo += c * sq
            for j in nz[k + 1:]:
                cc = abs(a[i] * b[j] + a[j] * b[i]) * rad[i] * rad[j]
                elo -= cc
                ehi += cc
        # remainder terms
        p1, p2 = self.poly_bound(), o.poly_bound()
        t1 = _imul(p1[0], p1[1], o.lo, o.hi)
        t2 = _imul(p2[0], p2[1], self.lo, self.hi)
        t3 = _imul(self.lo, self.hi, o.lo, o.hi)
        return TM1(a0 * b0, coef, elo + t1[0] + t2[0] + t3[0], ehi + t1[1] + t2[1] + t3[1])

    def __rmul__(self, o):
        return self * o

    def __truediv__(self, o):
        if isinstance(o, (int, float)):
            return self.scale(1.0 / o)
        raise NotImplementedError("division by a model is not implemented")

    def __pow__(self, k):
        if k == 1:
            return self
        from .unary import power
        return power(self, k)

    @property
    def slopes(self):
        """Coefficients of the linear part in the original variables."""
        return list(self.coef)

    def __repr__(self):
        return f"TM1({self.c0:g}; {self.coef}; [{self.lo:g}, {self.hi:g}])"


def tm_variables(n=None):
    """Taylor models of the variables of the current box: x_i = m_i + (x_i - m_i)."""
    box = get_box()
    n = len(box) if n is None else n
    out = []
    for i in range(n):
        c = [0.0] * n
        c[i] = 1.0
        out.append(TM1(0.5 * (box[i][0] + box[i][1]), c))
    return out
