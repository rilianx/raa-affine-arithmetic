"""AF2 affine arithmetic (Messine, 2002), used as the baseline.

A form is  x0 + sum_i x_i e_i + s e(+-) + p e(+) + n e(-),  with e_i, e(+-) in [-1, 1],
e(+) in [0, 1], e(-) in [-1, 0] and s, p, n >= 0. The three last terms collect the
errors of the non-affine operations.

This implementation follows the one used in an earlier version of our experiments,
with the following corrections (see the README):
  1. x + c no longer modifies x in place;
  2. c * x with c < 0, and -x, swap the e(+) and e(-) terms;
  3. in x * y, the e(+) and e(-) contributions of the centers are assigned by their
     sign in all cases (the case x0 < 0, y0 < 0 and the case x0 = 0 were missing);
  4. in x * y, the product of the two e(+-) terms is enclosed in [-1, 1] (e(+-)),
     since the two symbols are independent, instead of [0, 1].
Arithmetic is in floating point, without outward rounding.
"""
import math


class AF2:
    __slots__ = ("base", "eps", "s", "p", "n")

    def __init__(self, base, eps=None, s=0.0, p=0.0, n=0.0):
        self.base = float(base)
        self.eps = dict(eps) if eps else {}
        self.s, self.p, self.n = float(s), float(p), float(n)

    # ---- ranges -------------------------------------------------------------
    def get_range(self):
        dev = sum(abs(c) for c in self.eps.values()) + self.s
        return self.base - dev - self.n, self.base + dev + self.p

    def residual_width(self):
        """Width of the non-affine part, 2 s + p + n."""
        return 2.0 * self.s + self.p + self.n

    def envelopes(self, x, box):
        """Lower and upper envelopes at the point x of the box (e_i = (x_i - mid_i) / rad_i)."""
        val = self.base
        for k, c in self.eps.items():
            a, b = box[k - 1]
            rad = 0.5 * (b - a)
            if rad > 0:
                val += c * (x[k - 1] - 0.5 * (a + b)) / rad
        return val - self.s - self.n, val + self.s + self.p

    # ---- affine operations ----------------------------------------------------
    def _num(self, c):
        return isinstance(c, (int, float))

    def __add__(self, other):
        if isinstance(other, AF2):
            eps = dict(self.eps)
            for k, c in other.eps.items():
                eps[k] = eps.get(k, 0.0) + c
            return AF2(self.base + other.base, eps, self.s + other.s,
                       self.p + other.p, self.n + other.n)
        if self._num(other):
            return AF2(self.base + other, self.eps, self.s, self.p, self.n)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def scale(self, c):
        c = float(c)
        eps = {k: v * c for k, v in self.eps.items()}
        if c >= 0:
            return AF2(self.base * c, eps, self.s * c, self.p * c, self.n * c)
        a = -c   # c e(+) <= 0 goes to e(-), c e(-) >= 0 goes to e(+)
        return AF2(self.base * c, eps, self.s * a, self.n * a, self.p * a)

    def __neg__(self):
        return self.scale(-1.0)

    def __sub__(self, other):
        if isinstance(other, AF2):
            eps = dict(self.eps)
            for k, c in other.eps.items():
                eps[k] = eps.get(k, 0.0) - c
            return AF2(self.base - other.base, eps, self.s + other.s,
                       self.p + other.n, self.n + other.p)
        if self._num(other):
            return AF2(self.base - other, self.eps, self.s, self.p, self.n)
        return NotImplemented

    def __rsub__(self, other):
        if self._num(other):
            return self.scale(-1.0) + other
        return NotImplemented

    # ---- multiplication -------------------------------------------------------
    def _k1(self, o):
        """Symmetric terms, in e(+-)."""
        k = abs(self.base) * o.s + abs(o.base) * self.s
        for i, xi in self.eps.items():
            for j, yj in o.eps.items():
                if i != j:
                    k += abs(xi * yj)
            k += abs(xi) * (o.s + o.p + o.n)
        for yj in o.eps.values():
            k += abs(yj) * (self.s + self.p + self.n)
        k += self.s * (o.s + o.p + o.n)
        k += (self.p + self.n) * o.s
        k += self.p * o.n + self.n * o.p
        return k

    def _k2(self, o):
        """Terms in [0, 1], in e(+)."""
        x0, y0 = self.base, o.base
        k = x0 * o.p if x0 > 0 else -x0 * o.n
        k += y0 * self.p if y0 > 0 else -y0 * self.n
        for i, xi in self.eps.items():
            yi = o.eps.get(i, 0.0)
            if xi * yi > 0:
                k += xi * yi
        k += self.p * o.p + self.n * o.n
        return k

    def _k3(self, o):
        """Terms in [-1, 0], in e(-)."""
        x0, y0 = self.base, o.base
        k = x0 * o.n if x0 > 0 else -x0 * o.p
        k += y0 * self.n if y0 > 0 else -y0 * self.p
        for i, xi in self.eps.items():
            yi = o.eps.get(i, 0.0)
            if xi * yi < 0:
                k += abs(xi * yi)
        return k

    def __mul__(self, other):
        if self._num(other):
            return self.scale(other)
        if not isinstance(other, AF2):
            return NotImplemented
        eps = {}
        for k in set(self.eps) | set(other.eps):
            eps[k] = self.base * other.eps.get(k, 0.0) + self.eps.get(k, 0.0) * other.base
        return AF2(self.base * other.base, eps, self._k1(other), self._k2(other), self._k3(other))

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if self._num(other):
            return self.scale(1.0 / other)
        raise NotImplementedError("division by a form is not implemented")

    def __pow__(self, k):
        if k == 1:
            return self
        if k == 2:
            return self._square()
        from .unary import power
        return power(self, k)

    def _square(self):
        """Chebyshev approximation of x^2 over the range [a, b] of the form."""
        a, b = self.get_range()
        if b - a < 1e-15:
            return AF2(a * a)
        alpha = a + b
        u = alpha / 2.0
        r_u = alpha * u - alpha * a + a * a      # secant at u
        f_u = u * u
        res = self.scale(alpha) + (((f_u + r_u) / 2.0) - alpha * u)
        res.s += abs(f_u - r_u) / 2.0
        return res

    def __repr__(self):
        terms = " ".join(f"{c:+g}e{k}" for k, c in sorted(self.eps.items()))
        return f"AF2({self.base:g} {terms} +{self.s:g}e(+-) +{self.p:g}e(+) +{self.n:g}e(-))"


def af2_variables(box):
    """AF2 forms of the variables of the box: x_i = mid_i + rad_i e_i."""
    return [AF2(0.5 * (a + b), {i + 1: 0.5 * (b - a)}) for i, (a, b) in enumerate(box)]
