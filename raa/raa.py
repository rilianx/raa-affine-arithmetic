"""Revised affine arithmetic (RAA).

A form is f(x) = sum_i a_i x_i + r, with r in [lo, hi], over the box set with
raa.set_box. The residual includes the intercept and is not assumed to be
centered. Arithmetic is in floating point, without outward rounding.
"""
from .box import get_box


class RAAForm:
    """f(x) = sum(slopes[i] * x_i) + [lo, hi]."""
    __slots__ = ("slopes", "lo", "hi")

    def __init__(self, slopes, lo, hi=None):
        self.slopes = slopes if isinstance(slopes, list) else list(slopes)
        self.lo = float(lo)
        self.hi = float(lo if hi is None else hi)

    # ---- ranges -------------------------------------------------------------
    def linear_range(self):
        lo = hi = 0.0
        for s, (a, b) in zip(self.slopes, get_box()):
            if s > 0:
                lo += s * a
                hi += s * b
            elif s < 0:
                lo += s * b
                hi += s * a
        return lo, hi

    def get_range(self):
        """Enclosure of the range of the form over the box."""
        lo, hi = self.linear_range()
        return lo + self.lo, hi + self.hi

    def residual_width(self):
        """Vertical distance between the two parallel envelopes."""
        return self.hi - self.lo

    def __call__(self, x):
        """Lower and upper envelopes evaluated at the point x."""
        lin = sum(s * xi for s, xi in zip(self.slopes, x))
        return lin + self.lo, lin + self.hi

    # ---- affine operations (exact) ------------------------------------------
    def __add__(self, other):
        if isinstance(other, RAAForm):
            return RAAForm([s + o for s, o in zip(self.slopes, other.slopes)],
                           self.lo + other.lo, self.hi + other.hi)
        if isinstance(other, (int, float)):
            return RAAForm(self.slopes, self.lo + other, self.hi + other)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return RAAForm([-s for s in self.slopes], -self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def scale(self, c):
        c = float(c)
        lo, hi = (c * self.lo, c * self.hi) if c >= 0 else (c * self.hi, c * self.lo)
        return RAAForm([s * c for s in self.slopes], lo, hi)

    # ---- nonlinear operations -----------------------------------------------
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return self.scale(other)
        if isinstance(other, RAAForm):
            return mul_raa(self, other)
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return self.scale(1.0 / other)
        raise NotImplementedError("division by a form is not implemented")

    def __pow__(self, k):
        if k == 2:
            return pow_raa(self)
        from .unary import power
        return power(self, k)

    def __repr__(self):
        terms = " + ".join(f"{a:g}*x{i + 1}" for i, a in enumerate(self.slopes) if a != 0) or "0"
        return f"RAAForm({terms} + [{self.lo:g}, {self.hi:g}])"


def raa_variables(n=None):
    """RAA forms of the variables x_1..x_n of the current box."""
    n = len(get_box()) if n is None else n
    out = []
    for i in range(n):
        s = [0.0] * n
        s[i] = 1.0
        out.append(RAAForm(s, 0.0, 0.0))
    return out


def pow_raa(form):
    """Square of a form: Chebyshev slope a+b over its range [a, b] and residual
    [-(a+b)^2/4, -ab] (Section 3.2 of the paper)."""
    a, b = form.get_range()
    if b - a < 1e-15:
        v = a * a
        return RAAForm([0.0] * len(form.slopes), v, v)
    m = a + b
    rlo, rhi = -m * m / 4.0, -a * b
    elo, ehi = (m * form.lo, m * form.hi) if m >= 0 else (m * form.hi, m * form.lo)
    return RAAForm([s * m for s in form.slopes], rlo + elo, rhi + ehi)


def mul_raa(v, w):
    """Product of two RAA forms (Section 3.4, Eq. 17 of the paper).

    The midpoints of the residuals are absorbed into the linear part, every
    quadratic and bilinear monomial is enclosed by its own envelope, and the
    residual of the cross terms (Eq. 16) is evaluated at the four vertices.
    """
    box = get_box()
    a, b = v.slopes, w.slopes
    lv, hv, lw, hw = v.lo, v.hi, w.lo, w.hi
    mv, mw = 0.5 * (lv + hv), 0.5 * (lw + hw)
    nz = [i for i in range(len(a)) if a[i] != 0 or b[i] != 0]

    # 1. linear part (midpoints of the residuals absorbed)
    slopes = [mw * ai + mv * bi for ai, bi in zip(a, b)]
    elo = ehi = 0.0

    # 2. quadratic terms a_i b_i x_i^2 (Chebyshev envelope)
    for i in nz:
        c = a[i] * b[i]
        if c == 0:
            continue
        p, q = box[i]
        if q - p < 1e-15:
            elo += c * p * p
            ehi += c * p * p
            continue
        m = p + q
        rlo, rhi = -m * m / 4.0, -p * q
        slopes[i] += c * m
        if c > 0:
            elo += c * rlo
            ehi += c * rhi
        else:
            elo += c * rhi
            ehi += c * rlo

    # 3. bilinear terms (a_i b_j + a_j b_i) x_i x_j (tangent plane at the center)
    for k, i in enumerate(nz):
        pi, qi = box[i]
        for j in nz[k + 1:]:
            c = a[i] * b[j] + a[j] * b[i]
            if c == 0:
                continue
            pj, qj = box[j]
            mi, mj = 0.5 * (pj + qj), 0.5 * (pi + qi)
            r = (qi - pi) * (qj - pj) / 4.0
            shift = mi * mj
            slopes[i] += c * mi
            slopes[j] += c * mj
            if c > 0:
                elo += c * (-r - shift)
                ehi += c * (r - shift)
            else:
                elo += c * (r - shift)
                ehi += c * (-r - shift)

    # 4. residual of the cross terms, Eq. (16), evaluated at the four vertices
    glo, ghi = float("inf"), float("-inf")
    for rv in (lv, hv):
        dv = rv - mv
        for rw in (lw, hw):
            dw = rw - mw
            lo = hi = rv * rw
            for i in nz:
                c = dw * a[i] + dv * b[i]
                if c > 0:
                    lo += c * box[i][0]
                    hi += c * box[i][1]
                elif c < 0:
                    lo += c * box[i][1]
                    hi += c * box[i][0]
            glo = min(glo, lo)
            ghi = max(ghi, hi)

    return RAAForm(slopes, elo + glo, ehi + ghi)
