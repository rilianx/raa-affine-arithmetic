"""Univariate functions for ZCAA and AF2 forms.

ZCAA, AF2 and the first-order Taylor models use the same Chebyshev linearization u(v) ~ m v + [r_lo, r_hi]
over the range [a, b] of the argument (Section 3.2 of the paper), with the same
rule for the slope m, so that the comparisons isolate the effect of the
multiplication:
  - if u does not change curvature on [a, b]: secant slope;
  - otherwise: the better (narrower residual) of the secant slope and the slope of
    the tangent at the inflection point.
The residual is bounded by evaluating u(t) - m t at the endpoints and at the
critical points where u'(t) = m. Plain numbers are evaluated with the math module.
"""
import math

from .zcaa import ZCAAForm
from .af2 import AF2
from .tm import TM1

SAFETY = 1e-10


def _is_number(x):
    return isinstance(x, (int, float))


def _residual(f, m, a, b, crit):
    pts = [a, b] + [t for t in crit(m, a, b) if a <= t <= b]
    vals = [f(t) - m * t for t in pts]
    return min(vals), max(vals)


def chebyshev(form, f, df, crit, inflection=None):
    """Linear enclosure of f(form). crit(m, a, b) returns the points of [a, b] where f'(t) = m."""
    a, b = form.get_range()
    if b - a < 1e-15:                       # point range: constant result
        m = df(a)
        rlo = rhi = f(a) - m * a
    elif inflection is not None and a < inflection < b:
        m1 = (f(b) - f(a)) / (b - a)
        lo1, hi1 = _residual(f, m1, a, b, crit)
        m2 = df(inflection)
        lo2, hi2 = _residual(f, m2, a, b, crit)
        if hi2 - lo2 < hi1 - lo1:
            m, rlo, rhi = m2, lo2, hi2
        else:
            m, rlo, rhi = m1, lo1, hi1
        rlo, rhi = rlo - SAFETY, rhi + SAFETY
    else:
        m = (f(b) - f(a)) / (b - a)
        rlo, rhi = _residual(f, m, a, b, crit)
        if abs(rlo) < 1e-9:
            rlo = 0.0
        if abs(rhi) < 1e-9:
            rhi = 0.0
        rlo, rhi = rlo - SAFETY, rhi + SAFETY

    if isinstance(form, ZCAAForm):
        elo, ehi = (m * form.lo, m * form.hi) if m >= 0 else (m * form.hi, m * form.lo)
        return ZCAAForm([s * m for s in form.slopes], rlo + elo, rhi + ehi)
    if isinstance(form, AF2):
        res = form.scale(m) + 0.5 * (rlo + rhi)
        res.s += 0.5 * (rhi - rlo)
        return res
    if isinstance(form, TM1):
        # constant part of the residual to the polynomial, symmetric remainder
        res = form.scale(m) + 0.5 * (rlo + rhi)
        return res.add_interval(-0.5 * (rhi - rlo), 0.5 * (rhi - rlo))
    raise TypeError(f"unsupported form: {type(form).__name__}")


# ---- functions ---------------------------------------------------------------
def exp(x):
    if _is_number(x):
        return math.exp(x)
    return chebyshev(x, math.exp, math.exp, lambda m, a, b: [math.log(m)] if m > 0 else [])


def log(x):
    if _is_number(x):
        return math.log(x)
    if x.get_range()[0] <= 0:
        raise ValueError("log: the range of the argument is not positive")
    return chebyshev(x, math.log, lambda t: 1.0 / t, lambda m, a, b: [1.0 / m] if m > 0 else [])


def sqrt(x):
    if _is_number(x):
        return math.sqrt(x)
    if x.get_range()[0] < 0:
        raise ValueError("sqrt: the range of the argument contains negative values")
    return chebyshev(x, math.sqrt, lambda t: 0.5 / math.sqrt(t) if t > 0 else math.inf,
                     lambda m, a, b: [1.0 / (4.0 * m * m)] if m > 0 else [])


def sinh(x):
    if _is_number(x):
        return math.sinh(x)
    crit = lambda m, a, b: [-math.acosh(m), math.acosh(m)] if m >= 1 else []
    return chebyshev(x, math.sinh, math.cosh, crit, inflection=0.0)


def power(x, k):
    """x^k for an integer k >= 2."""
    if _is_number(x):
        return x ** k
    if not (isinstance(k, int) and k >= 2):
        raise NotImplementedError(f"power: exponent {k} is not supported")

    def crit(m, a, b):
        q = m / k
        if k % 2 == 0:                     # k - 1 odd: one real root
            return [math.copysign(abs(q) ** (1.0 / (k - 1)), q)]
        if q < 0:
            return []
        t = q ** (1.0 / (k - 1))
        return [-t, t]

    return chebyshev(x, lambda t: t ** k, lambda t: k * t ** (k - 1), crit,
                     inflection=0.0 if k % 2 == 1 else None)


def _periods(a, b):
    return range(math.floor(a / (2 * math.pi)) - 1, math.ceil(b / (2 * math.pi)) + 2)


def _sin_crit(m, a, b):
    if abs(m) > 1:
        return []
    base = math.acos(m)                      # sin'(t) = cos(t) = m
    return [s * base + 2 * math.pi * j for j in _periods(a, b) for s in (1, -1)]


def _cos_crit(m, a, b):
    if abs(m) > 1:
        return []
    base = math.asin(-m)                     # cos'(t) = -sin(t) = m
    return [p + 2 * math.pi * j for j in _periods(a, b) for p in (base, math.pi - base)]


def sin(x):
    if _is_number(x):
        return math.sin(x)
    a = x.get_range()[0]
    infl = math.ceil(a / math.pi) * math.pi                    # first inflection >= a
    return chebyshev(x, math.sin, math.cos, _sin_crit, inflection=infl)


def cos(x):
    if _is_number(x):
        return math.cos(x)
    a = x.get_range()[0]
    infl = math.ceil((a - math.pi / 2) / math.pi) * math.pi + math.pi / 2
    return chebyshev(x, math.cos, lambda t: -math.sin(t), _cos_crit, inflection=infl)
