"""Sampling checks of the enclosures and regression tests.

    python -m pytest tests          # or: python tests/test_soundness.py

For random boxes, every test expression is evaluated at random points of the box,
and the value must lie between the lower and upper envelopes of ZCAA and AF2.
This is a numerical check, not a proof (the arithmetic has no outward rounding).
"""
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from zcaa import set_box, zcaa_variables, af2_variables, tm_variables, AF2, ZCAAForm, mul_zcaa   # noqa: E402
from benchmarks.functions import ALL, SIZES                                     # noqa: E402

BOXES, POINTS = 20, 200


def _box(f, size):
    wmin, wmax = SIZES[size]
    lo, hi = f.domain
    out = []
    for _ in range(3):
        w = min(random.uniform(wmin, wmax), hi - lo)
        s = random.uniform(lo, hi - w)
        out.append((s, s + w))
    return out


def test_enclosures():
    random.seed(0)
    for f in ALL:
        for size in f.sizes:
            for _ in range(BOXES):
                box = _box(f, size)
                set_box(box)
                r = f(*zcaa_variables(3))
                a = f(*af2_variables(box))
                t = f(*tm_variables(3))
                for _ in range(POINTS):
                    x = [random.uniform(lo, hi) for lo, hi in box]
                    v = f(*x)
                    tol = 1e-7 * max(1.0, abs(v))
                    lo, hi = r(x)
                    assert lo - tol <= v <= hi + tol, (f.name, "ZCAA", box, x)
                    lo, hi = a.envelopes(x, box)
                    assert lo - tol <= v <= hi + tol, (f.name, "AF2", box, x)
                    lo, hi = t(x)
                    assert lo - tol <= v <= hi + tol, (f.name, "TM", box, x)


def test_zcaa_never_looser_than_af2_and_tm():
    random.seed(1)
    for f in ALL:
        for size in f.sizes:
            for _ in range(BOXES):
                box = _box(f, size)
                set_box(box)
                wr = f(*zcaa_variables(3)).residual_width()
                wa = f(*af2_variables(box)).residual_width()
                wt = f(*tm_variables(3)).residual_width()
                assert wr <= wa * (1 + 1e-7) + 1e-9, (f.name, size, wr, wa)
                assert wr <= wt * (1 + 1e-7) + 1e-9, (f.name, size, wr, wt)


def test_af2_independent_error_symbols():
    # (x1^2 - 0.5)(0.5 - x2^2) over [-1, 1]^2 ranges over [-0.25, 0.25]; the product of
    # the two e(+-) terms must be enclosed in [-1, 1], not in [0, 1].
    x1 = AF2(0.0, {1: 1.0})
    x2 = AF2(0.0, {2: 1.0})
    lo, hi = ((x1 ** 2 - 0.5) * (0.5 - x2 ** 2)).get_range()
    assert lo <= -0.25 + 1e-12 and hi >= 0.25 - 1e-12


def test_af2_negative_scaling_swaps_one_sided_terms():
    x = AF2(1.0, {1: 1.0}, s=0.0, p=0.5, n=0.0)        # range [0, 2.5]
    lo, hi = (-2 * x).get_range()
    assert abs(lo + 5.0) < 1e-12 and abs(hi) < 1e-12


def test_af2_add_constant_is_not_in_place():
    x = AF2(1.0, {1: 1.0})
    _ = x + 3.0
    assert x.base == 1.0


def test_lemma_example():
    # Example of Section 3.4: v = x1 + 2 x2 + [-1, 1], w = 2 x1 + 3 x2 + [0, 2], x in [0, 2]^2
    set_box([(0, 2), (0, 2)])
    u = mul_zcaa(ZCAAForm([1.0, 2.0], -1, 1), ZCAAForm([2.0, 3.0], 0, 2))
    assert [round(s, 12) for s in u.slopes] == [12.0, 21.0]
    assert abs(u.lo + 38) < 1e-12 and abs(u.hi - 18) < 1e-12


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
