# ZCAA: a symbolic-free affine arithmetic with zero-centered residuals

This repository contains the implementation and the experiments of the paper

> I. Araya, A. Alvear Tapia, D. Miranda, J. Aramayo, V. Reyes.
> *A Symbolic-Free Affine Arithmetic with Zero-Centered Residuals.* Submitted.

**Zero-centered affine arithmetic (ZCAA)** represents every subexpression of a factorable
function as a linear function of the input variables plus a bounded residual interval,

    f(x) ∈ a_1 x_1 + ... + a_n x_n + [r_lo, r_hi]     for all x in the box,

without the symbolic noise variables of classical affine arithmetic. The two parallel
hyperplanes give a range enclosure and a linear relaxation of f over the box. Products of
forms are decomposed into linear, quadratic and bilinear terms, each enclosed by its own
envelope, and the midpoints of the residuals are absorbed into the linear part at every
multiplication, so that the residuals propagated to the next operations are centered at zero.
The experiments compare ZCAA with the AF2 variant of affine arithmetic (Messine, 2002) and with
first-order Taylor models (Makino and Berz, 2003).

## Repository layout

```
zcaa/                   the three methods
  zcaa.py               ZCAA forms, square and product (Sections 3.1-3.4 of the paper)
  af2.py                AF2 baseline (see "Notes on the AF2 baseline" below)
  tm.py                 first-order Taylor model baseline (TM)
  unary.py              sin, cos, exp, log, sqrt, sinh and integer powers, with the same
                        Chebyshev linearization for the three methods
  box.py                input box shared by the ZCAA forms
benchmarks/
  functions.py          test expressions: Cases 1-4 (Table 1) and the benchmark constraints
  instances/            original benchmark files (Ibex/minibex format) of these constraints
experiments/
  run_experiments.py    residual widths of ZCAA, AF2 and TM and optimal width with IbexOpt
  make_tables.py        Tables 2 and 3, win/tie/loss counts and Wilcoxon tests
  make_figure.py        Figure 2
  timing.py             running time of ZCAA, AF2 and TM
tests/
  test_soundness.py     sampling checks of the enclosures and regression tests
results/                results used in the paper (JSON), tables (LaTeX) and figure
```

## Installation

```
git clone https://github.com/rilianx/raa-affine-arithmetic.git
cd raa-affine-arithmetic
pip install -r requirements.txt
```

The scripts are run from the root of the repository; no installation of the package is needed.

## Requirements

- Python 3.9 or later with `numpy`, `scipy` and `matplotlib`.
  The arithmetics themselves use only the Python standard library.
- [Ibex](https://github.com/ibex-team/ibex-lib) with the IbexOpt optimizer, only to compute
  the optimal widths used as reference in the tables (`run_experiments.py`). Install Ibex
  following its documentation and set the path to the `ibexopt` executable:

  ```
  export IBEXOPT=/path/to/ibexopt
  ```

  The results used in the paper are included in `results/`, so the tables and the figure can
  be regenerated without Ibex.

## Quick start

```python
from zcaa import set_box, zcaa_variables, af2_variables, sin

box = [(0.0, 2.0), (1.0, 3.0)]
set_box(box)
x1, x2 = zcaa_variables()

f = (x1 + 2 * x1 * x2) * (3 * x2 - x1 * x2) + sin(x1 * x2)
print(f)                     # linear part and residual interval
print(f.get_range())         # enclosure of the range: (-43.53, 81.35)
print(f.residual_width())    # distance between the two envelopes: 68.41

y1, y2 = af2_variables(box)  # the same expression with AF2
g = (y1 + 2 * y1 * y2) * (3 * y2 - y1 * y2) + sin(y1 * y2)
print(g.get_range())         # (-63.53, 97.35)
print(g.residual_width())    # 104.41
```

## Reproducing the results of the paper

```
python tests/test_soundness.py                      # checks (about 1 min)
python experiments/run_experiments.py controlled    # Table 2 data (about 1 min)
python experiments/run_experiments.py instances     # Table 3 data (1-2 h, mostly IbexOpt)
python experiments/make_tables.py --latex           # Tables 2 and 3 and the statistics
python experiments/make_figure.py                   # Figure 2
python experiments/timing.py                        # running times
```

`run_experiments.py --runs N` uses only the first N boxes of each configuration, with the same
random sequence as the full run, for a quick check. The seeds, numbers of runs and IbexOpt
settings of the paper are fixed in `run_experiments.py`:

| Set | Expressions | Boxes per size | Sizes (width of each interval) | IbexOpt precision |
|---|---|---|---|---|
| controlled | Cases 1-4 | 100 | small [0.01, 0.1], medium [0.1, 1], large [1, 10] | 1e-8 |
| instances | Trigexp1 | 100 | small, medium, large | 1e-8 |
| instances | Trigexp2 (1), Trigexp2 (2) | 50 | small, medium, large | 1e-6 |
| instances | Pramanik | 50 | small, medium | 1e-6 |

The boxes are sampled in [-10, 10]^3, except for Pramanik, sampled in its domain [0.06, 1]^3.
For every box, the gap of a method M is (W_M - W_OPT) / W_OPT, where W_M is the residual width
of M (the vertical distance between its two parallel envelopes) and W_OPT is the exact width
of the range of f(x) - L(x), L being the linear part of the ZCAA enclosure, computed with
IbexOpt. The same W_OPT is used for both methods. Two widths that differ by less than a
relative 1e-9 are counted as equal.

Running times depend on the machine and on the Python version; the relevant quantity is the
ratio between the methods.

## Benchmark constraints

The constraints of `benchmarks/functions.py` come from the benchmark of systems of
equations used to test interval solvers (COPRIN examples, distributed with Ibex); the original
files are in `benchmarks/instances/`. We selected constraints with three variables that
contain products of two factors depending on at least two variables each:

- `Trigexp1-020.bch`, constraint 2 (variables x1, x2, x3);
- `Trigexp2-3.bch`, constraints 1 and 2;
- `Pramanik.bch`, constraint 1, with the constants f0, f1, y0, y1 replaced by their values.

## Notes on the Taylor model baseline

`zcaa/tm.py` implements first-order Taylor models T(x) = c0 + sum_i c_i (x_i - m_i) + I, centered
at the midpoint m of the box. The product keeps the first-order part and encloses the rest by
B(P_e) + B(P_1) I_2 + B(P_2) I_1 + I_1 I_2 (Makino and Berz, 2003, Def. 2), bounding the
second-order monomials one by one, so that the squares (x_i - m_i)^2 are enclosed in [0, r_i^2].
Univariate functions use the same Chebyshev linearization as ZCAA and AF2, with the constant part
of the residual moved to c0, so that the remainder they produce is centered. The remainders
produced by the products are not re-centered.

## Notes on the AF2 baseline

`zcaa/af2.py` implements AF2 as described by Messine (2002). It corrects four issues of the
implementation used in an earlier version of our experiments, which could produce enclosures
that do not contain the range of the expression:

1. `x + c` modified `x` in place;
2. `c * x` with `c < 0`, and `-x`, did not swap the e(+) and e(-) terms;
3. in `x * y`, the e(+) contributions of the centers were omitted when both centers were
   negative (and all of them when a center was zero);
4. in `x * y`, the product of the two e(+-) terms was enclosed in [0, 1]; since the two
   symbols are independent, it is enclosed in [-1, 1].

For example, over x1, x2 in [-1, 1], the product (x1^2 - 0.5)(0.5 - x2^2) ranges over
[-0.25, 0.25]; the earlier implementation returned [0, 0.25]. `tests/test_soundness.py`
contains regression tests for these cases.

## Numerical remarks

Both arithmetics are implemented in floating point **without outward rounding**: the
enclosures are correct up to rounding errors, and they are checked numerically by sampling
(`tests/test_soundness.py`), not rigorously. A rigorous implementation would round every
bound outward, for example with an interval library.

## Citation

If you use this code, please cite the paper (the reference will be updated upon publication):

```bibtex
@unpublished{araya2026zcaa,
  title  = {A Symbolic-Free Affine Arithmetic with Zero-Centered Residuals},
  author = {Araya, Ignacio and Alvear Tapia, Angelo and Miranda, Daniel and
            Aramayo, Juan and Reyes, Victor},
  note   = {Submitted},
  year   = {2026}
}
```

## Acknowledgements

This work was supported by ANID, FONDECYT Regular 2026, project No. 1260159.

## License

MIT, see `LICENSE`.
