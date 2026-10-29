"""Test expressions of the paper.

Each expression is written once, in Ibex syntax (^ for powers, x_1, x_2, x_3 for
the variables), and evaluated in Python with RAA forms, AF2 forms or numbers.

  - Cases 1-4: controlled products (Table 1 of the paper).
  - Constraints with products of nonlinear factors taken from the benchmark of
    nonlinear systems used to test interval solvers (instances/*.bch).
"""
import re

from raa import sin, cos, exp, log, sqrt, sinh, power

_NAMESPACE = {"sin": sin, "cos": cos, "exp": exp, "log": log, "sqrt": sqrt,
              "sinh": sinh, "power": power}


class TestFunction:
    def __init__(self, name, expr, domain=(-10.0, 10.0), sizes=("small", "medium", "large"),
                 description=""):
        self.name = name
        self.expr = expr                     # Ibex syntax
        self.domain = domain                 # boxes are sampled in domain^3
        self.sizes = sizes
        self.description = description
        self._code = compile(expr.replace("^", "**"), name, "eval")

    def __call__(self, x_1, x_2, x_3):
        return eval(self._code, _NAMESPACE, {"x_1": x_1, "x_2": x_2, "x_3": x_3})

    def ibex_expr(self):
        """Expression with Ibex vector indexing x[0], x[1], x[2]."""
        return re.sub(r"x_(\d+)", lambda m: f"x[{int(m.group(1)) - 1}]", self.expr)


# ---- controlled cases (Table 1) ------------------------------------------------
CASES = [
    TestFunction("Case 1", "(x_1+2*x_2+3*x_3)*(4*x_1+5*x_2+6*x_3)",
                 description="product of two linear expressions"),
    TestFunction("Case 2", "(x_1+2*x_2^2+3*x_3)*(4*x_1+5*x_2+6*x_3^2)",
                 description="product of expressions with quadratic terms"),
    TestFunction("Case 3", "(x_1 + 2*x_2*x_3)*(3*x_2 + 4*x_1*x_3)",
                 description="product of expressions with bilinear terms"),
    TestFunction("Case 4", "(x_1+2*x_2)*(3*x_2+4*x_3)*(5*x_3+6*x_1)",
                 description="product of three linear expressions"),
]

# ---- constraints from the benchmark of nonlinear systems -------------------------
_PRAMANIK_CONST = {
    "f0": "1.7461756494150842271", "f1": "2.0364691127919609051",
    "y0": "1.3954170041747090114", "y1": "1.7444828545735749268",
}
_pramanik = (
    "-(-(x_1*(-x_3 + x_2*cos(y1))*(-x_3 + x_2*sin(f1))) + x_1*(1 + x_2*cos(f1))*(-x_3 + x_2*sin(y1)))^2"
    " + ((x_2*(-cos(f0) + cos(f1)) - x_2*x_3*(-sin(f0) + sin(f1)) - x_1*(-x_3 + x_2*sin(f1)))*(-x_3 + x_2*sin(y1))"
    " - (sin(f1)*x_2 - x_3)*(x_1*x_3 + x_2*cos(y0) - x_2*cos(y1) + x_2*(-x_1 + x_3)*sin(y0) - x_2*x_3*sin(y1)))^2"
    " + (-((-1 + x_2*cos(y1))*(x_2*(-cos(f0) + cos(f1)) - x_2*x_3*(-sin(f0) + sin(f1)) - x_1*(-x_3 + x_2*sin(f1))))"
    " + (1 + x_2*cos(f1))*(x_1*x_3 + x_2*cos(y0) - x_2*cos(y1) + x_2*(-x_1 + x_3)*sin(y0) - x_2*x_3*sin(y1)))^2"
)
for _k, _v in _PRAMANIK_CONST.items():
    _pramanik = re.sub(r"\b" + _k + r"\b", _v, _pramanik)

INSTANCES = [
    TestFunction("Trigexp1",
                 "-8 + sin(x_2 - x_3)*sin(x_2 + x_3) - exp(x_1 - x_2)*x_1 - 4*x_2 + 3*x_2^3 + 2*x_3",
                 description="Trigexp1-020.bch, constraint 2 (variables x1, x2, x3)"),
    TestFunction("Trigexp2 (1)",
                 "-5 + sin(x_1 - x_2 - x_3)*sin(x_1 + x_2 - x_3) + 2*x_2 + 3*(x_1 - x_3)^3",
                 description="Trigexp2-3.bch, constraint 1"),
    TestFunction("Trigexp2 (2)",
                 "-3 + 4*x_2 - exp(x_1 - x_2 - x_3)*(x_1 - x_3)",
                 description="Trigexp2-3.bch, constraint 2"),
    TestFunction("Pramanik", _pramanik, domain=(0.06, 1.0), sizes=("small", "medium"),
                 description="Pramanik.bch, constraint 1 (f0, f1, y0, y1 replaced by their values)"),
]

SIZES = {"small": (0.01, 0.1), "medium": (0.1, 1.0), "large": (1.0, 10.0)}

ALL = CASES + INSTANCES
