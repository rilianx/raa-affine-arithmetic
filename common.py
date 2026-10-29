"""Shared utilities for the experiments."""
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from zcaa import set_box, zcaa_variables, af2_variables, tm_variables   # noqa: E402

RESULTS = os.path.join(ROOT, "results")


def random_interval(min_val, max_val, min_width, max_width):
    """Random interval inside [min_val, max_val]; the width is clipped to the domain."""
    width = min(random.uniform(min_width, max_width), max_val - min_val)
    start = random.uniform(min_val, max_val - width)
    return (start, start + width)


def random_box(function, size, sizes):
    wmin, wmax = sizes[size]
    lo, hi = function.domain
    return [random_interval(lo, hi, wmin, wmax) for _ in range(3)]


def evaluate(function, box):
    """Enclosures of the function over the box with ZCAA, AF2 and first-order Taylor models."""
    set_box(box)
    r = function(*zcaa_variables(3))
    a = function(*af2_variables(box))
    t = function(*tm_variables(3))
    return r, a, t


def ibexopt_path():
    path = os.environ.get("IBEXOPT") or shutil.which("ibexopt")
    if not path:
        sys.exit("IbexOpt not found: install Ibex with the optim plugin and set IBEXOPT=/path/to/ibexopt")
    return path


def _ibex_best(expr, box, prec, timeout):
    """Best value found by IbexOpt for min expr over the box (upper bound of f*)."""
    dom = "(" + " ; ".join(f"[{a!r}, {b!r}]" for a, b in box) + ")"
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(f"variables\nx[{len(box)}] in {dom};\n\nminimize {expr};")
        name = fh.name
    try:
        out = subprocess.run([ibexopt_path(), "-a", str(prec), "-r", str(prec), "-t", str(timeout), name],
                             capture_output=True, text=True).stdout
    finally:
        for p in (name, name[:-4] + ".cov", name[:-4] + ".cov~"):
            if os.path.exists(p):
                os.remove(p)
    for line in out.splitlines():
        if "f* in" in line:
            return float(re.search(r"\[(.*),(.*)\]", line).group(2))
    raise RuntimeError("IbexOpt failed:\n" + out)


def optimal_residual_width(function, slopes, box, prec=1e-8, timeout=30):
    """Exact width of the range of f(x) - sum(slopes_i x_i) over the box (W_OPT)."""
    lin = "+".join(f"({s!r})*x[{i}]" for i, s in enumerate(slopes))
    g = f"{function.ibex_expr()}-({lin})"
    return -_ibex_best(f"-({g})", box, prec, timeout) - _ibex_best(g, box, prec, timeout)
