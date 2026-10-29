"""Running time of one evaluation of each test expression with RAA and AF2.

    python experiments/timing.py     # writes results/timing.json

For every expression, 30 random boxes of each size are drawn and the expression is
evaluated 30 times per box with each arithmetic. Times are in milliseconds and
depend on the machine; the ratio RAA/AF2 is the relevant quantity.
"""
import json
import os
import random
import time

from common import RESULTS, random_box
from raa import set_box, raa_variables, af2_variables
from benchmarks.functions import ALL, SIZES

BOXES, REPS = 30, 30

if __name__ == "__main__":
    random.seed(11)
    out = {}
    for f in ALL:
        t_raa = t_af2 = 0.0
        n = 0
        for size in ("small", "medium", "large"):
            for _ in range(BOXES):
                box = random_box(f, size, SIZES)
                set_box(box)
                xr, xa = raa_variables(3), af2_variables(box)
                t0 = time.perf_counter()
                for _ in range(REPS):
                    f(*xr)
                t1 = time.perf_counter()
                for _ in range(REPS):
                    f(*xa)
                t2 = time.perf_counter()
                t_raa += t1 - t0
                t_af2 += t2 - t1
                n += REPS
        out[f.name] = {"raa_ms": t_raa / n * 1e3, "af2_ms": t_af2 / n * 1e3, "ratio": t_raa / t_af2}
        print(f"{f.name:<14} RAA {out[f.name]['raa_ms']:.4f} ms   AF2 {out[f.name]['af2_ms']:.4f} ms   "
              f"ratio {out[f.name]['ratio']:.2f}", flush=True)
    with open(os.path.join(RESULTS, "timing.json"), "w") as fh:
        json.dump(out, fh, indent=1)
