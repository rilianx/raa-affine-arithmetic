"""Running time of one evaluation of each test expression with ZCAA, AF2 and TM.

    python experiments/timing.py     # writes results/timing.json

For every expression, 30 random boxes of each size are drawn and the expression is
evaluated 30 times per box with each arithmetic. Times are in milliseconds and
depend on the machine; the ratio ZCAA/AF2 is the relevant quantity.
"""
import json
import os
import random
import time

from common import RESULTS, random_box
from zcaa import set_box, zcaa_variables, af2_variables, tm_variables
from benchmarks.functions import ALL, SIZES

BOXES, REPS = 30, 30

if __name__ == "__main__":
    random.seed(11)
    out = {}
    for f in ALL:
        t_zcaa = t_af2 = t_tm = 0.0
        n = 0
        for size in ("small", "medium", "large"):
            for _ in range(BOXES):
                box = random_box(f, size, SIZES)
                set_box(box)
                xr, xa, xt = zcaa_variables(3), af2_variables(box), tm_variables(3)
                t0 = time.perf_counter()
                for _ in range(REPS):
                    f(*xr)
                t1 = time.perf_counter()
                for _ in range(REPS):
                    f(*xa)
                t2 = time.perf_counter()
                for _ in range(REPS):
                    f(*xt)
                t3 = time.perf_counter()
                t_zcaa += t1 - t0
                t_af2 += t2 - t1
                t_tm += t3 - t2
                n += REPS
        out[f.name] = {"zcaa_ms": t_zcaa / n * 1e3, "af2_ms": t_af2 / n * 1e3, "tm_ms": t_tm / n * 1e3,
                       "zcaa/af2": t_zcaa / t_af2, "zcaa/tm": t_zcaa / t_tm}
        o = out[f.name]
        print(f"{f.name:<14} ZCAA {o['zcaa_ms']:.4f} ms   AF2 {o['af2_ms']:.4f} ms   TM {o['tm_ms']:.4f} ms   "
              f"ZCAA/AF2 {o['zcaa/af2']:.2f}   ZCAA/TM {o['zcaa/tm']:.2f}", flush=True)
    with open(os.path.join(RESULTS, "timing.json"), "w") as fh:
        json.dump(out, fh, indent=1)
