"""Compute the residual widths of RAA and AF2 and the optimal width (IbexOpt).

    python experiments/run_experiments.py controlled   # Table 2 of the paper
    python experiments/run_experiments.py instances    # Table 3 of the paper

Results are written to results/<set>.json. The seeds, numbers of runs and IbexOpt
settings below are those used in the paper. Use --runs N for a quick check.
"""
import argparse
import json
import os
import random
import time

from common import RESULTS, random_box, evaluate, optimal_residual_width
from benchmarks.functions import CASES, INSTANCES, SIZES

# (function, seed, runs per size, IbexOpt precision, IbexOpt timeout in s)
CONFIG = {
    "controlled": {"seed": 2026, "functions": [(f, 100, 1e-8, 30) for f in CASES]},
    "instances": {"functions": [
        (INSTANCES[0], 2030, 100, 1e-8, 30),
        (INSTANCES[1], 2031, 50, 1e-6, 10),
        (INSTANCES[2], 2032, 50, 1e-6, 10),
        (INSTANCES[3], 2033, 50, 1e-6, 10),
    ]},
}


def run(name, max_runs=None):
    cfg = CONFIG[name]
    if "seed" in cfg:                      # one seed for the whole set
        random.seed(cfg["seed"])
        plan = [(f, None, n, p, t) for f, n, p, t in cfg["functions"]]
    else:                                  # one seed per function
        plan = cfg["functions"]
    records = []
    for f, seed, runs, prec, timeout in plan:
        if seed is not None:
            random.seed(seed)
        t0 = time.time()
        for size in f.sizes:
            for k in range(runs):
                box = random_box(f, size, SIZES)
                if max_runs is not None and k >= max_runs:
                    continue               # keep the random sequence of the full run
                r, a = evaluate(f, box)
                w_opt = optimal_residual_width(f, r.slopes, box, prec, timeout)
                records.append({"function": f.name, "size": size, "run": k, "box": box,
                                "raa_width": r.residual_width(), "af2_width": a.residual_width(),
                                "opt_width": w_opt})
        print(f"{f.name}: done in {time.time() - t0:.0f} s", flush=True)
    return records


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("set", choices=sorted(CONFIG))
    ap.add_argument("--runs", type=int, default=None, help="only the first RUNS boxes of each size")
    args = ap.parse_args()
    recs = run(args.set, args.runs)
    os.makedirs(RESULTS, exist_ok=True)
    suffix = "" if args.runs is None else f"_first{args.runs}"
    out = os.path.join(RESULTS, f"{args.set}{suffix}.json")
    with open(out, "w") as fh:
        json.dump(recs, fh, indent=0)
    print("written", out)
