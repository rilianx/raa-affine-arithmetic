"""Tables 2 and 3 of the paper and the paired statistics, from results/*.json.

    python experiments/make_tables.py            # prints the tables and the statistics
    python experiments/make_tables.py --latex    # also writes results/table_<set>.tex

Gap(M) = (W_M - W_OPT) / W_OPT * 100, where W_M is the residual width of method M
(the vertical distance between its two parallel envelopes) and W_OPT the width of the
range of f(x) - L(x) over the box, L being the linear part of the ZCAA enclosure.
"""
import argparse
import json
import os

import numpy as np
from scipy import stats

from common import RESULTS

REL_TOL = 1e-9      # widths that differ by less than this relative amount count as equal


def load(name):
    path = os.path.join(RESULTS, f"{name}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def _z(v):
    return 0.0 if abs(v) < 0.005 else v


def summarize(records, with_runs=False, with_reduction=False):
    rows, groups = [], {}
    for r in records:
        groups.setdefault((r["function"], r["size"]), []).append(r)
    funcs = list(dict.fromkeys(r["function"] for r in records))
    per_function = {}
    for f in funcs:
        all_af, all_ra, wins, ties, losses = [], [], 0, 0, 0
        wins_tm = losses_tm = 0
        for size in ("small", "medium", "large"):
            sel = groups.get((f, size))
            if not sel:
                continue
            wa = np.array([r["af2_width"] for r in sel])
            wr = np.array([r["zcaa_width"] for r in sel])
            wt = np.array([r["tm_width"] for r in sel])
            wo = np.array([r["opt_width"] for r in sel])
            tol = REL_TOL * np.maximum(wa, wr)
            wr = np.where(np.abs(wr - wa) <= tol, wa, wr)     # equal up to rounding
            ga, gr = (wa - wo) / wo * 100, (wr - wo) / wo * 100
            gt = (wt - wo) / wo * 100
            d = ga - gr
            red = (wa - wr) / wa * 100
            red_tm = (wt - wr) / wt * 100
            rows.append(dict(function=f, size=size, runs=len(sel),
                             af2=(ga.mean(), ga.std()), zcaa=(gr.mean(), gr.std()), tm=(gt.mean(), gt.std()),
                             median=np.median(d), p10=np.percentile(d, 10), p90=np.percentile(d, 90),
                             reduction=red.mean(), reduction_tm=red_tm.mean()))
            tol_t = REL_TOL * np.maximum(wt, wr)
            wins_tm += int(np.sum(wr < wt - tol_t))
            losses_tm += int(np.sum(wr > wt + tol_t))
            wins += int(np.sum(wr < wa - tol))
            losses += int(np.sum(wr > wa + tol))
            ties += int(np.sum(np.abs(wr - wa) <= tol))
            all_af += list(ga)
            all_ra += list(gr)
        all_af, all_ra = np.array(all_af), np.array(all_ra)
        p = (stats.wilcoxon(all_af, all_ra, alternative="greater").pvalue
             if np.any(all_af != all_ra) else float("nan"))
        per_function[f] = dict(wins=wins, ties=ties, losses=losses, wilcoxon_p=p, gaps=(all_af, all_ra),
                               wins_tm=wins_tm, losses_tm=losses_tm)
    return rows, per_function


def latex(rows, with_runs, with_reduction):
    lines, last = [], None
    for r in rows:
        first = r["function"] != last
        if last is not None and first:
            lines.append(r"\hline")
        last = r["function"]
        if with_reduction:
            # compact table (Table 3): function name once, size as S/M/L, mean gaps only
            name = r["function"] if first else ""
            cells = [name, r["size"][0].upper(),
                     f"{r['tm'][0]:.2f}", f"{r['af2'][0]:.2f}", f"{r['zcaa'][0]:.2f}",
                     f"{_z(r['reduction_tm']):.1f}", f"{_z(r['reduction']):.1f}"]
        else:
            cells = [f"{r['function']} ({r['size']})"]
            if with_runs:
                cells.append(str(r["runs"]))
            cells += [f"${r['tm'][0]:.2f} \\pm {r['tm'][1]:.2f}$",
                      f"${r['af2'][0]:.2f} \\pm {r['af2'][1]:.2f}$", f"${r['zcaa'][0]:.2f} \\pm {r['zcaa'][1]:.2f}$",
                      f"${_z(r['median']):.2f}$ $[{_z(r['p10']):.2f}, {_z(r['p90']):.2f}]$"]
        lines.append(" & ".join(cells) + r" \\")
    lines[-1] = lines[-1][:-3]
    return "\n".join(lines)


def report(name, records, with_runs, with_reduction, write_latex):
    rows, per_f = summarize(records)
    print(f"\n=== {name} ===")
    print(f"{'Function':<22}{'Runs':>5}{'TM gap (%)':>20}{'AF2 gap (%)':>20}{'ZCAA gap (%)':>20}"
          f"{'AF2-ZCAA median [P10, P90]':>30}{'Red./TM':>9}{'Red./AF2':>9}")
    for r in rows:
        print(f"{r['function'] + ' (' + r['size'] + ')':<22}{r['runs']:>5}"
              f"{r['tm'][0]:>11.2f} ± {r['tm'][1]:<6.2f}"
              f"{r['af2'][0]:>11.2f} ± {r['af2'][1]:<6.2f}{r['zcaa'][0]:>11.2f} ± {r['zcaa'][1]:<6.2f}"
              f"{_z(r['median']):>12.2f} [{_z(r['p10']):.2f}, {_z(r['p90']):.2f}]{'':>3}"
              f"{_z(r['reduction_tm']):>8.1f}{_z(r['reduction']):>9.1f}")
    print(f"\n{'Function':<16}{'ZCAA tighter':>12}{'equal':>8}{'AF2 tighter':>13}{'Wilcoxon p':>14}"
          f"{'ZCAA<TM':>9}{'ZCAA>TM':>9}")
    tot = [0, 0, 0]
    all_af, all_ra = [], []
    for f, s in per_f.items():
        print(f"{f:<16}{s['wins']:>12}{s['ties']:>8}{s['losses']:>13}{s['wilcoxon_p']:>14.1e}"
              f"{s['wins_tm']:>9}{s['losses_tm']:>9}")
        tot[0] += s["wins"]; tot[1] += s["ties"]; tot[2] += s["losses"]
        all_af += list(s["gaps"][0]); all_ra += list(s["gaps"][1])
    all_af, all_ra = np.array(all_af), np.array(all_ra)
    p = stats.wilcoxon(all_af, all_ra, alternative="greater").pvalue
    d = all_af - all_ra
    better = d[d > 1e-7]
    print(f"{'Total':<16}{tot[0]:>12}{tot[1]:>8}{tot[2]:>13}{p:>14.1e}")
    if len(better):
        print(f"When ZCAA is tighter, its gap is lower by {better.mean():.1f} points on average "
              f"and by up to {better.max():.0f} points.")
    nontrivial = [r for r in records if r["function"] != "Case 1"]
    wz = np.array([r["zcaa_width"] for r in nontrivial])
    for label, key in (("AF2", "af2_width"), ("TM", "tm_width")):
        w = np.array([r[key] for r in nontrivial])
        red = (w - wz) / w * 100
        print(f"Width reduction of ZCAA vs {label} (excluding Case 1): mean {red.mean():.1f}%, "
              f"median {np.median(red):.1f}%, P90 {np.percentile(red, 90):.1f}%, max {red.max():.1f}%")
    if write_latex:
        out = os.path.join(RESULTS, f"table_{name}.tex")
        with open(out, "w") as fh:
            fh.write(latex(rows, with_runs, with_reduction) + "\n")
        print("written", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--latex", action="store_true")
    args = ap.parse_args()
    for name, runs, red in (("controlled", False, False), ("instances", True, True)):
        recs = load(name)
        if recs is None:
            print(f"results/{name}.json not found: run experiments/run_experiments.py {name}")
            continue
        report(name, recs, runs, red, args.latex)
