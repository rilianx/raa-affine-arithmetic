"""Figure 2 of the paper: gap of ZCAA against gap of AF2 for every run of results/controlled.json.

    python experiments/make_figure.py      # writes results/ZCAAvsAF2.png
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402
import numpy as np                                                     # noqa: E402
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset  # noqa: E402

from common import RESULTS                                             # noqa: E402

with open(os.path.join(RESULTS, "controlled.json")) as fh:
    recs = json.load(fh)
wo = np.array([r["opt_width"] for r in recs])
ga = (np.array([r["af2_width"] for r in recs]) - wo) / wo
gr = (np.array([r["zcaa_width"] for r in recs]) - wo) / wo

plt.rcParams.update({"font.size": 12})
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(ga, gr, s=10, alpha=0.35, color="#2a9d8f")
m = max(ga.max(), gr.max()) * 1.02
ax.plot([0, m], [0, m], "k--", lw=1)
ax.set_xlim(0, m); ax.set_ylim(0, m)
ax.set_xlabel("AF2 Gap"); ax.set_ylabel("ZCAA Gap"); ax.set_title("ZCAA vs AF2"); ax.grid(alpha=0.3)
ins = inset_axes(ax, width="42%", height="42%", loc="lower right", borderpad=2)
ins.scatter(ga, gr, s=8, alpha=0.35, color="#2a9d8f")
ins.plot([0.2, 0.8], [0.2, 0.8], "k--", lw=1)
ins.set_xlim(0.2, 0.8); ins.set_ylim(0.2, 0.8); ins.tick_params(labelsize=9)
mark_inset(ax, ins, loc1=2, loc2=4, fc="none", ec="0.6")
out = os.path.join(RESULTS, "ZCAAvsAF2.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
print("written", out)
