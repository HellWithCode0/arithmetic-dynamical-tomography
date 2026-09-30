"""
Figures for paper/ADT_foundations.md.

  fig1_single_probe_collisions.png   prime collisions for one probe vs X
                                     (recomputed from build/data/primes_c{c}.txt)
  fig2_collision_vs_depth.png        Pr[depth-K fingerprints of two primes agree]
  fig3_least_depth.png               least injective depth vs number of primes
                                     (both read results/depth_law.json, written by
                                     exp_depth_law.py from primes below 1e6)
"""
import json
import os
from collections import Counter
from math import log, sqrt
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from adt import SCRATCH, load_primes, zset_canon

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
os.makedirs(OUT, exist_ok=True)

# reference palette, slots 1-4 in fixed order (validated light mode on SURF;
# slots 3-4 are below 3:1 contrast, so every series also has a marker shape,
# a legend entry and a table in results/)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
MARKERS = ["o", "s", "^", "D"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10.5,
    "axes.titlesize": 12, "axes.titleweight": "semibold", "axes.titlecolor": INK,
    "legend.frameon": False, "legend.labelcolor": INK,
})


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(length=0)


with open(os.path.join(HERE, "..", "results", "depth_law.json")) as fh:
    LAW = json.load(fh)
INF = LAW["INF"]

# ---------------------------------------------------------------- fig 2
KS = LAW["KS"]
ANCHOR = {1: 64, 2: 32, 3: 16}          # reference line anchored in the asymptotic range
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for i, rec in enumerate(LAW["coll"][:3]):
    S, s = rec["S"], len(rec["S"])
    pts = [(K, f) for K, f in zip(KS, rec["frac"]) if f > 0]
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=SERIES[i], lw=2, marker=MARKERS[i], ms=6, mec=SURF, mew=1.5,
            label=f"|S| = {s}   (S = {{{', '.join(map(str, S))}}})")
    k0 = ANCHOR[s]
    y0 = dict(pts)[k0]
    k1 = min(256, xs[-1] * 2)
    ref = [y0 * (k / k0) ** (-2 * s) for k in (k0 / 2, k1)]
    ax.plot([k0 / 2, k1], ref, color=INK2, lw=1, ls=":", zorder=0)
    ax.annotate(f"slope −{2 * s}", (k1, ref[1]), textcoords="offset points", xytext=(6, 0),
                ha="left", va="center", color=INK2, fontsize=9)
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlim(1.8, 520)
style(ax, "Two primes' depth-K fingerprints agree w.p. ≍ K^(−2|S|)",
      "depth K", "fraction of prime pairs with equal fingerprint")
ax.legend(loc="lower left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2_collision_vs_depth.png"), dpi=160)

# ---------------------------------------------------------------- fig 3
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for i, s in enumerate((2, 3, 4, 5)):
    recs = [r for r in LAW["depth"] if len(r["S"]) == s]
    ns = recs[0]["n"]
    xs, ys = [], []
    for j, n in enumerate(ns):
        Kv = [r["K"][j] for r in recs if r["K"][j] < INF]
        if len(Kv) == len(recs):
            xs.append(n)
            ys.append(median(Kv))
    slopes = sorted(r["slope"] for r in recs if r["slope"] is not None)
    ax.plot(xs, ys, color=SERIES[i], lw=2, marker=MARKERS[i], ms=6, mec=SURF, mew=1.5,
            label=f"|S| = {s}: median of {len(recs)} sets, slope {median(slopes):.2f}")
    ref = [ys[-1] * (x / xs[-1]) ** (1 / s) for x in (xs[0], xs[-1])]
    ax.plot([xs[0], xs[-1]], ref, color=INK2, lw=1, ls=":", zorder=0)
    ax.annotate(f"1/{s}", (xs[-1], ref[1]), textcoords="offset points", xytext=(8, 0),
                ha="left", va="center", color=INK2, fontsize=9)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(100, 80000)
style(ax, "Depth needed to separate primes in (X/2, X]",
      "number of primes in the window", "least injective depth K")
ax.set_ylim(3.5, 1500)
ax.legend(loc="upper left", fontsize=9.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3_least_depth.png"), dpi=160)

# ---------------------------------------------------------------- fig 1
GEN = [1, 2, 3]
loc = load_primes(os.path.join(SCRATCH, "data", "primes_c{c}.txt"), GEN)
primes = sorted({p for (p, c) in loc})
Xs = [1000, 2000, 5000, 10000, 20000, 50000, 100000]
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for i, c in enumerate(GEN):
    ys = []
    for X in Xs:
        b = Counter(zset_canon(loc[(p, c)]) for p in primes if p <= X)
        ys.append(sum(v * (v - 1) // 2 for v in b.values()))
    ax.plot(Xs, ys, color=SERIES[i], lw=2, marker=MARKERS[i], ms=6, mec=SURF, mew=1.5, label=f"c = {c}")
h = [sqrt(X) / log(X) ** 2 for X in Xs]
scale = 120 / h[-1]
ax.plot(Xs, [scale * v for v in h], color=INK2, lw=1, ls=":", zorder=0)
ax.annotate("∝ X^(1/2) / log²X", (Xs[-1], scale * h[-1]), textcoords="offset points",
            xytext=(-6, -14), ha="right", color=INK2, fontsize=9)
ax.set_xscale("log")
ax.set_yscale("log")
style(ax, "One probe: prime collisions keep appearing",
      "X", "pairs p < q ≤ X with Per_c(p) ≅ Per_c(q)")
ax.legend(loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_single_probe_collisions.png"), dpi=160)
print("figures written to", os.path.abspath(OUT))
