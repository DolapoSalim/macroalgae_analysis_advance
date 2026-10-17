"""Step 7 — Early-warning signals: can we see a collapse coming?

Story: a healthy reef is slowly overfished. Grazing falls from 0.40 to 0.10
over 300 years. The reef is surveyed once a year (coral cover only), with
small everyday random kicks (noise 0.01).

For every reef that collapses (coral falls below 50%) we look ONLY at the
survey data from before the collapse (minus 5 years) and ask:
were the wobbles getting wider (variance) and slower (autocorrelation)?

Fair test: we run the same analysis on reefs that are NOT being overfished
(grazing stays 0.40). Any 'warning' there is a false alarm.

Run from the project root:
    python -m steps.step07_early_warnings
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.parameters import make_parameters
from src.stochastic import simulate_noisy
from src.early_warning import rolling_indicators, trend_strength

YEARS = 300
G_START, G_END = 0.40, 0.10
SIGMA = 0.01
N_REEFS = 200
COLLAPSE_CORAL = 0.5     # "collapsed" = coral below 50%
STOP_BEFORE = 5          # analyse data up to 5 years before the collapse
WARNING_TAU = 0.3        # a clear rising trend counts as a warning
CORAL_HOME = (0.809, 0.0, 0.467)
OUT = Path("figures")


def grazing_at(t):
    return G_START + (G_END - G_START) * t / YEARS


def main():
    rng = np.random.default_rng(7)
    t, fished = simulate_noisy(CORAL_HOME, lambda t: make_parameters(G_max=grazing_at(t)),
                               YEARS, SIGMA, N_REEFS, rng)
    _, steady = simulate_noisy(CORAL_HOME, make_parameters(G_max=G_START),
                               YEARS, SIGMA, N_REEFS, rng)

    results = {"overfished (collapsing)": [], "not fished (steady)": []}
    collapse_grazing = []
    example = None
    for r in range(N_REEFS):
        below = np.where(fished[:, 0, r] < COLLAPSE_CORAL)[0]
        if len(below) == 0:
            continue
        k = below[0]
        collapse_grazing.append(grazing_at(t[k]))
        end = k - STOP_BEFORE
        for name, series in [("overfished (collapsing)", fished[:end, 0, r]),
                             ("not fished (steady)", steady[:end, 0, r])]:
            var, ac = rolling_indicators(series, window=len(series) // 2)
            results[name].append((trend_strength(var), trend_strength(ac)))
        if example is None:
            example = (r, k)

    print(f"{len(collapse_grazing)} of {N_REEFS} overfished reefs collapsed "
          f"(coral < {COLLAPSE_CORAL:.0%}).")
    print(f"Grazing when they collapsed: typically {np.median(collapse_grazing):.2f} "
          f"(range {np.percentile(collapse_grazing, 10):.2f}-{np.percentile(collapse_grazing, 90):.2f})")
    print(f"\nWarning = BOTH wobble width and slowness clearly rising (trend > {WARNING_TAU}):")
    share = {}
    for name, v in results.items():
        v = np.array(v)
        share[name] = np.mean((v[:, 0] > WARNING_TAU) & (v[:, 1] > WARNING_TAU))
        print(f"  {name:25s} warning in {share[name]:4.0%} of reefs   "
              f"(typical trend: width {np.median(v[:, 0]):+.2f}, slowness {np.median(v[:, 1]):+.2f})")

    # ---- Figure 1: one example reef
    r, k = example
    end = k - STOP_BEFORE
    series = fished[:end, 0, r]
    var, ac = rolling_indicators(series, window=len(series) // 2)
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    axes[0].plot(t, fished[:, 0, r], color="#1f77b4")
    tipping_year = (G_START - 0.18) / (G_START - G_END) * YEARS
    axes[0].axvline(tipping_year, color="black", ls=":")
    axes[0].text(tipping_year - 3, 0.3, "grazing passes the\ntipping point (0.18)",
                 fontsize=8, ha="right")
    axes[0].axvline(t[k], color="#d62728", ls="--")
    axes[0].text(t[k] + 2, 0.75, "collapse", color="#d62728")
    axes[0].axvspan(0, t[end - 1], color="grey", alpha=0.12)
    axes[0].text(3, 0.12, "data used for the warning (before the collapse)", fontsize=8)
    axes[0].set_ylabel("Coral cover")
    axes[0].set_ylim(0, 1)
    axes[1].plot(t[:end], var, color="#ff7f0e")
    axes[1].set_ylabel("Wobble width\n(variance)")
    axes[2].plot(t[:end], ac, color="#9467bd")
    axes[2].set_ylabel("Slowness\n(autocorrelation)")
    axes[2].set_xlabel("Years   (grazing falls from 0.40 at year 0 to 0.10 at year 300)")
    for ax in axes:
        ax.grid(alpha=0.3)
    axes[0].set_title("One overfished reef: the warning signs before collapse")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "step07_example_warning.png", dpi=120)

    # ---- Figure 2: how reliable is the warning?
    fig, ax = plt.subplots(figsize=(6, 4))
    names = list(share)
    ax.bar(["Overfished reefs\n(really collapsing)", "Steady reefs\n(not collapsing)"],
           [share[n] * 100 for n in names], color=["#d62728", "#7f7f7f"])
    for i, n in enumerate(names):
        ax.text(i, share[n] * 100 + 2, f"{share[n]:.0%}", ha="center", fontsize=12)
    ax.set_ylabel("% of reefs showing a warning")
    ax.set_ylim(0, 100)
    ax.set_title("Warning signs: real alarms vs false alarms")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "step07_warning_reliability.png", dpi=120)
    print(f"\nFigures saved in {OUT}/")


if __name__ == "__main__":
    main()
