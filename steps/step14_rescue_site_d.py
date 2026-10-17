"""Step 14 — Can we rescue Site D? (management experiments, guide §23)

Site D (fishing 0.04, nutrients 1.5) tipped to seaweed after the year-8 storm.
We start from its true state at the end of monitoring and try six plans.

  policy changes (permanent):   protect fish (fishing -> 0), cut pollution (nutrients -> 0.5)
  hands-on work (one-off):      replant coral (+20% cover), remove seaweed (-70%)

For each plan:
  - safety margin: today's grazing (0.40) minus the tipping point (positive = coral valley exists)
  - years until coral is back above 50% (no shaking, true biology)
  - chance of being a coral reef after 50 years (with everyday shaking, 300 repeats)
  - does the manager's FITTED model (Step 11) reach the same verdict?

Run from the project root:
    python -m steps.step14_rescue_site_d
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from src.estimation import FITTED
from src.parameters import make_parameters
from src.sensitivity import tipping_grazing
from src.stochastic import settle, simulate_noisy

SITE = "D Fished and polluted"
FISHING, NUTRIENTS = 0.04, 1.5
YEARS = 50
SIGMA = 0.01
N_REPEATS = 300
OUT = Path("figures")

PLANS = {
    #                       protect fish  cut pollution  replant coral  remove seaweed
    "Do nothing":           (False,       False,         False,         False),
    "Replant coral":        (False,       False,         True,          False),
    "Remove seaweed":       (False,       False,         False,         True),
    "Cut pollution":        (False,       True,          False,         False),
    "Protect fish":         (True,        False,         False,         False),
    "Everything together":  (True,        True,          True,          True),
}


def plan_setup(plan, x_now, biology=None):
    fish, pollution, replant, remove = PLANS[plan]
    p = make_parameters(**(biology or {}), F_H=0.0 if fish else FISHING, N=0.5 if pollution else NUTRIENTS)
    C, M, H = x_now
    if remove:
        M *= 0.3
    if replant:
        C = min(C + 0.20, 1 - M)
    return p, np.array([C, M, H])


def years_to_recover(t, coral, level=0.5):
    above = coral >= level
    if not above[-1]:
        return None
    below = np.where(~above)[0]
    return 0.0 if len(below) == 0 else t[below[-1] + 1]


def main():
    rng = np.random.default_rng(14)
    clean = pd.read_csv("data/synthetic_clean/true_states.csv")
    last = clean[clean.site == SITE].iloc[-1]
    x_now = np.array([last.coral, last.seaweed, last.herbivores])
    saved = json.loads(Path("data/fits/step11_fits.json").read_text())
    fitted = dict(zip(FITTED, min(saved["fits"], key=lambda f: f["misfit"])["theta"]))
    print(f"Site D today: coral {x_now[0]:.0%}, seaweed {x_now[1]:.0%}, fish {x_now[2]:.2f}\n")

    print(f"{'plan':20s} {'safety margin':>14s} {'coral >50% after':>17s} {'chance of coral reef':>21s} "
          f"{'fitted model says':>18s}")
    runs, chance, margin = {}, {}, {}
    for plan in PLANS:
        p, x0 = plan_setup(plan, x_now)
        margin[plan] = p["G_max"] - float(tipping_grazing(p))

        t, path = simulate_noisy(x0, p, 150, 0.0, 1, rng, record_every=0.25)
        runs[plan] = (t, path[:, :, 0])
        yrs = years_to_recover(t, path[:, 0, 0])

        _, noisy = simulate_noisy(x0, p, YEARS, SIGMA, N_REPEATS, rng)
        chance[plan] = np.mean(settle(noisy[-1], p)[0] > 0.5)

        pf, x0f = plan_setup(plan, x_now, biology=fitted)
        Xf = settle(x0f[:, None], pf, years=150)
        fitted_verdict = "recovers" if Xf[0, 0] > 0.5 else "stays seaweed"

        yrs_txt = f"{yrs:.0f} years" if yrs is not None else "never"
        print(f"{plan:20s} {margin[plan]:+14.2f} {yrs_txt:>17s} {chance[plan]:21.0%} {fitted_verdict:>18s}")

    # ---- Figure 1: coral over time under each plan
    fig, ax = plt.subplots(figsize=(9, 5))
    colours = ["#7f7f7f", "#8c564b", "#e377c2", "#ff7f0e", "#1f77b4", "#2ca02c"]
    for (plan, (t, path)), c in zip(runs.items(), colours):
        ax.plot(t, path[:, 0], color=c, lw=2.2, label=plan)
    ax.axhline(0.5, color="black", ls=":", lw=1)
    ax.set_xlabel("Years after the rescue starts")
    ax.set_ylabel("Coral cover")
    ax.set_title("Site D: coral cover under each rescue plan (no shaking)")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "step14_rescue_over_time.png", dpi=120)

    # ---- Figure 2: margin and chance side by side
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    names = list(PLANS)
    axes[0].barh(names, [margin[n] for n in names],
                 color=["#2ca02c" if margin[n] > 0 else "#d62728" for n in names])
    axes[0].axvline(0, color="black", lw=1)
    axes[0].set_xlabel("safety margin = grazing - tipping point\n(negative = no coral valley exists)")
    axes[0].set_title("Does the plan create a safe coral valley?")
    axes[1].barh(names, [chance[n] * 100 for n in names], color="#1f77b4")
    for i, n in enumerate(names):
        axes[1].text(chance[n] * 100 + 1, i, f"{chance[n]:.0%}", va="center")
    axes[1].set_xlim(0, 110)
    axes[1].set_xlabel(f"% of {N_REPEATS} shaken repeats that are coral reefs after {YEARS} years")
    axes[1].set_title("Chance of recovery (with everyday shaking)")
    for ax in axes:
        ax.grid(axis="x", alpha=0.3)
        ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUT / "step14_rescue_summary.png", dpi=120)

    # ---- 3D: path of each plan through (coral, seaweed, fish)
    fig3 = go.Figure()
    for (plan, (t, path)), c in zip(runs.items(), colours):
        fig3.add_trace(go.Scatter3d(x=path[:, 0], y=path[:, 1], z=path[:, 2], mode="lines",
                                    line=dict(width=6, color=c), name=plan))
    fig3.add_trace(go.Scatter3d(x=[x_now[0]], y=[x_now[1]], z=[x_now[2]], mode="markers+text",
                                text=["Site D today"], marker=dict(size=7, color="black"), name="Site D today"))
    fig3.update_layout(title="Site D: where each rescue plan takes the reef (drag to rotate)",
                       scene=dict(xaxis_title="Coral (C)", yaxis_title="Seaweed (M)", zaxis_title="Fish (H)"))
    fig3.write_html(OUT / "step14_rescue_3d.html", include_plotlyjs="cdn")
    print(f"\nFigures saved in {OUT}/")


if __name__ == "__main__":
    main()
