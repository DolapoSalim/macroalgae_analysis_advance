"""Step 11 — Work the biology numbers back out from the messy survey data.

We pretend we don't know the 13 biology parameters. Each fit starts from a
different rough first guess (each number off by up to 2x, like values borrowed
from the literature) and lets the computer search for the numbers that make the
model match the recorded data best. Then we open the answer key.

Run from the project root:
    python -m steps.step11_fit_parameters
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.estimation import FITTED, build_segments, fit, predict
from src.parameters import BASELINE

N_STARTS = 3
STORM_YEAR = 8
OUT = Path("figures")


def main():
    observed = pd.read_csv("data/synthetic_observed/observations.csv")
    clean = pd.read_csv("data/synthetic_clean/true_states.csv")
    meta = json.loads(Path("data/metadata/sites.json").read_text())
    segments = build_segments(observed, meta["sites"], STORM_YEAR)
    truth = np.array([BASELINE[k] for k in FITTED])
    n_values = sum(int((~np.isnan(s["obs"])).sum()) for s in segments)
    print(f"Fitting {len(FITTED)} biology numbers + {3 * len(segments)} start values "
          f"to {n_values} recorded values...\n")

    rng = np.random.default_rng(11)
    fits = []
    for k in range(N_STARTS):
        guess = truth * np.exp(rng.uniform(-np.log(2), np.log(2), len(truth)))
        theta, starts, result = fit(segments, guess)
        fits.append({"guess": guess, "theta": theta, "starts": starts, "misfit": 2 * result.cost})
        print(f"  fit {k + 1}: misfit = {2 * result.cost:7.1f}   ({result.nfev} model runs)")
    print(f"  (a perfect model would score about {n_values} - one unit per recorded value)\n")

    print(f"{'parameter':8s} {'TRUE':>7s}" + "".join(f"{'fit ' + str(k + 1):>9s}" for k in range(N_STARTS))
          + f"{'off by':>9s}{'fits agree?':>13s}")
    for i, name in enumerate(FITTED):
        vals = np.array([f["theta"][i] for f in fits])
        off = best_off = vals[np.argmin([f["misfit"] for f in fits])] / truth[i] - 1
        agree = "yes" if vals.max() / vals.min() < 1.02 else "NO"
        print(f"{name:8s} {truth[i]:7.3f}" + "".join(f"{v:9.3f}" for v in vals)
              + f"{off:+9.0%}{agree:>13s}")

    Path("data/fits").mkdir(parents=True, exist_ok=True)
    Path("data/fits/step11_fits.json").write_text(json.dumps({
        "parameters": FITTED, "truth": truth.tolist(),
        "fits": [{"guess": f["guess"].tolist(), "theta": f["theta"].tolist(),
                  "starts": f["starts"].tolist(), "misfit": f["misfit"]} for f in fits]}, indent=2))

    # ---- Figure 1: best fit vs data, per site
    best = min(fits, key=lambda f: f["misfit"])
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), sharex=True, sharey=True)
    colours = ["#1f77b4", "#d62728", "#2ca02c"]
    names = ["coral", "seaweed", "herbivores"]
    sites = list(meta["sites"])
    for ax, site in zip(axes.ravel(), sites):
        d = observed[observed.site == site]
        c = clean[clean.site == site]
        for seg, x0 in zip(segments, best["starts"]):
            if seg["site"] != site:
                continue
            t_end = STORM_YEAR if seg["t_start"] == 0 else 20
            tt = np.linspace(seg["t_start"], t_end - 1e-6, 200)
            pred = predict(best["theta"], x0, seg["settings"], tt, seg["t_start"])
            for j in range(3):
                ax.plot(tt, pred[:, j], color=colours[j], lw=2.2)
        for j, var in enumerate(names):
            ax.plot(c.year, c[var], color=colours[j], lw=1, ls=":", alpha=0.9)
            ax.plot(d.year, d[var], "o", ms=2.5, color=colours[j], alpha=0.6)
        ax.set_title(site, fontsize=10)
        ax.set_ylim(0, 1.2)
        ax.grid(alpha=0.3)
    for j, var in enumerate(names):
        axes[0, 0].plot([], [], color=colours[j], lw=2.2, label=f"fitted model: {var}")
    axes[0, 0].plot([], [], "k:", label="truth")
    axes[0, 0].plot([], [], "ko", ms=3, label="recorded data")
    axes[0, 0].legend(fontsize=7, loc="upper right")
    for ax in axes[1]:
        ax.set_xlabel("Year")
    fig.suptitle("Best fitted model (thick lines) vs recorded data (dots) and the hidden truth (dotted)")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "step11_fit_vs_data.png", dpi=120)

    # ---- Figure 2: how close is each fitted number to the truth?
    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(len(FITTED))
    ax.axhspan(0.8, 1.25, color="#2ca02c", alpha=0.12, label="within 25% of the truth")
    ax.axhline(1, color="black", lw=1)
    for k, (f, m) in enumerate(zip(fits, ["o", "s", "^"])):
        ax.plot(x + (k - 1) * 0.15, f["theta"] / truth, m, ms=7, label=f"fit {k + 1}")
        ax.plot(x + (k - 1) * 0.15, f["guess"] / truth, m, ms=4, mfc="none", color="grey",
                label="first guesses" if k == 0 else None)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(FITTED, rotation=30)
    ax.set_ylabel("fitted value / true value\n(1 = perfect)")
    ax.set_title("How close did each fit get to the true numbers?")
    ax.legend(fontsize=8, ncol=5, loc="upper center")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT / "step11_fitted_vs_true.png", dpi=120)
    print(f"\nResults saved in data/fits/, figures saved in {OUT}/")


if __name__ == "__main__":
    main()
