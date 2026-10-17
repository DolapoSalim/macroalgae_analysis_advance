"""Step 12 — Which fitted numbers can we trust?

Part A  Profile likelihood: lock each number at 1/3x ... 3x its fitted value,
        refit everything else, and see how much worse the fit gets.
Part B  Fisher information: which COMBINATIONS of numbers the data pins down
        ('stiff') and which it hardly constrains ('sloppy').
Part C  3D pictures of the misfit 'landscape' for two pairs of numbers.

Run from the project root (takes ~10 minutes):
    python -m steps.step12_identifiability
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.estimation import FITTED, build_segments, residuals
from src.identifiability import biology_information, pack, profile, residual_jacobian

FACTORS = np.array([1 / 3, 1 / 2, 1 / 1.5, 1 / 1.2, 1.0, 1.2, 1.5, 2.0, 3.0])
STORM_YEAR = 8
OUT = Path("figures")


def crossing(factors, rise, threshold, side):
    """Where the misfit rise first exceeds the threshold, going left (-1) or right (+1)."""
    centre = int(np.argmin(np.abs(np.log(factors))))
    idx = range(centre, len(factors)) if side > 0 else range(centre, -1, -1)
    prev = centre
    for k in idx:
        if rise[k] > threshold:
            x0, x1 = np.log(factors[prev]), np.log(factors[k])
            y0, y1 = rise[prev], rise[k]
            return np.exp(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))
        prev = k
    return None


def main():
    observed = pd.read_csv("data/synthetic_observed/observations.csv")
    meta = json.loads(Path("data/metadata/sites.json").read_text())
    segments = build_segments(observed, meta["sites"], STORM_YEAR)
    saved = json.loads(Path("data/fits/step11_fits.json").read_text())
    truth = np.array(saved["truth"])
    best = min(saved["fits"], key=lambda f: f["misfit"])
    z_best = pack(np.array(best["theta"]), best["starts"])
    r = residuals(z_best, segments)
    best_misfit = float(np.sum(r ** 2))
    scale = best_misfit / (len(r) - len(z_best))      # misfit per degree of freedom
    threshold = 3.84 * scale                         # 95% confidence, adjusted for model imperfection
    print(f"Best misfit {best_misfit:.1f}. A number is 'pinned down' where moving it raises "
          f"the misfit by more than {threshold:.1f}.\n")

    # ---------------- Part A: profiles
    profiles = {}
    for i, name in enumerate(FITTED):
        profiles[name] = profile(segments, z_best, i, FACTORS) - best_misfit
        print(f"  profiled {name}")
    rows = []
    for i, name in enumerate(FITTED):
        rise = profiles[name]
        lo = crossing(FACTORS, rise, threshold, -1)
        hi = crossing(FACTORS, rise, threshold, +1)
        verdict = ("PINNED DOWN" if lo and hi else "ONE-SIDED" if (lo or hi) else "NOT PINNED DOWN")
        fitted = best["theta"][i]
        rng_txt = (f"{fitted * lo if lo else 0:.3f} to " + (f"{fitted * hi:.3f}" if hi else ">" + f"{fitted * 3:.2f}")) \
            if (lo or hi) else "anything in 1/3x-3x"
        contains_truth = (lo is None or truth[i] >= fitted * lo) and (hi is None or truth[i] <= fitted * hi)
        rows.append((name, truth[i], fitted, rng_txt, verdict, contains_truth))
    print(f"\n{'number':8s} {'TRUE':>7s} {'fitted':>8s}   {'95% plausible range':24s} {'verdict':16s} truth inside?")
    for name, t, f, rng_txt, verdict, inside in rows:
        print(f"{name:8s} {t:7.3f} {f:8.3f}   {rng_txt:24s} {verdict:16s} {'yes' if inside else 'NO'}")

    fig, axes = plt.subplots(3, 5, figsize=(16, 9), sharey=True)
    colours = {"PINNED DOWN": "#2ca02c", "ONE-SIDED": "#ff7f0e", "NOT PINNED DOWN": "#d62728"}
    for ax, (i, name), row in zip(axes.ravel(), enumerate(FITTED), rows):
        values = best["theta"][i] * FACTORS / truth[i]
        ax.plot(values, profiles[name], "o-", color=colours[row[4]])
        ax.axhline(threshold, color="grey", ls="--", lw=1)
        ax.axvline(1, color="black", lw=1)
        ax.set_xscale("log")
        ax.set_ylim(-1, 60)
        ax.set_title(f"{name}: {row[4].lower()}", color=colours[row[4]], fontsize=10)
        ax.grid(alpha=0.3, which="both")
    for ax in axes.ravel()[len(FITTED):]:
        ax.axis("off")
    axes.ravel()[-1].text(0, 0.5, "x-axis: locked value / TRUE value\n(black line = truth)\n\n"
                          "y-axis: how much worse the best fit gets\n(dashed = 'pinned down' threshold)\n\n"
                          "Steep U  = data pins it down\nFlat line = data does not know it",
                          fontsize=10, va="center")
    for ax in axes[:, 0]:
        ax.set_ylabel("misfit increase")
    fig.suptitle("Profile likelihood: lock one number, refit the rest, watch the fit get worse (or not)")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "step12_profiles.png", dpi=120)

    # ---------------- Part B: Fisher information
    J = residual_jacobian(segments, z_best)
    info = biology_information(J) / scale
    w, v = np.linalg.eigh(info)
    print("\nFisher information: directions in 'number space', from sloppiest to stiffest")
    print("(std. error = how uncertain that combination is, on a log scale; >1 means 'unknown')")
    for k in range(len(w)):
        top = np.argsort(-np.abs(v[:, k]))[:3]
        combo = " ".join(f"{'+' if v[j, k] > 0 else '-'}{abs(v[j, k]):.2f}*log({FITTED[j]})" for j in top)
        se = 1 / np.sqrt(w[k]) if w[k] > 1e-9 else np.inf
        print(f"  {k + 1:2d}. std. error {se:8.3f}   mainly: {combo}")

    fig, ax = plt.subplots(figsize=(10, 4.5))
    shown = np.clip(w, 1e-2, None)
    labels = []
    for k in range(len(w)):
        top = np.argsort(-np.abs(v[:, k]))[:2]
        labels.append(" & ".join(FITTED[j] for j in top))
    ax.bar(range(len(w)), shown, color=["#d62728" if x < 1 else "#ff7f0e" if x < 100 else "#2ca02c" for x in w])
    if w[0] < 1e-6:
        ax.text(0, 0.02, "exactly\n0", ha="center", va="bottom", fontsize=8)
    ax.set_yscale("log")
    ax.set_xticks(range(len(w)))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("information (higher = better known)")
    ax.set_title("Fisher information: red = sloppy (data can't tell), green = stiff (well pinned down)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "step12_fisher_information.png", dpi=120)

    # ---------------- Part C: 3D misfit landscapes (other numbers held at the best fit)
    def landscape(i, j, keep_ratio_with=None, n=40, span=2.0):
        a = np.exp(np.linspace(-np.log(span), np.log(span), n))
        Z = np.zeros((n, n))
        for p_, fa in enumerate(a):
            for q_, fb in enumerate(a):
                z = z_best.copy()
                z[i] += np.log(fa)
                z[j] += np.log(fb)
                if keep_ratio_with is not None:      # keep r_H / K_H fixed
                    z[keep_ratio_with] += np.log(fa)
                Z[q_, p_] = np.sum(residuals(z, segments) ** 2) - best_misfit
        return best["theta"][i] * a / truth[i], best["theta"][j] * a / truth[j], Z

    iC, iD = FITTED.index("r_C"), FITTED.index("mu_C")
    iH, iM, iK = FITTED.index("r_H"), FITTED.index("mu_H"), FITTED.index("K_H")
    print("\nComputing 3D landscapes...")
    land1 = landscape(iC, iD)
    land2 = landscape(iH, iM, keep_ratio_with=iK)

    fig3 = make_subplots(rows=1, cols=2, specs=[[{"type": "surface"}, {"type": "surface"}]],
                         subplot_titles=["coral growth vs coral death: a long tilted valley",
                                         "fish birth vs fish death: a perfectly flat valley floor"])
    for col, (x, y, Z), (xl, yl) in [(1, land1, ("r_C / true", "mu_C / true")),
                                     (2, land2, ("r_H / true", "mu_H / true"))]:
        fig3.add_trace(go.Surface(x=x, y=y, z=np.minimum(Z, 300), colorscale="Viridis",
                                  showscale=False), row=1, col=col)
        fig3.update_scenes(dict(xaxis_title=xl, yaxis_title=yl, zaxis_title="misfit increase",
                                xaxis_type="log", yaxis_type="log"), row=1, col=col)
    fig3.update_layout(title="Misfit landscapes (lower = fits the data better) - drag to rotate")
    fig3.write_html(OUT / "step12_landscapes_3d.html", include_plotlyjs="cdn")

    fig = plt.figure(figsize=(13, 5.5))
    for k, ((x, y, Z), title, (xl, yl)) in enumerate([
            (land1, "Coral growth vs coral death\na long, tilted valley", ("r_C / true", "mu_C / true")),
            (land2, "Fish birth vs fish death\na perfectly flat valley floor", ("r_H / true", "mu_H / true"))]):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        X, Y = np.meshgrid(np.log10(x), np.log10(y))
        ax.plot_surface(X, Y, np.minimum(Z, 300), cmap="viridis", alpha=0.9, linewidth=0)
        ticks = [0.5, 1, 2]
        ax.set_xticks(np.log10(ticks))
        ax.set_xticklabels(ticks)
        ax.set_yticks(np.log10(ticks))
        ax.set_yticklabels(ticks)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_zlabel("misfit increase")
        ax.set_title(title)
        ax.view_init(elev=35, azim=-50)
    fig.tight_layout()
    fig.savefig(OUT / "step12_landscapes_3d.png", dpi=120)
    print(f"\nFigures saved in {OUT}/")


if __name__ == "__main__":
    main()
