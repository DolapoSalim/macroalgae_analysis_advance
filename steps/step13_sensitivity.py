"""Step 13 — Which numbers matter most for the reef's resilience?

Every biology number is allowed to be anywhere from half to double its baseline
value (like our Step 11 first guesses). We build ~30,000 versions of the reef
and measure two resilience outcomes for each:

  1. TIPPING POINT: the grazing level below which seaweed can invade a clean reef
  2. STORM TEST: coral cover 50 years after a strong storm hits the reef's own clean
     coral home (fishing 0, nutrients 0.5)

Grazing strength G_max is itself one of the varied numbers, so each version is
compared with its OWN grazing: safety margin = G_max - tipping point.

Uncertainty: how spread out are the outcomes?
Sensitivity: which numbers cause that spread? (Sobol indices)
Priorities:  combine with Step 12's trust ratings.

Run from the project root:
    python -m steps.step13_sensitivity
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from scipy.stats import qmc

from src.estimation import FITTED
from src.parameters import BASELINE, make_parameters
from src.sensitivity import coral_after_storm, coral_home, sobol_indices, tipping_grazing

N = 2048
SPREAD = 2.0          # each number anywhere from baseline/2 to baseline*2
OUT = Path("figures")

# Trust ratings from Step 12 (profile likelihood)
TRUST = {"r_C": "loose", "r_rec": "pinned", "mu_C": "loose", "a_CM": "pinned", "r_M": "pinned",
         "beta_N": "pinned", "mu_M": "pinned", "G_max": "pinned", "K_G": "pinned",
         "r_H": "half known", "K_H": "half known", "q_M": "half known", "mu_H": "unknown"}
TRUST_LEVEL = {"pinned": 0, "loose": 1, "half known": 2, "unknown": 3}


def to_params(U):
    """Rows of numbers in [0, 1] -> dict of parameter arrays (log-uniform 1/SPREAD..SPREAD)."""
    p = {k: np.full(len(U), v, dtype=float) for k, v in BASELINE.items()}
    for j, name in enumerate(FITTED):
        p[name] = BASELINE[name] * SPREAD ** (2 * U[:, j] - 1)
    return p


def main():
    k = len(FITTED)
    U = qmc.Sobol(d=2 * k, scramble=True, seed=13).random(N)
    A, B = U[:, :k], U[:, k:]
    AB = [A.copy() for _ in range(k)]
    for i in range(k):
        AB[i][:, i] = B[:, i]
    all_U = np.vstack([A, B] + AB)
    print(f"Evaluating {len(all_U)} versions of the reef...")
    p_all = to_params(all_U)

    outputs = {"Tipping point (grazing)": tipping_grazing(p_all),
               "Coral 50 yrs after storm": coral_after_storm(p_all)}

    # ---------------- uncertainty (use the A block: an even spread of versions)
    tip = outputs["Tipping point (grazing)"][:N]
    storm = outputs["Coral 50 yrs after storm"][:N]
    pA = to_params(A)
    margin = pA["G_max"] - tip
    valley = margin > 0
    print("\nUNCERTAINTY (if every number could be off by up to 2x):")
    print(f"  tipping point: middle value {np.median(tip):.2f}, 90% of versions between "
          f"{np.percentile(tip, 5):.2f} and {np.percentile(tip, 95):.2f}  (baseline 0.17)")
    print(f"  versions with NO coral valley (own grazing below own tipping point): {np.mean(~valley):.0%}")
    print(f"  versions where fish death >= fish birth (fish can't survive at all): "
          f"{np.mean(pA['mu_H'] >= pA['r_H']):.0%}")
    print(f"  strong storm on the clean coral home -> coral below 50% after 50 yrs:")
    print(f"      all versions {np.mean(storm < 0.5):.0%};  versions WITH a coral valley "
          f"{np.mean(storm[valley] < 0.5):.0%};  versions WITHOUT {np.mean(storm[~valley] < 0.5):.0%}")

    # What do the Step 11 fits (which disagreed on every fish number) say?
    saved = json.loads(Path("data/fits/step11_fits.json").read_text())
    print("\nTIPPING POINT FROM THE STEP 11 FITS (fish numbers differ, but...):")
    for j, fit in enumerate(saved["fits"], 1):
        p = make_parameters(**dict(zip(FITTED, fit["theta"])))
        print(f"  fit {j}: fish birth {p['r_H']:.3f}, death {p['mu_H']:.3f}, max {p['K_H']:.3f} "
              f"-> fish at coral home {coral_home(p)[1]:.3f}, tipping point {tipping_grazing(p):.3f}")
    print(f"  truth: tipping point {tipping_grazing(make_parameters()):.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].hist(np.clip(margin, -0.6, 0.8), bins=50, color="#9467bd", alpha=0.8)
    axes[0].axvline(0.40 - 0.172, color="black", ls="--", label="baseline (+0.23)")
    axes[0].axvline(0, color="#d62728", lw=2, label="0 = no safety margin left")
    axes[0].set_xlabel("Safety margin = own grazing - tipping point")
    axes[0].set_ylabel("number of reef versions")
    axes[0].set_title("Is there a coral valley?")
    axes[0].legend(fontsize=8)
    bins = np.linspace(0, 1, 51)
    axes[1].hist(storm[valley], bins=bins, color="#1f77b4", alpha=0.8, label="with a coral valley")
    axes[1].hist(storm[~valley], bins=bins, color="#d62728", alpha=0.6, label="without a coral valley")
    axes[1].axvline(0.5, color="black", ls="--", lw=1)
    axes[1].set_xlabel("Coral cover 50 years after a strong storm (from the clean coral home)")
    axes[1].set_title("Does the reef recover from a strong storm?")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle("Uncertainty: the same reef, if each biology number could be off by up to 2x")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "step13_uncertainty.png", dpi=120)

    # ---------------- sensitivity (Sobol)
    results = {}
    for name, f in outputs.items():
        f_A, f_B = f[:N], f[N:2 * N]
        f_AB = f[2 * N:].reshape(k, N)
        results[name] = sobol_indices(f_A, f_B, f_AB)
    print("\nSENSITIVITY (total Sobol index: share of the spread caused by each number):")
    print(f"{'number':8s} {'tipping point':>15s} {'storm test':>12s}   trust (Step 12)")
    order = np.argsort(-results["Tipping point (grazing)"][1])
    for i in order:
        print(f"{FITTED[i]:8s} {results['Tipping point (grazing)'][1][i]:15.2f} "
              f"{results['Coral 50 yrs after storm'][1][i]:12.2f}   {TRUST[FITTED[i]]}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (name, (s1, st, s1_ci, st_ci)) in zip(axes, results.items()):
        o = np.argsort(st)
        y = np.arange(k)
        ax.barh(y + 0.2, st[o], height=0.4, color="#1f77b4", label="total (alone + teamwork)",
                xerr=[st[o] - st_ci[0][o], st_ci[1][o] - st[o]], capsize=2)
        ax.barh(y - 0.2, np.clip(s1[o], 0, None), height=0.4, color="#aec7e8", label="alone")
        ax.set_yticks(y)
        ax.set_yticklabels([FITTED[i] for i in o])
        ax.set_xlabel("share of the outcome's spread")
        ax.set_title(name)
        ax.grid(axis="x", alpha=0.3)
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("Sensitivity: which numbers drive the outcome?")
    fig.tight_layout()
    fig.savefig(OUT / "step13_sensitivity.png", dpi=120)

    # ---------------- priorities: important AND poorly known?
    importance = np.maximum(results["Tipping point (grazing)"][1], results["Coral 50 yrs after storm"][1])
    trust = np.array([TRUST_LEVEL[TRUST[n]] for n in FITTED])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axvspan(0.1, 1, ymin=0.37, ymax=1, color="#d62728", alpha=0.1)
    ax.text(0.12, 3.35, "RESEARCH PRIORITY\n(matters a lot, poorly known)", color="#d62728", fontsize=9)
    jitter = np.random.default_rng(1).uniform(-0.15, 0.15, k)
    ax.scatter(importance, trust + jitter, s=60, color="#1f77b4")
    for i, n in enumerate(FITTED):
        ax.annotate(n, (importance[i], trust[i] + jitter[i]), xytext=(5, 3), textcoords="offset points", fontsize=9)
    ax.set_yticks(range(4))
    ax.set_yticklabels(["pinned down", "loosely known", "half known", "unknown"])
    ax.set_xscale("log")
    ax.set_xlim(0.002, 1)
    ax.set_ylim(-0.5, 3.8)
    ax.set_xlabel("importance (largest total Sobol index, log scale)")
    ax.set_ylabel("how well the data knows it (Step 12)")
    ax.set_title("Where should research effort go?")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "step13_priorities.png", dpi=120)

    # ---------------- 3D: tipping point vs the two most important numbers
    top2 = [FITTED[i] for i in np.argsort(-results["Tipping point (grazing)"][1])[:2]]
    f = np.linspace(0.5, 2, 40)
    F1, F2 = np.meshgrid(f, f)
    p = {kk: np.full(F1.size, v, dtype=float) for kk, v in BASELINE.items()}
    p[top2[0]] = BASELINE[top2[0]] * F1.ravel()
    p[top2[1]] = BASELINE[top2[1]] * F2.ravel()
    Z = tipping_grazing(p).reshape(F1.shape)
    fig3 = go.Figure(go.Surface(x=F1, y=F2, z=Z, colorscale="RdYlBu_r",
                                colorbar=dict(title="tipping<br>point")))
    fig3.add_trace(go.Surface(x=F1, y=F2, z=np.full_like(Z, 0.40), opacity=0.35, showscale=False,
                              colorscale=[[0, "grey"], [1, "grey"]], name="today's grazing 0.40"))
    fig3.update_layout(title=f"Tipping point vs {top2[0]} and {top2[1]} (grey sheet = today's grazing 0.40)",
                       scene=dict(xaxis_title=f"{top2[0]} (x baseline)", yaxis_title=f"{top2[1]} (x baseline)",
                                  zaxis_title="tipping point"))
    fig3.write_html(OUT / "step13_tipping_3d.html", include_plotlyjs="cdn")

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(F1, F2, Z, cmap="RdYlBu_r", alpha=0.9, linewidth=0)
    ax.plot_surface(F1, F2, np.full_like(Z, 0.40), color="grey", alpha=0.3)
    ax.set_xlabel(f"{top2[0]} (x baseline)")
    ax.set_ylabel(f"{top2[1]} (x baseline)")
    ax.set_zlabel("tipping point (grazing)")
    ax.set_title(f"Tipping point vs the two most important numbers\n"
                 f"grey sheet = today's grazing (0.40): above it = reef in danger")
    ax.view_init(elev=25, azim=-130)
    fig.tight_layout()
    fig.savefig(OUT / "step13_tipping_3d.png", dpi=120)
    print(f"\nFigures saved in {OUT}/")


if __name__ == "__main__":
    main()
