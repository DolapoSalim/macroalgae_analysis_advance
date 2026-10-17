"""Step 6 — Shake the ground: can random noise alone tip a healthy reef?

Every reef starts AT the coral home (the bottom of the coral valley).
No storm. Only small random kicks every time step, for 50 years.
Then we switch the noise off and see which home each reef settles in.

Run from the project root:
    python -m steps.step06_noise
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

from src.equilibria import find_equilibria
from src.parameters import make_parameters
from src.stability import classify
from src.stochastic import settle, simulate_noisy, which_home_many

GRAZING_VALUES = [0.185, 0.19, 0.195]        # all inside the two-homes zone
NOISE_LEVELS = [0.01, 0.02, 0.03, 0.04, 0.06, 0.08]
N_REPLICATES = 500
YEARS = 50
OUT = Path("figures")


def homes_and_pass(p):
    pts = [x for x in find_equilibria(p) if x[2] > 1e-6]
    homes = [x for x in pts if classify(x, p)[1].startswith("STABLE")]       # [coral, seaweed]
    saddle = [x for x in pts if classify(x, p)[1].startswith("UNSTABLE")][0]
    return homes, saddle


def main():
    rng = np.random.default_rng(6)
    OUT.mkdir(exist_ok=True)

    # ---- Experiment 1: tipping chance vs noise strength, for 3 grazing levels
    print(f"Share of {N_REPLICATES} healthy reefs that end in the SEAWEED home after "
          f"{YEARS} years of shaking:\n")
    print("noise   " + "   ".join(f"grazing {G}" for G in GRAZING_VALUES))
    table = {G: [] for G in GRAZING_VALUES}
    for sigma in NOISE_LEVELS:
        row = []
        for G in GRAZING_VALUES:
            p = make_parameters(G_max=G)
            homes, _ = homes_and_pass(p)
            _, path = simulate_noisy(homes[0], p, YEARS, sigma, N_REPLICATES, rng)
            ends = which_home_many(settle(path[-1], p), homes)
            share = np.mean(ends == 1)
            table[G].append(share)
            row.append(f"{share:13.0%}")
        print(f"{sigma:5.2f}  " + "".join(row))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for G, colour in zip(GRAZING_VALUES, ["#d62728", "#ff7f0e", "#2ca02c"]):
        ax.plot(NOISE_LEVELS, np.array(table[G]) * 100, "o-", color=colour, label=f"grazing {G}")
    ax.set_xlabel("Noise strength (how hard the ground shakes)")
    ax.set_ylabel(f"% of healthy reefs tipped\nto seaweed within {YEARS} years")
    ax.set_title("Can shaking alone tip a healthy reef?")
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "step06_tipping_vs_noise.png", dpi=120)

    # ---- Experiment 2: 12 example reefs, same start, same conditions, different fates
    G, sigma = 0.19, 0.04
    p = make_parameters(G_max=G)
    homes, saddle = homes_and_pass(p)
    t, path = simulate_noisy(homes[0], p, YEARS, sigma, 12, rng)
    fate = which_home_many(settle(path[-1], p), homes)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for r in range(path.shape[2]):
        ax.plot(t, path[:, 0, r], color="#1f77b4" if fate[r] == 0 else "#d62728", lw=1.2, alpha=0.8)
    ax.axhline(saddle[0], color="black", ls="--", lw=1)
    ax.text(1, saddle[0] + 0.01, "coral level of the mountain pass", fontsize=8)
    ax.plot([], [], color="#1f77b4", label=f"stayed healthy ({np.sum(fate == 0)})")
    ax.plot([], [], color="#d62728", label=f"tipped to seaweed ({np.sum(fate == 1)})")
    ax.set_xlabel("Years")
    ax.set_ylabel("Coral cover")
    ax.set_title(f"12 identical healthy reefs, shaken randomly (grazing {G}, noise {sigma})")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "step06_example_reefs.png", dpi=120)

    # interactive 3D: the same 12 reefs wandering through (coral, seaweed, herbivores)
    fig3 = go.Figure()
    for r in range(path.shape[2]):
        fig3.add_trace(go.Scatter3d(
            x=path[:, 0, r], y=path[:, 1, r], z=path[:, 2, r], mode="lines",
            line=dict(color="#1f77b4" if fate[r] == 0 else "#d62728", width=3),
            showlegend=False))
    for x, label, colour in [(homes[0], "coral home", "black"), (homes[1], "seaweed home", "black"),
                             (saddle, "mountain pass", "white")]:
        fig3.add_trace(go.Scatter3d(
            x=[x[0]], y=[x[1]], z=[x[2]], mode="markers+text", text=[label], name=label,
            marker=dict(size=7, color=colour, line=dict(color="black", width=3))))
    fig3.update_layout(
        title="12 shaken reefs: blue stayed healthy, red tipped to seaweed (drag to rotate)",
        scene=dict(xaxis_title="Coral (C)", yaxis_title="Seaweed (M)", zaxis_title="Herbivores (H)"))
    fig3.write_html(OUT / "step06_example_reefs_3d.html", include_plotlyjs="cdn")

    print(f"\nFigures saved in {OUT}/")


if __name__ == "__main__":
    main()
