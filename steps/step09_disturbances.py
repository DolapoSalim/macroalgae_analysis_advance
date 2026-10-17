"""Step 9 — Hit a healthy reef with real-world events. Does it recover?

Two reefs, both starting at their coral home:
  - well-grazed reef   (grazing 0.40, only one home)
  - degraded reef      (grazing 0.19, two homes: the danger zone)

Part A (no shaking): what happens after each event, and how long recovery takes.
Part B (with everyday shaking, degraded reef only): chance of recovery, from
         300 random repeats of each event.

Run from the project root:
    python -m steps.step09_disturbances
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

from src.disturbances import SCENARIOS, apply_shock, conditions_over_time
from src.equilibria import find_equilibria
from src.parameters import make_parameters
from src.stability import classify
from src.stochastic import settle, simulate_noisy, which_home_many

REEFS = {"well-grazed (0.40)": 0.40, "degraded (0.19)": 0.19}
YEARS = 60          # Part B
YEARS_A = 200       # Part A (slow recoveries need time)
SIGMA = 0.01
N_REPEATS = 300
OUT = Path("figures")


def homes_for(G):
    p = make_parameters(G_max=G)
    pts = [x for x in find_equilibria(p) if x[2] > 1e-6]
    homes = [x for x in pts if classify(x, p)[1].startswith("STABLE")]
    return p, homes  # homes[0] = coral home


def recovery_time(t, path, home, tol=0.02):
    """First year after which the reef stays within 2% of the coral home."""
    far = np.max(np.abs(path - home[None, :]), axis=1) > tol
    if far[-1]:
        return None
    last_far = np.where(far)[0]
    return 0.0 if len(last_far) == 0 else t[last_far[-1] + 1]


def main():
    rng = np.random.default_rng(9)
    OUT.mkdir(exist_ok=True)
    events = [s for s in SCENARIOS if s != "No event"]
    runs = {}

    # ---------------- Part A: no shaking
    print("PART A - no shaking. Outcome (recovery time) and highest seaweed cover:\n")
    print(f"{'event':18s}" + "".join(f"   {name:>33s}" for name in REEFS))
    for event in events:
        row = f"{event:18s}"
        for reef, G in REEFS.items():
            p_static, homes = homes_for(G)
            p_t = conditions_over_time(event, G_max=G)
            x0 = apply_shock(homes[0], event)
            t, path = simulate_noisy(x0, p_t, YEARS_A, 0.0, 1, rng, record_every=0.1)
            path = path[:, :, 0]
            runs[(event, reef)] = (t, path)
            end = which_home_many(settle(path[-1][:, None], p_static), homes)[0]
            rt = recovery_time(t, path, homes[0])
            outcome = (f"recovers in {rt:.0f} yrs" if rt is not None else "still recovering") \
                if end == 0 else "TIPS to seaweed"
            row += f"   {outcome:>17s}, max seaweed {path[:, 1].max():4.0%}"
        print(row)

    # ---------------- Part B: with shaking, degraded reef
    G = 0.19
    p_static, homes = homes_for(G)
    print(f"\nPART B - degraded reef (grazing {G}) with everyday shaking (noise {SIGMA}).")
    print(f"Chance of being back in the coral home {YEARS} years after the event "
          f"({N_REPEATS} repeats each):\n")
    chance = {}
    for event in SCENARIOS:
        p_t = conditions_over_time(event, G_max=G)
        x0 = apply_shock(homes[0], event)
        _, path = simulate_noisy(x0, p_t, YEARS, SIGMA, N_REPEATS, rng)
        chance[event] = np.mean(which_home_many(settle(path[-1], p_static), homes) == 0)
        print(f"  {event:18s} {chance[event]:5.0%}")

    # ---------------- Figure 1: coral over time for every event
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
    for ax, event in zip(axes.ravel(), events):
        for (reef, colour, ls) in [("well-grazed (0.40)", "#1f77b4", "-"), ("degraded (0.19)", "#d62728", "--")]:
            t, path = runs[(event, reef)]
            ax.plot(t, path[:, 0], color=colour, ls=ls, lw=2, label=f"coral, {reef}")
            ax.plot(t, path[:, 1], color=colour, ls=ls, lw=1, alpha=0.5, label=f"seaweed, {reef}")
        ax.set_title(event)
        ax.set_xlim(0, 120)
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("Years after the event")
    for ax in axes[:, 0]:
        ax.set_ylabel("Cover")
    axes[0, 0].legend(fontsize=7, loc="center right")
    fig.suptitle("After each event (no shaking): thick = coral, thin = seaweed")
    fig.tight_layout()
    fig.savefig(OUT / "step09_events_over_time.png", dpi=120)

    # ---------------- Figure 2: chance of recovery with shaking
    fig, ax = plt.subplots(figsize=(9, 4.5))
    names = list(chance)
    colours = ["#7f7f7f"] + ["#1f77b4"] * (len(names) - 1)
    ax.bar(names, [chance[n] * 100 for n in names], color=colours)
    for i, n in enumerate(names):
        ax.text(i, chance[n] * 100 + 2, f"{chance[n]:.0%}", ha="center")
    ax.set_ylabel(f"% back in coral home after {YEARS} years")
    ax.set_title(f"Degraded reef (grazing {G}) with everyday shaking: chance of recovery")
    ax.set_ylim(0, 100)
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "step09_recovery_chance.png", dpi=120)

    # ---------------- 3D: paths of the degraded reef after each event
    fig3 = go.Figure()
    for event in events:
        _, path = runs[(event, "degraded (0.19)")]
        fig3.add_trace(go.Scatter3d(x=path[:, 0], y=path[:, 1], z=path[:, 2], mode="lines",
                                    line=dict(width=5), name=event))
    for x, label in [(homes[0], "coral home"), (homes[1], "seaweed home")]:
        fig3.add_trace(go.Scatter3d(x=[x[0]], y=[x[1]], z=[x[2]], mode="markers+text", text=[label],
                                    marker=dict(size=7, color="black"), name=label))
    fig3.update_layout(title="Degraded reef (grazing 0.19): path after each event - drag to rotate",
                       scene=dict(xaxis_title="Coral (C)", yaxis_title="Seaweed (M)",
                                  zaxis_title="Herbivores (H)"))
    fig3.write_html(OUT / "step09_events_3d.html", include_plotlyjs="cdn")
    print(f"\nFigures saved in {OUT}/")


if __name__ == "__main__":
    main()
