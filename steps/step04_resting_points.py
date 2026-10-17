"""Step 4 — Find the resting points directly and test whether each is stable.

Two tests for every resting point:
  1. Maths test (Jacobian eigenvalues): do tiny nudges shrink or grow?
  2. Simulation test: actually nudge the reef 20 times and watch where it goes.

Run from the project root:
    python -m steps.step04_resting_points
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.equilibria import find_equilibria
from src.parameters import make_parameters
from src.simulate import simulate
from src.stability import classify

GRAZING_CASES = {0.40: "baseline (one home in Step 3)", 0.19: "inside the two-homes zone"}
N_NUDGES = 20
NUDGE_SIZE = 0.01
FIGURE = Path("figures/step04_resting_points.png")


def nudge(x, rng):
    """Push the reef a tiny distance in a random direction (staying feasible)."""
    direction = rng.normal(size=3)
    y = x + NUDGE_SIZE * direction / np.linalg.norm(direction)
    y = np.clip(y, 0, None)
    if y[0] + y[1] > 1:
        y[:2] *= 0.999 / (y[0] + y[1])
    return y


def nearest(x, points):
    return int(np.argmin([np.max(np.abs(x - q)) for q in points]))


def main():
    rng = np.random.default_rng(4)
    trails = []  # for the figure (G = 0.19 only)

    for G, note in GRAZING_CASES.items():
        p = make_parameters(G_max=G)
        points = find_equilibria(p)
        with_fish = [x for x in points if x[2] > 1e-6]
        no_fish = [x for x in points if x[2] <= 1e-6]
        homes = [x for x in with_fish if classify(x, p)[1].startswith("STABLE")]

        print(f"\n================ Grazing G_max = {G}  ({note}) ================")
        for i, x in enumerate(with_fish, 1):
            eig, label = classify(x, p)
            print(f"\nResting point {i}: coral {x[0]:.2f}, seaweed {x[1]:.2f}, herbivores {x[2]:.2f}")
            print(f"  Maths test   -> eigenvalues (real part): "
                  f"{', '.join(f'{v:+.3f}' for v in eig.real)}  => {label}")

            # Simulation test: 20 tiny random nudges, 1000 years each
            ends = []
            for _ in range(N_NUDGES):
                y0 = nudge(x, rng)
                t, X = simulate(y0, p, years=1000, points_per_year=2)
                ends.append(nearest(X[-1], homes))
                if G == 0.19 and label.startswith("UNSTABLE"):
                    trails.append((X, ends[-1]))
            counts = [ends.count(k) for k in range(len(homes))]
            summary = ", ".join(
                f"{c} went to the home with {homes[k][0]:.0%} coral" for k, c in enumerate(counts) if c)
            print(f"  Nudge test   -> of {N_NUDGES} nudges: {summary}")

        if no_fish:
            print(f"\n  (+ {len(no_fish)} resting point(s) with ZERO herbivores -- all unstable: "
                  "any fish arriving would multiply. Ignored from here on.)")

    # Figure: the G = 0.19 landscape seen from above (coral vs seaweed)
    p = make_parameters(G_max=0.19)
    points = [x for x in find_equilibria(p) if x[2] > 1e-6]
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    colours = {0: "#1f77b4", 1: "#d62728"}  # home 0 = coral home, home 1 = seaweed home
    for X, k in trails:
        ax.plot(X[:, 0], X[:, 1], color=colours.get(k, "grey"), lw=1, alpha=0.7)
    for x in points:
        stable = classify(x, p)[1].startswith("STABLE")
        ax.plot(x[0], x[1], "o", ms=11, mfc="black" if stable else "white", mec="black", mew=2)
        name = "coral home" if (stable and x[0] > 0.5) else "seaweed home" if stable else "balance point"
        ax.annotate(name, (x[0], x[1]), textcoords="offset points", xytext=(8, 8), fontsize=10)
    ax.plot([], [], color=colours[0], label="nudged balance point -> coral home")
    ax.plot([], [], color=colours[1], label="nudged balance point -> seaweed home")
    ax.plot([], [], "o", mfc="black", mec="black", label="stable (home)")
    ax.plot([], [], "o", mfc="white", mec="black", mew=2, label="unstable (balance point)")
    ax.set_xlabel("Coral cover (C)")
    ax.set_ylabel("Seaweed cover (M)")
    ax.set_title("Grazing = 0.19: two homes and the balance point between them")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    FIGURE.parent.mkdir(exist_ok=True)
    fig.savefig(FIGURE, dpi=120)
    print(f"\nFigure saved to {FIGURE}")


if __name__ == "__main__":
    main()
