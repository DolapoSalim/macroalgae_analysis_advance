"""Step 2 — First simulation: follow two reefs through 100 years.

Run from the project root:
    python -m steps.step02_first_simulation
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # save figures to file without opening a window
import matplotlib.pyplot as plt

from src.parameters import make_parameters
from src.simulate import simulate

SCENARIOS = {
    # name: ((C, M, H) at the start, nutrient level N)
    "Healthy reef": ((0.60, 0.15, 0.70), 0.5),
    "After storm":  ((0.25, 0.15, 0.50), 1.2),
}
REPORT_YEARS = [0, 5, 10, 25, 50, 100]
FIGURE = Path("figures/step02_trajectories.png")


def main():
    fig, axes = plt.subplots(1, len(SCENARIOS), figsize=(11, 4), sharey=True)

    for ax, (name, (x0, N)) in zip(axes, SCENARIOS.items()):
        p = make_parameters(N=N)
        t, X = simulate(x0, p, years=100)

        print(f"\n=== {name} (nutrients N={N}) ===")
        print("  year   coral  seaweed  herbivores")
        for year in REPORT_YEARS:
            i = year * 10  # 10 points per year
            C, M, H = X[i]
            print(f"  {year:>4}   {C:5.2f}   {M:5.2f}     {H:5.2f}")

        ax.plot(t, X[:, 0], color="#1f77b4", label="Coral (C)")
        ax.plot(t, X[:, 1], color="#d62728", label="Seaweed (M)")
        ax.plot(t, X[:, 2], color="#2ca02c", label="Herbivores (H)")
        ax.set_title(f"{name} (N = {N})")
        ax.set_xlabel("Years")
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("Cover (C, M) or biomass (H)")
    axes[0].legend()
    fig.tight_layout()
    FIGURE.parent.mkdir(exist_ok=True)
    fig.savefig(FIGURE, dpi=120)
    print(f"\nFigure saved to {FIGURE}")


if __name__ == "__main__":
    main()
