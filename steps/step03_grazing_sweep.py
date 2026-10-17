"""Step 3 — Search for a tipping point by turning down the grazing dial.

For each grazing strength G_max we start TWO reefs:
  - a coral reef   (lots of coral, little seaweed)
  - a seaweed reef (little coral, lots of seaweed)
run each for a very long time, and record where it finally settles.

If both reefs always end in the same place -> one "home", no tipping point.
If, for some grazing values, they end in DIFFERENT places -> two "homes"
(the reef's past decides its future): evidence of a tipping point.

Run from the project root:
    python -m steps.step03_grazing_sweep
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.model import rhs
from src.parameters import make_parameters
from src.simulate import simulate

STARTS = {
    "starts as coral reef":   (0.70, 0.02, 0.50),
    "starts as seaweed reef": (0.05, 0.80, 0.50),
}
GRAZING_VALUES = np.round(np.arange(0.0, 0.501, 0.005), 3)
YEARS = 3000  # long enough for every reef to settle
FIGURE = Path("figures/step03_grazing_sweep.png")


def final_state(x0, p):
    _, X = simulate(x0, p, years=YEARS, points_per_year=1)
    x_end = X[-1]
    still_changing = np.max(np.abs(rhs(0, x_end, p)))
    return x_end, still_changing


def main():
    results = {name: [] for name in STARTS}
    worst_change = 0.0
    for G in GRAZING_VALUES:
        p = make_parameters(G_max=G)
        for name, x0 in STARTS.items():
            x_end, change = final_state(x0, p)
            results[name].append(x_end)
            worst_change = max(worst_change, change)
    results = {name: np.array(v) for name, v in results.items()}

    coral_from_coral = results["starts as coral reef"][:, 0]
    coral_from_seaweed = results["starts as seaweed reef"][:, 0]
    differ = np.abs(coral_from_coral - coral_from_seaweed) > 0.05

    print("grazing   final coral        final coral")
    print("G_max     (coral start)      (seaweed start)")
    for G, a, b, d in zip(GRAZING_VALUES, coral_from_coral, coral_from_seaweed, differ):
        if round(G * 1000) % 50 == 0 or d:  # print every 0.05, plus all 'two homes' rows
            flag = "  <-- TWO HOMES" if d else ""
            print(f"{G:5.3f}     {a:5.2f}              {b:5.2f}{flag}")

    print(f"\nLargest rate of change at the end of any run: {worst_change:.1e} per year "
          "(near zero = every reef had settled)")
    if differ.any():
        lo, hi = GRAZING_VALUES[differ].min(), GRAZING_VALUES[differ].max()
        print(f"Two homes exist for grazing between {lo:.3f} and {hi:.3f}")
    else:
        print("The two starting reefs always ended in the same place: no two homes found")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    if differ.any():
        ax.axvspan(lo - 0.0025, hi + 0.0025, color="#ffcc66", alpha=0.4,
                   label="two homes (depends on the past)")
    ax.plot(GRAZING_VALUES, coral_from_coral, "o", ms=4, color="#1f77b4",
            label="starts as coral reef")
    ax.plot(GRAZING_VALUES, coral_from_seaweed, "x", ms=5, color="#d62728",
            label="starts as seaweed reef")
    ax.axvline(0.40, color="grey", ls="--", lw=1)
    ax.text(0.402, 0.05, "our baseline\n(Steps 1–2)", fontsize=8, color="grey")
    ax.set_xlabel("Grazing strength  G_max  (how much seaweed herbivores can eat)")
    ax.set_ylabel("Coral cover where the reef finally settles")
    ax.set_title("Turning down the grazing dial")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(loc="center right")
    fig.tight_layout()
    FIGURE.parent.mkdir(exist_ok=True)
    fig.savefig(FIGURE, dpi=120)
    print(f"Figure saved to {FIGURE}")


if __name__ == "__main__":
    main()
