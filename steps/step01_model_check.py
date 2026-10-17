"""Step 1 — Check that the model equations behave sensibly.

No simulation yet. We only ask: at a given reef state, which way is each
variable heading, and which processes are responsible?

Run from the project root:
    python -m steps.step01_model_check
"""

import numpy as np

from src.model import process_rates, rhs
from src.parameters import make_parameters

STATES = {
    # (C, M, H) and the nutrient level N
    "healthy reef": ((0.60, 0.15, 0.70), 0.5),
    "after storm":  ((0.25, 0.15, 0.50), 1.2),
}


def show_budget(name, x, N):
    p = make_parameters(N=N)
    r = process_rates(x, p)
    dC, dM, dH = rhs(0, x, p)
    C, M, H = x
    print(f"\n=== {name}: C={C}, M={M}, H={H}, N={N}  (free space B={1 - C - M:.2f}) ===")
    print("  Coral:")
    print(f"    + growth into free space   {r['coral_growth']:+.4f}")
    print(f"    + recruitment              {r['coral_recruitment']:+.4f}")
    print(f"    - overgrown by algae       {-r['coral_overgrown_by_algae']:+.4f}")
    print(f"    - mortality                {-r['coral_mortality']:+.4f}")
    print(f"    = dC/dt                    {dC:+.4f} per year")
    print("  Macroalgae:")
    print(f"    + growth into free space   {r['algae_growth']:+.4f}")
    print(f"    + overgrowing coral        {r['algae_overgrowing_coral']:+.4f}")
    print(f"    - grazed by herbivores     {-r['algae_grazed']:+.4f}")
    print(f"    - background loss          {-r['algae_loss']:+.4f}")
    print(f"    = dM/dt                    {dM:+.4f} per year")
    print("  Herbivores:")
    print(f"    + growth                   {r['herbivore_growth']:+.4f}")
    print(f"    + food boost               {r['herbivore_food_boost']:+.4f}")
    print(f"    - mortality                {-r['herbivore_mortality']:+.4f}")
    print(f"    - fishing                  {-r['herbivore_fishing']:+.4f}")
    print(f"    = dH/dt                    {dH:+.4f} per year")


def check_boundaries():
    """The model must never push the reef outside the feasible region:
    C >= 0, M >= 0, H >= 0 and C + M <= 1."""
    p = make_parameters()
    rng = np.random.default_rng(1)
    worst = {"C=0 edge: dC/dt": np.inf, "M=0 edge: dM/dt": np.inf,
             "H=0 edge: dH/dt": np.inf, "C+M=1 edge: d(C+M)/dt": -np.inf}
    for _ in range(10_000):
        a, h = rng.uniform(0, 1), rng.uniform(0, 1.5)
        worst["C=0 edge: dC/dt"] = min(worst["C=0 edge: dC/dt"], rhs(0, (0, a, h), p)[0])
        worst["M=0 edge: dM/dt"] = min(worst["M=0 edge: dM/dt"], rhs(0, (a, 0, h), p)[1])
        worst["H=0 edge: dH/dt"] = min(worst["H=0 edge: dH/dt"], rhs(0, (a, 1 - a, 0), p)[2])
        d = rhs(0, (a, 1 - a, h), p)
        worst["C+M=1 edge: d(C+M)/dt"] = max(worst["C+M=1 edge: d(C+M)/dt"], d[0] + d[1])

    print("\n=== Boundary check (10,000 random points on each edge) ===")
    print(f"  smallest dC/dt where C=0        {worst['C=0 edge: dC/dt']:+.4f}  (must be >= 0)")
    print(f"  smallest dM/dt where M=0        {worst['M=0 edge: dM/dt']:+.4f}  (must be >= 0)")
    print(f"  smallest dH/dt where H=0        {worst['H=0 edge: dH/dt']:+.4f}  (must be >= 0)")
    print(f"  largest d(C+M)/dt where C+M=1   {worst['C+M=1 edge: d(C+M)/dt']:+.4f}  (must be <= 0)")


if __name__ == "__main__":
    for name, (x, N) in STATES.items():
        show_budget(name, x, N)
    check_boundaries()
