"""Step 10 — Make realistic, messy field data from six virtual reef sites.

Every site follows the SAME biology (baseline parameters). Sites differ only in
things a manager could measure: fishing pressure F_H and nutrient level N,
plus where the reef starts. A strong regional storm hits every site in year 8.

Outputs:
  data/synthetic_clean/true_states.csv          the truth at each survey date
  data/synthetic_observed/observations.csv      what the 'divers' recorded
  data/metadata/sites.json                      site settings + true parameters
  figures/step10_fake_field_data.png

Run from the project root:
    python -m steps.step10_fake_field_data
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.parameters import BASELINE, make_parameters
from src.synthetic_data import knock_out, observe, simulate_truth, survey_times, true_at

YEARS = 20
STORM_YEAR = 8
SIGMA = 0.01
SEED = 10

SITES = {
    # name: start (C, M, H), fishing F_H, nutrients N, funding gap (years) or None
    "A Protected, clean water": ((0.60, 0.10, 0.70), 0.00, 0.5, None),
    "B Protected, polluted":    ((0.45, 0.25, 0.60), 0.00, 1.5, None),
    "C Lightly fished":         ((0.60, 0.10, 0.40), 0.03, 0.5, None),
    "D Fished and polluted":    ((0.60, 0.10, 0.30), 0.04, 1.5, (13, 14.5)),
    "E Heavily fished":         ((0.50, 0.10, 0.30), 0.06, 0.5, (13, 14.5)),
    "F Seaweed reef, protected": ((0.10, 0.60, 0.30), 0.00, 1.0, None),
}


def main():
    rng = np.random.default_rng(SEED)
    clean_rows, obs_rows, truths = [], [], {}

    for site, (x0, F, N, gap) in SITES.items():
        p = make_parameters(F_H=F, N=N)
        t, X = simulate_truth(x0, p, YEARS, STORM_YEAR, SIGMA, rng)
        times = survey_times(YEARS, rng)
        true_vals = true_at(times, t, X)
        obs = knock_out(times, observe(true_vals, rng), rng, gap=gap)
        truths[site] = (t, X, times, obs)
        for ti, tv, ov in zip(times, true_vals, obs):
            clean_rows.append({"site": site, "year": ti, "coral": tv[0], "seaweed": tv[1], "herbivores": tv[2]})
            obs_rows.append({"site": site, "year": ti, "coral": ov[0], "seaweed": ov[1], "herbivores": ov[2]})

    clean = pd.DataFrame(clean_rows)
    observed = pd.DataFrame(obs_rows)
    for folder in ["data/synthetic_clean", "data/synthetic_observed", "data/metadata"]:
        Path(folder).mkdir(parents=True, exist_ok=True)
    clean.to_csv("data/synthetic_clean/true_states.csv", index=False, float_format="%.5f")
    observed.to_csv("data/synthetic_observed/observations.csv", index=False, float_format="%.5f")
    meta = {
        "description": "Synthetic reef monitoring data. Same biology at every site.",
        "true_parameters": BASELINE,
        "sites": {s: {"start_CMH": v[0], "F_H": v[1], "N": v[2], "funding_gap_years": v[3]}
                  for s, v in SITES.items()},
        "storm": {"year": STORM_YEAR, "type": "Strong storm (see src/disturbances.py)"},
        "everyday_noise_sigma": SIGMA,
        "surveys": "every 0.25 yr +/- 0.04 yr; 15% missed at random; herbivores every 2nd survey",
        "measurement_error": "cover ~ beta(precision 150); herbivores ~ lognormal(sd 0.2)",
        "seed": SEED,
    }
    Path("data/metadata/sites.json").write_text(json.dumps(meta, indent=2))

    # ---- summary
    print(f"{'site':28s} {'surveys':>8s} {'coral':>7s} {'seaweed':>8s} {'fish':>6s}   (number of values recorded)")
    for site in SITES:
        o = observed[observed.site == site]
        print(f"{site:28s} {len(o):8d} {o.coral.notna().sum():7d} {o.seaweed.notna().sum():8d} "
              f"{o.herbivores.notna().sum():6d}")
    both = clean.join(observed[["coral", "seaweed", "herbivores"]], rsuffix="_obs")
    for var in ["coral", "seaweed", "herbivores"]:
        err = (both[f"{var}_obs"] - both[var]).abs().mean()
        print(f"Average measurement error, {var:10s}: {err:.3f}")

    # ---- figure
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), sharex=True, sharey=True)
    colours = {"coral": "#1f77b4", "seaweed": "#d62728", "herbivores": "#2ca02c"}
    for ax, site in zip(axes.ravel(), SITES):
        t, X, times, obs = truths[site]
        for j, var in enumerate(["coral", "seaweed", "herbivores"]):
            ax.plot(t, X[:, j], color=colours[var], lw=1.5, alpha=0.6, label=f"true {var}")
            ax.plot(times, obs[:, j], "o", ms=3, color=colours[var], label=f"recorded {var}")
        ax.axvline(STORM_YEAR, color="grey", ls=":")
        gap = SITES[site][3]
        if gap:
            ax.axvspan(*gap, color="grey", alpha=0.15)
        ax.set_title(f"{site}\n(fishing {SITES[site][1]}, nutrients {SITES[site][2]})", fontsize=10)
        ax.set_ylim(0, 1.2)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("Year of monitoring")
    axes[0, 0].legend(fontsize=7, loc="upper right", ncol=2)
    fig.suptitle("Fake field data: lines = the truth, dots = what divers recorded "
                 "(dotted line = storm, grey band = no funding)")
    fig.tight_layout()
    Path("figures").mkdir(exist_ok=True)
    fig.savefig("figures/step10_fake_field_data.png", dpi=120)
    print("\nData saved in data/, figure saved in figures/")


if __name__ == "__main__":
    main()
