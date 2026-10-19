"""Export model results the website draws itself (charts that follow the page's theme).

Writes website/data/reef_data.js, which sets window.REEF_DATA. (A .js file, not
.json, so the page also works when opened straight from disk.)

Run from the project root (about 1.5 minutes):
    python -m tools.export_web_data
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.basins import seaweed_tolerance
from src.equilibria import find_equilibria
from src.parameters import BASELINE, make_parameters
from src.sensitivity import tipping_grazing
from src.simulate import simulate
from src.stability import classify
from steps.step14_rescue_site_d import PLANS, plan_setup

OUT = Path("website/data/reef_data.js")


def grazing_sweep():
    """Step 3: final coral cover from a coral start and a seaweed start."""
    starts = {"coral_start": (0.70, 0.02, 0.50), "seaweed_start": (0.05, 0.80, 0.50)}
    G = np.round(np.arange(0.0, 0.501, 0.005), 3)
    out = {"grazing": G.tolist()}
    for name, x0 in starts.items():
        out[name] = [round(float(simulate(x0, make_parameters(G_max=g), years=3000,
                                          points_per_year=1)[1][-1, 0]), 4) for g in G]
    return out


def tolerance_wall():
    """Step 5: max starting seaweed that still recovers, grazing 0.19."""
    p = make_parameters(G_max=0.19)
    pts = [x for x in find_equilibria(p) if x[2] > 1e-6]
    homes = [x for x in pts if classify(x, p)[1].startswith("STABLE")]
    coral = np.round(np.arange(0.0, 0.951, 0.05), 3)
    herb = np.round(np.arange(0.1, 1.51, 0.1), 3)
    wall = [[round(seaweed_tolerance(C, H, p, homes), 4) for C in coral] for H in herb]
    return {"coral": coral.tolist(), "herbivores": herb.tolist(), "max_seaweed": wall}


def rescue_paths():
    """Step 14: coral cover over 150 years under each plan (no shaking)."""
    clean = pd.read_csv("data/synthetic_clean/true_states.csv")
    last = clean[clean.site == "D Fished and polluted"].iloc[-1]
    x_now = np.array([last.coral, last.seaweed, last.herbivores])
    out = {"years": None, "plans": {}}
    for plan in PLANS:
        p, x0 = plan_setup(plan, x_now)
        t, X = simulate(x0, p, years=150, points_per_year=2)
        out["years"] = t.tolist()
        out["plans"][plan] = {"coral": np.round(X[:, 0], 4).tolist(),
                              "margin": round(p["G_max"] - float(tipping_grazing(p)), 3)}
    return out


def main():
    data = {
        "baseline": BASELINE,
        "tipping_point_baseline": round(float(tipping_grazing(make_parameters())), 4),
        "sweep": grazing_sweep(),
        "wall": tolerance_wall(),
        "rescue": rescue_paths(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("window.REEF_DATA = " + json.dumps(data) + ";\n", encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
