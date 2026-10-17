"""The fitting machinery is consistent: with perfect data and the true numbers,
the misfit is zero, and wrong numbers make it worse."""

import numpy as np
import pandas as pd

from src.disturbances import apply_shock
from src.estimation import FITTED, build_segments, residuals
from src.identifiability import pack
from src.parameters import BASELINE, make_parameters
from src.simulate import simulate

SITES = {"A": {"F_H": 0.0, "N": 0.5}, "D": {"F_H": 0.04, "N": 1.5}}
STORM = 8


def perfect_data():
    """Noise-free surveys every 3 months, true biology, storm in year 8."""
    rows, starts = [], []
    for site, s in SITES.items():
        p = make_parameters(**s)
        x0 = np.array([0.6, 0.1, 0.4])
        t1, X1 = simulate(x0, p, years=STORM, points_per_year=4)
        x_storm = apply_shock(X1[-1], "Strong storm")
        t2, X2 = simulate(x_storm, p, years=12, points_per_year=4)
        starts += [x0, x_storm]
        for t, x in list(zip(t1[:-1], X1[:-1])) + list(zip(STORM + t2, X2)):
            rows.append({"site": site, "year": t, "coral": x[0], "seaweed": x[1], "herbivores": x[2]})
    return pd.DataFrame(rows), np.array(starts)


def test_true_numbers_fit_perfect_data():
    data, starts = perfect_data()
    segments = build_segments(data, SITES, STORM)
    truth = np.array([BASELINE[k] for k in FITTED])
    r = residuals(pack(truth, starts), segments)
    assert np.max(np.abs(r)) < 1e-3


def test_wrong_numbers_fit_worse():
    data, starts = perfect_data()
    segments = build_segments(data, SITES, STORM)
    truth = np.array([BASELINE[k] for k in FITTED])
    good = np.sum(residuals(pack(truth, starts), segments) ** 2)
    wrong = truth.copy()
    wrong[FITTED.index("G_max")] *= 1.3
    bad = np.sum(residuals(pack(wrong, starts), segments) ** 2)
    assert bad > good + 100
