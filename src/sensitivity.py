"""Uncertainty and sensitivity analysis (guide §21–22).

Every function here works on MANY parameter sets at once: each parameter in p
can be an array (one value per parameter set).
"""

import numpy as np

from src.disturbances import SCENARIOS
from src.stochastic import keep_feasible
from src.model import rhs


def coral_home(p):
    """The clean coral home (M = 0), worked out exactly.
    Coral:  r_C C (1-C) + r_rec (1-C) - mu_C C = 0   ->  a quadratic in C
    Fish:   H = K_H (1 - (mu_H + F_H) / r_H)"""
    a, b, c = -p["r_C"], p["r_C"] - p["r_rec"] - p["mu_C"], p["r_rec"]
    C = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    H = p["K_H"] * (1 - (p["mu_H"] + p["F_H"]) / p["r_H"])
    return C, np.clip(H, 0, None)


def tipping_grazing(p, cap=1.0):
    """Grazing strength G_max below which seaweed can invade a clean coral reef.

    At the coral home, a tiny patch of seaweed grows if
        (r_M + beta_N N)(1 - C) + a_CM C - mu_M  >  G_max H / K_G
    so the tipping point is G_max = K_G * [left side] / H.
    0 = seaweed can never invade; `cap` = no realistic grazing can stop it
    (including the case where the fish population cannot survive)."""
    C, H = coral_home(p)
    seaweed_push = (p["r_M"] + p["beta_N"] * p["N"]) * (1 - C) + p["a_CM"] * C - p["mu_M"]
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(H > 1e-3, p["K_G"] * seaweed_push / H, cap)
    return np.clip(g, 0, cap)


def coral_after_storm(p, years_after=50, dt=0.05):
    """Coral cover `years_after` years after a strong storm hits the CLEAN CORAL HOME
    (no shaking). Starting every reef version from its own coral home means the
    storm, not an already-degraded starting state, is what is being tested."""
    C, H = coral_home(p)
    coral_loss, seaweed_arrival, herbivore_loss, _ = SCENARIOS["Strong storm"]
    C = (1 - coral_loss) * C
    X = np.vstack([C, np.minimum(seaweed_arrival, 1 - C), (1 - herbivore_loss) * H])
    for _ in range(int(round(years_after / dt))):
        X = keep_feasible(X + rhs(0, X, p) * dt)
    return X[0]


def sobol_indices(f_A, f_B, f_AB, n_boot=200, rng=None):
    """First-order (S1) and total (ST) Sobol indices with bootstrap 95% ranges.

    S1_i: share of the output's variation caused by parameter i ON ITS OWN.
    ST_i: share caused by parameter i INCLUDING its teamwork with others.
    Estimators: Saltelli (2010) for S1, Jansen (1999) for ST."""
    rng = rng or np.random.default_rng(0)

    def compute(idx):
        a, b, ab = f_A[idx], f_B[idx], f_AB[:, idx]
        var = np.var(np.concatenate([a, b]))
        s1 = np.mean(b * (ab - a), axis=1) / var
        st = 0.5 * np.mean((a - ab) ** 2, axis=1) / var
        return s1, st

    n = len(f_A)
    s1, st = compute(np.arange(n))
    boots = [compute(rng.integers(0, n, n)) for _ in range(n_boot)]
    s1_ci = np.percentile([b[0] for b in boots], [2.5, 97.5], axis=0)
    st_ci = np.percentile([b[1] for b in boots], [2.5, 97.5], axis=0)
    return s1, st, s1_ci, st_ci
