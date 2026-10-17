"""Fit the model's biology parameters to survey data (guide §19).

What is known to the 'scientist':   each site's fishing F_H and nutrients N,
                                    the storm year, and the survey data.
What must be estimated:             13 biology parameters (shared by all sites),
                                    plus each site's reef state at the start
                                    and just after the storm (the storm damage
                                    itself is not known).

Fit = make the model's predictions as close as possible to the recorded values,
measured in units of typical measurement error:
    cover:      (recorded - predicted) / 0.03
    herbivores: (log recorded - log predicted) / 0.2
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from src.model import rhs
from src.parameters import make_parameters

FITTED = ["r_C", "r_rec", "mu_C", "a_CM", "r_M", "beta_N", "mu_M",
          "G_max", "K_G", "r_H", "K_H", "q_M", "mu_H"]
COVER_SD = 0.03
HERBIVORE_SD = 0.2
LOWER = 1e-3
UPPER = 5.0


def predict(theta, x0, site_settings, times, t_start):
    """Deterministic model prediction at the given times."""
    p = make_parameters(**dict(zip(FITTED, theta)), **site_settings)
    if len(times) == 0:
        return np.empty((0, 3))
    sol = solve_ivp(rhs, (t_start, times.max() + 1e-9), x0, args=(p,), t_eval=times,
                    method="LSODA", rtol=1e-6, atol=1e-8)
    if not sol.success or sol.y.shape[1] != len(times):
        return np.full((len(times), 3), np.nan)
    return sol.y.T


def unpack(z, n_segments):
    """z = [log biology params (13), then 3 start values per segment (C, M, log H)]."""
    theta = np.exp(z[:len(FITTED)])
    starts = z[len(FITTED):].reshape(n_segments, 3).copy()
    starts[:, 2] = np.exp(starts[:, 2])
    return theta, starts


def residuals(z, segments):
    theta, starts = unpack(z, len(segments))
    out = []
    for seg, x0 in zip(segments, starts):
        pred = predict(theta, x0, seg["settings"], seg["times"], seg["t_start"])
        if np.isnan(pred).any():
            pred = np.full_like(pred, 0.5)          # failed run -> big penalty
        obs = seg["obs"]
        r_cover = (obs[:, :2] - pred[:, :2]) / COVER_SD
        with np.errstate(divide="ignore", invalid="ignore"):
            r_herb = (np.log(obs[:, 2]) - np.log(np.clip(pred[:, 2], 1e-6, None))) / HERBIVORE_SD
        r = np.column_stack([r_cover, r_herb])
        out.append(r[~np.isnan(obs)])
    return np.concatenate(out)


def build_segments(observed, sites_meta, storm_year):
    """Split each site's record into 'before storm' and 'after storm' pieces."""
    segments = []
    for site, info in sites_meta.items():
        d = observed[observed.site == site]
        for t0, t1 in [(0.0, storm_year), (storm_year, np.inf)]:
            part = d[(d.year >= t0) & (d.year < t1)] if t0 == 0 else d[d.year >= t0]
            segments.append({
                "site": site, "t_start": t0,
                "settings": {"F_H": info["F_H"], "N": info["N"]},
                "times": part.year.to_numpy(),
                "obs": part[["coral", "seaweed", "herbivores"]].to_numpy(),
            })
    return segments


def first_guess_starts(segments):
    """Start state of each segment = average of its first few recorded values
    (for each variable, the first 3 values that were actually recorded)."""
    z = []
    for seg in segments:
        C, M, H = (np.mean(col[~np.isnan(col)][:3]) for col in seg["obs"].T)
        z += [C, M, np.log(H)]
    return np.array(z)


def fit(segments, theta_guess, max_evaluations=200):
    z0 = np.concatenate([np.log(theta_guess), first_guess_starts(segments)])
    n_seg = len(segments)
    lower = np.concatenate([np.full(len(FITTED), np.log(LOWER)), np.tile([0, 0, np.log(0.01)], n_seg)])
    upper = np.concatenate([np.full(len(FITTED), np.log(UPPER)), np.tile([1, 1, np.log(3.0)], n_seg)])
    z0 = np.clip(z0, lower + 1e-6, upper - 1e-6)
    result = least_squares(residuals, z0, args=(segments,), bounds=(lower, upper),
                           x_scale="jac", max_nfev=max_evaluations, verbose=0)
    theta, starts = unpack(result.x, n_seg)
    return theta, starts, result
