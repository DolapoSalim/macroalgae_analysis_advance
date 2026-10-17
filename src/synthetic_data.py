"""Turn 'perfect' simulated reefs into realistic, messy survey data (guide §17–18).

Truth:        noisy simulation (everyday shaking), with a regional storm.
Observation:  surveys every ~3 months, with
                - measurement error (beta for cover, lognormal for herbivores)
                - surveys a few days/weeks off schedule
                - randomly missed surveys
                - herbivores counted only on every second survey
                - an optional funding gap with no surveys at all
"""

import numpy as np

from src.disturbances import apply_shock
from src.stochastic import simulate_noisy


def simulate_truth(x0, p, years, storm_year, sigma, rng, dt=0.05):
    """True reef path (one replicate) with a strong storm at storm_year.
    Returns t (every dt) and X with columns C, M, H."""
    t1, path1 = simulate_noisy(x0, p, storm_year, sigma, 1, rng, dt=dt, record_every=dt)
    x_storm = apply_shock(path1[-1, :, 0], "Strong storm")
    t2, path2 = simulate_noisy(x_storm, p, years - storm_year, sigma, 1, rng, dt=dt, record_every=dt)
    t = np.concatenate([t1, storm_year + t2[1:]])
    X = np.vstack([path1[:, :, 0], path2[1:, :, 0]])
    return t, X


def survey_times(years, rng, every=0.25, jitter=0.04):
    nominal = np.arange(0, years + 1e-9, every)
    return np.clip(nominal + rng.uniform(-jitter, jitter, nominal.size), 0, years)


def true_at(times, t, X):
    return np.column_stack([np.interp(times, t, X[:, j]) for j in range(3)])


def observe(true_values, rng, cover_precision=150, herbivore_error=0.2):
    """Add measurement error.
    Cover: beta distribution centred on the true value (stays between 0 and 1).
    Herbivores: lognormal, i.e. about +/-20% multiplicative error (never negative)."""
    obs = np.empty_like(true_values)
    for j in (0, 1):
        m = np.clip(true_values[:, j], 1e-3, 1 - 1e-3)
        obs[:, j] = rng.beta(m * cover_precision, (1 - m) * cover_precision)
    obs[:, 2] = true_values[:, 2] * rng.lognormal(-herbivore_error ** 2 / 2, herbivore_error,
                                                   len(true_values))
    return obs


def knock_out(times, obs, rng, missing_share=0.15, gap=None):
    """Remove surveys: random misses, fish counted every 2nd survey, optional gap (start, end)."""
    obs = obs.copy()
    missed = rng.uniform(size=len(times)) < missing_share
    if gap is not None:
        missed |= (times >= gap[0]) & (times < gap[1])
    obs[missed, :] = np.nan
    obs[1::2, 2] = np.nan
    return obs
