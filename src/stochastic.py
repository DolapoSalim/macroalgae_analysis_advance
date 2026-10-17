"""Noise ("shaking"): the model plus small random kicks at every time step.

Method: Euler–Maruyama. Each small time step dt does
    x_next = x + f(x) * dt                     (the normal rules)
               + sigma * sqrt(dt) * random     (a random kick)
where `random` is drawn from a standard normal (bell-curve) distribution.

All replicates are simulated at once: X has shape (3, n_replicates).
"""

import numpy as np

from src.model import rhs


def keep_feasible(X):
    """Clip back into the possible region: no negatives, coral + seaweed <= 1."""
    X = np.clip(X, 0, None)
    total = X[0] + X[1]
    over = total > 1
    X[:2, over] /= total[over]
    return X


def simulate_noisy(x0, p, years, sigma, n_replicates, rng, dt=0.05, record_every=1.0):
    """Run n_replicates noisy copies of the same reef. Returns (t, path),
    path has shape (n_records, 3, n_replicates)."""
    X = np.repeat(np.asarray(x0, dtype=float)[:, None], n_replicates, axis=1)
    n_steps = int(round(years / dt))
    every = int(round(record_every / dt))
    t, path = [0.0], [X.copy()]
    for i in range(1, n_steps + 1):
        kick = sigma * np.sqrt(dt) * rng.normal(size=X.shape)
        X = keep_feasible(X + rhs(0, X, p) * dt + kick)
        if i % every == 0:
            t.append(i * dt)
            path.append(X.copy())
    return np.array(t), np.array(path)


def settle(X, p, years=1500, dt=0.1):
    """Switch the noise off and let every replicate roll to its resting point."""
    X = X.copy()
    for _ in range(int(round(years / dt))):
        X = keep_feasible(X + rhs(0, X, p) * dt)
    return X


def which_home_many(X, homes):
    """Index of the nearest home for each column of X."""
    d = np.array([np.max(np.abs(X - h[:, None]), axis=0) for h in homes])
    return np.argmin(d, axis=0)
