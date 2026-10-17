"""Lyapunov 'safety certificate' for the coral home.

Score:     V(x) = (x - x*)^T P (x - x*)        ("distance to home", always >= 0)
Its trend: Vdot(x) = 2 (x - x*)^T P f(x)        (negative = moving closer to home)

If Vdot < 0 everywhere inside the zone {V <= c} (within the possible region),
then any reef starting in that zone can never leave it and must end at home.
"""

import numpy as np
from scipy.linalg import solve_continuous_lyapunov

from src.model import rhs


def P_from_model(J, Q):
    """Solve the Lyapunov equation J^T P + P J = -Q for the score's shape P."""
    return solve_continuous_lyapunov(J.T, -Q)


def V(X, x_star, P):
    E = X - x_star[:, None]
    return np.einsum("in,ij,jn->n", E, P, E)


def Vdot(X, x_star, P, p):
    E = X - x_star[:, None]
    return 2 * np.einsum("in,ij,jn->n", E, P, rhs(0, X, p))


def sample_zone(x_star, P, c, n, rng):
    """Random points inside {V <= c} that are ecologically possible.
    Half of them lie exactly on the zone's edge, where problems usually appear."""
    L = np.linalg.cholesky(P)
    u = rng.normal(size=(3, n))
    u /= np.linalg.norm(u, axis=0)
    radius = rng.uniform(0, 1, n) ** (1 / 3)
    radius[: n // 2] = 1.0
    X = x_star[:, None] + np.linalg.solve(L.T, u * radius) * np.sqrt(c)
    ok = (X[0] >= 0) & (X[1] >= 0) & (X[2] >= 0) & (X[0] + X[1] <= 1)
    return X[:, ok]


def zone_is_safe(x_star, P, c, p, rng, n=60_000):
    X = sample_zone(x_star, P, c, n, rng)
    v, vd = V(X, x_star, P), Vdot(X, x_star, P, p)
    return np.all(vd[v > 1e-6 * c] < 0)


def largest_safe_level(x_star, P, p, rng, iterations=40):
    """Biggest c for which the score always decreases inside {V <= c}
    (found by repeatedly halving the search range, on a log scale)."""
    lo, hi = 1e-8, 10.0
    if not zone_is_safe(x_star, P, lo, p, rng):
        return 0.0
    for _ in range(iterations):
        mid = np.sqrt(lo * hi)
        if zone_is_safe(x_star, P, mid, p, rng):
            lo = mid
        else:
            hi = mid
    return lo
