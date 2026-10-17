"""Which fitted numbers can the data actually pin down? (guide §19.5)

1. Profile likelihood: lock one parameter at a series of values, refit all the
   others each time, and record how much worse the best possible fit becomes.
2. Fisher information: from how sharply the misfit responds to each parameter
   (and to combinations of parameters) right at the best fit.
"""

import numpy as np
from scipy.optimize import least_squares

from src.estimation import FITTED, LOWER, UPPER, residuals


def pack(theta, starts):
    """Inverse of estimation.unpack: back to the search vector z."""
    s = np.array(starts, dtype=float).copy()
    s[:, 2] = np.log(s[:, 2])
    return np.concatenate([np.log(theta), s.ravel()])


def bounds(n_segments):
    lower = np.concatenate([np.full(len(FITTED), np.log(LOWER)), np.tile([0, 0, np.log(0.01)], n_segments)])
    upper = np.concatenate([np.full(len(FITTED), np.log(UPPER)), np.tile([1, 1, np.log(3.0)], n_segments)])
    return lower, upper


def misfit_with_locked(segments, z_start, i, log_value):
    """Best misfit when parameter i is locked at exp(log_value) and everything else is refitted.
    Returns (misfit, z of that refit)."""
    free = np.ones(len(z_start), dtype=bool)
    free[i] = False
    lower, upper = bounds(len(segments))

    def res_free(zf):
        z = z_start.copy()
        z[free] = zf
        z[i] = log_value
        return residuals(z, segments)

    zf0 = np.clip(z_start[free], lower[free] + 1e-9, upper[free] - 1e-9)
    out = least_squares(res_free, zf0, bounds=(lower[free], upper[free]), x_scale="jac", max_nfev=100)
    z = z_start.copy()
    z[free] = out.x
    z[i] = log_value
    return 2 * out.cost, z


def profile(segments, z_best, i, factors):
    """Misfit across locked values best_value * factor, walking outwards from the
    best fit so each refit starts from its neighbour (faster and more reliable)."""
    factors = np.asarray(factors, dtype=float)
    misfits = np.full(len(factors), np.nan)
    centre = int(np.argmin(np.abs(np.log(factors))))
    for direction in (range(centre, len(factors)), range(centre, -1, -1)):
        z = z_best.copy()
        for k in direction:
            misfits[k], z = misfit_with_locked(segments, z, i, z_best[i] + np.log(factors[k]))
    return misfits


def residual_jacobian(segments, z, h=1e-4):
    """How each residual changes when each entry of z is nudged (central differences)."""
    r0 = residuals(z, segments)
    J = np.zeros((len(r0), len(z)))
    for j in range(len(z)):
        dz = np.zeros(len(z))
        dz[j] = h
        J[:, j] = (residuals(z + dz, segments) - residuals(z - dz, segments)) / (2 * h)
    return J


def biology_information(J, n_bio=len(FITTED)):
    """Fisher information for the biology parameters (log scale), after allowing the
    nuisance start values to adjust (Schur complement)."""
    F = J.T @ J
    A, B, D = F[:n_bio, :n_bio], F[:n_bio, n_bio:], F[n_bio:, n_bio:]
    return A - B @ np.linalg.solve(D, B.T)
