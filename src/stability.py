"""Stability of resting points: does a tiny nudge die out or grow?"""

import numpy as np

from src.model import rhs


def numerical_jacobian(x, p, h=1e-6):
    """The Jacobian: a 3x3 table of how each rate (dC/dt, dM/dt, dH/dt)
    responds to a tiny change in each variable (C, M, H).

    Column j = (rates after nudging variable j up - rates after nudging it down) / (2h).
    """
    x = np.asarray(x, dtype=float)
    J = np.zeros((3, 3))
    for j in range(3):
        step = np.zeros(3)
        step[j] = h
        J[:, j] = (rhs(0, x + step, p) - rhs(0, x - step, p)) / (2 * h)
    return J


def classify(x, p):
    """Return (eigenvalues, label).

    Each eigenvalue is the growth (+) or shrink (-) speed, per year, of a
    nudge in one particular direction. All negative -> every nudge dies out.
    """
    eig = np.linalg.eigvals(numerical_jacobian(x, p))
    real = eig.real
    if np.all(real < -1e-9):
        label = "STABLE (a home)"
    elif np.any(real > 1e-9):
        label = "UNSTABLE (a balance point)"
    else:
        label = "INCONCLUSIVE (on a knife edge)"
    return eig, label
