"""Find resting points (equilibria): reef states where every rate is zero."""

import numpy as np
from scipy.optimize import root

from src.model import rhs


def is_feasible(x, tol=1e-9):
    C, M, H = x
    return C >= -tol and M >= -tol and H >= -tol and C + M <= 1 + tol


def find_equilibria(p, n_guesses=400, seed=0, tol=1e-10):
    """Start an equation solver from many random reef states and collect
    every distinct, ecologically possible point where dC/dt = dM/dt = dH/dt = 0.
    """
    rng = np.random.default_rng(seed)
    found = []
    for _ in range(n_guesses):
        C, M = rng.dirichlet([1, 1, 1])[:2]   # random C, M with C + M <= 1
        H = rng.uniform(0, 1.5)
        with np.errstate(all="ignore"):  # the solver may try impossible states on the way
            sol = root(lambda x: rhs(0, x, p), (C, M, H), method="hybr", tol=tol)
        x = sol.x
        if not sol.success or not is_feasible(x):
            continue
        if np.max(np.abs(rhs(0, x, p))) > 1e-9:
            continue
        x = np.clip(x, 0, None)
        if all(np.max(np.abs(x - y)) > 1e-5 for y in found):
            found.append(x)
    return sorted(found, key=lambda x: -x[0])  # most coral first
