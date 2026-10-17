"""Basin mapping: which home does each starting reef end up in?"""

import numpy as np

from src.simulate import simulate

UNRESOLVED = -1


def which_home(x0, p, homes, years=1500, tol=0.02):
    """Simulate from x0 and return the index of the home it settles in,
    or UNRESOLVED if it is not close to any home by the end."""
    _, X = simulate(x0, p, years=years, points_per_year=1)
    distances = [np.max(np.abs(X[-1] - h)) for h in homes]
    k = int(np.argmin(distances))
    return k if distances[k] < tol else UNRESOLVED


def grid_of_reefs(step, herbivore_levels):
    """Every (C, M, H) on a grid, keeping only reefs with C + M <= 1."""
    values = np.round(np.arange(0, 1 + 1e-9, step), 6)
    return np.array([(C, M, H) for H in herbivore_levels
                     for C in values for M in values if C + M <= 1 + 1e-9])


def map_basins(starts, p, homes):
    return np.array([which_home(x0, p, homes) for x0 in starts])


def seaweed_tolerance(C, H, p, homes, coral_home=0, tol=0.002):
    """The most seaweed a reef with coral C and herbivores H can start with
    and still reach the coral home (found by repeatedly halving the interval)."""
    lo, hi = 0.0, 1.0 - C
    if which_home((C, lo, H), p, homes) != coral_home:
        return 0.0
    if which_home((C, hi, H), p, homes) == coral_home:
        return hi
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if which_home((C, mid, H), p, homes) == coral_home:
            lo = mid
        else:
            hi = mid
    return lo
