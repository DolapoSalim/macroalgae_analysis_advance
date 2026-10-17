"""Run the model forward in time."""

import numpy as np
from scipy.integrate import solve_ivp

from src.model import rhs


def simulate(x0, p, years=100, points_per_year=10):
    """Follow the reef from starting state x0 = (C, M, H) for a number of years.

    Returns (t, X): t is an array of times, X has one row per time and
    columns C, M, H.
    """
    t_eval = np.linspace(0, years, years * points_per_year + 1)
    sol = solve_ivp(rhs, (0, years), x0, args=(p,), t_eval=t_eval,
                    method="RK45", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"Simulation failed: {sol.message}")
    return sol.t, sol.y.T
