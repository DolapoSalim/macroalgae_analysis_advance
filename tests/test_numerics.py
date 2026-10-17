"""Results don't depend on computing choices (step sizes, tolerances) and are repeatable."""

import numpy as np
from scipy.integrate import solve_ivp

from src.model import rhs
from src.parameters import make_parameters
from src.simulate import simulate
from src.stochastic import simulate_noisy

X0 = (0.25, 0.15, 0.5)


def test_tighter_solver_tolerance_gives_same_trajectory():
    p = make_parameters(N=1.2)
    t, X = simulate(X0, p, years=100)
    tight = solve_ivp(rhs, (0, 100), X0, args=(p,), t_eval=t, rtol=1e-11, atol=1e-13)
    assert np.max(np.abs(X - tight.y.T)) < 1e-6


def test_simple_stepping_matches_accurate_solver():
    """The small-step method used for noisy runs (with the noise switched off)
    agrees with the accurate solver."""
    p = make_parameters(N=1.2)
    t, X = simulate(X0, p, years=50, points_per_year=1)
    for dt, tol in [(0.05, 5e-3), (0.01, 1e-3)]:
        _, path = simulate_noisy(X0, p, 50, 0.0, 1, np.random.default_rng(0), dt=dt)
        assert np.max(np.abs(path[:, :, 0] - X)) < tol


def test_same_seed_same_result():
    p = make_parameters(G_max=0.19)
    a = simulate_noisy(X0, p, 20, 0.02, 10, np.random.default_rng(42))[1]
    b = simulate_noisy(X0, p, 20, 0.02, 10, np.random.default_rng(42))[1]
    assert np.array_equal(a, b)
