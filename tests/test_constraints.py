"""The reef can never reach an impossible state (negative cover, or coral + seaweed > 100%)."""

import numpy as np

from src.model import rhs
from src.parameters import make_parameters
from src.simulate import simulate
from src.stochastic import simulate_noisy

TOL = 1e-9


def feasible(X):
    C, M, H = X[..., 0], X[..., 1], X[..., 2]
    return np.all(C >= -TOL) and np.all(M >= -TOL) and np.all(H >= -TOL) and np.all(C + M <= 1 + TOL)


def test_rules_point_inwards_at_every_edge():
    p = make_parameters(G_max=0.19, N=1.5, F_H=0.03)
    rng = np.random.default_rng(0)
    for _ in range(2000):
        a, h = rng.uniform(0, 1), rng.uniform(0, 1.5)
        assert rhs(0, (0, a, h), p)[0] >= 0           # no coral: coral can't fall
        assert rhs(0, (a, 0, h), p)[1] >= 0           # no seaweed: seaweed can't fall
        assert rhs(0, (a, 1 - a, 0), p)[2] >= 0       # no fish: fish can't fall
        d = rhs(0, (a, 1 - a, h), p)
        assert d[0] + d[1] <= TOL                     # full reef: can't get fuller


def test_simulations_stay_possible():
    rng = np.random.default_rng(1)
    for _ in range(30):
        C, M = rng.dirichlet([1, 1, 1])[:2]
        p = make_parameters(G_max=rng.uniform(0.05, 0.6), N=rng.uniform(0.5, 2), F_H=rng.uniform(0, 0.06))
        _, X = simulate((C, M, rng.uniform(0, 1.5)), p, years=200, points_per_year=2)
        assert feasible(X)


def test_noisy_simulations_stay_possible():
    p = make_parameters(G_max=0.19)
    _, path = simulate_noisy((0.5, 0.4, 0.5), p, 50, 0.08, 200, np.random.default_rng(2))
    assert feasible(np.moveaxis(path, 1, -1))
