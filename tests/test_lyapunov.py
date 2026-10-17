"""The safety certificate (Lyapunov function) is built correctly and actually holds."""

import numpy as np

from src.lyapunov import P_from_model, V, Vdot, largest_safe_level, sample_zone
from src.parameters import make_parameters
from src.sensitivity import coral_home
from src.stability import numerical_jacobian
from src.stochastic import settle

P_BASE = make_parameters(G_max=0.19)
C, H = coral_home(P_BASE)
HOME = np.array([C, 0.0, H])


def test_lyapunov_equation_is_solved():
    J = numerical_jacobian(HOME, P_BASE)
    Q = np.diag([1.0, 30.0, 10.0])
    P = P_from_model(J, Q)
    assert np.allclose(J.T @ P + P @ J, -Q, atol=1e-8)
    assert np.all(np.linalg.eigvalsh(P) > 0)          # score is positive away from home


def test_score_falls_inside_the_certified_zone():
    rng = np.random.default_rng(0)
    P = P_from_model(numerical_jacobian(HOME, P_BASE), np.diag([1.0, 30.0, 10.0]))
    c = largest_safe_level(HOME, P, P_BASE, rng)
    assert c > 0
    X = sample_zone(HOME, P, c, 20_000, rng)
    v = V(X, HOME, P)
    assert np.all(Vdot(X, HOME, P, P_BASE)[v > 1e-6 * c] < 0)


def test_reefs_in_the_certified_zone_recover():
    rng = np.random.default_rng(1)
    P = P_from_model(numerical_jacobian(HOME, P_BASE), np.diag([1.0, 30.0, 10.0]))
    c = largest_safe_level(HOME, P, P_BASE, rng)
    X = sample_zone(HOME, P, c, 400, rng)
    end = settle(X, P_BASE, years=1500)
    assert np.all(np.max(np.abs(end - HOME[:, None]), axis=0) < 0.02)
