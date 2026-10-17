"""The model's rules and resting points are internally consistent."""

import numpy as np

from src.equilibria import find_equilibria
from src.model import process_rates, rhs
from src.parameters import make_parameters
from src.sensitivity import coral_home, tipping_grazing
from src.stability import classify, numerical_jacobian


def test_rates_add_up():
    """dC/dt, dM/dt, dH/dt equal the sum of their individual processes."""
    p = make_parameters()
    x = (0.4, 0.3, 0.6)
    r = process_rates(x, p)
    dC, dM, dH = rhs(0, x, p)
    assert np.isclose(dC, r["coral_growth"] + r["coral_recruitment"]
                      - r["coral_overgrown_by_algae"] - r["coral_mortality"])
    assert np.isclose(dM, r["algae_growth"] + r["algae_overgrowing_coral"]
                      - r["algae_grazed"] - r["algae_loss"])
    assert np.isclose(dH, r["herbivore_growth"] + r["herbivore_food_boost"]
                      - r["herbivore_mortality"] - r["herbivore_fishing"])


def test_resting_points_really_rest():
    """Every resting point found has (almost) zero rates of change."""
    p = make_parameters(G_max=0.19)
    for x in find_equilibria(p):
        assert np.max(np.abs(rhs(0, x, p))) < 1e-9


def test_coral_home_formula_matches_solver():
    """The exact coral-home formula agrees with the numerical equation solver."""
    p = make_parameters()
    C, H = coral_home(p)
    found = [x for x in find_equilibria(p) if x[2] > 1e-6 and x[1] < 1e-9]
    assert any(np.allclose(x, [C, 0, H], atol=1e-6) for x in found)


def test_tipping_formula_matches_stability_test():
    """Just above the tipping point the coral home is stable, just below it is not."""
    g = float(tipping_grazing(make_parameters()))
    for G, expected in [(g + 0.005, "STABLE"), (g - 0.005, "UNSTABLE")]:
        p = make_parameters(G_max=G)
        C, H = coral_home(p)
        assert classify(np.array([C, 0.0, H]), p)[1].startswith(expected)


def test_jacobian_does_not_depend_on_step_size():
    p = make_parameters(G_max=0.19)
    x = np.array([0.62, 0.09, 0.53])
    assert np.allclose(numerical_jacobian(x, p, h=1e-5), numerical_jacobian(x, p, h=1e-7), atol=1e-5)


def test_bistability_window():
    """Two stable homes at grazing 0.19 (Step 3/4), one at the baseline 0.40."""
    for G, n_homes in [(0.19, 2), (0.40, 1)]:
        p = make_parameters(G_max=G)
        homes = [x for x in find_equilibria(p) if x[2] > 1e-6 and classify(x, p)[1].startswith("STABLE")]
        assert len(homes) == n_homes
