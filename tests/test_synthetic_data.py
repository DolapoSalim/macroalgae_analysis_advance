"""The synthetic-data generator records the storm correctly and simulate() accepts any duration."""

import numpy as np

from src.disturbances import apply_shock
from src.parameters import make_parameters
from src.simulate import simulate
from src.synthetic_data import simulate_truth, true_at


def test_survey_just_after_storm_sees_damaged_reef():
    """Regression test: a survey a few days after the storm must not blend the
    pre-storm and post-storm reefs."""
    p = make_parameters()
    t, X = simulate_truth((0.6, 0.1, 0.5), p, years=12, storm_year=8, sigma=0.0,
                          rng=np.random.default_rng(0))
    before = true_at(np.array([8.0]), t, X)[0]
    after = true_at(np.array([8.01]), t, X)[0]
    expected = apply_shock(before, "Strong storm")
    assert np.allclose(after, expected, atol=2e-3)


def test_simulate_accepts_fractional_years():
    t, X = simulate((0.5, 0.1, 0.5), make_parameters(), years=2.5, points_per_year=4)
    assert np.isclose(t[-1], 2.5) and len(t) == 11
