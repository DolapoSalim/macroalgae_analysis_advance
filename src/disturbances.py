"""Disturbance events (guide §16).

Instant shock (at the moment of the event):
    coral      C -> (1 - coral_loss) * C
    seaweed    M -> M + seaweed_arrival        (spores settle on freshly freed rock)
    herbivores H -> (1 - herbivore_loss) * H

Nutrient pulse (over time, bell-shaped):
    N(t) = N0 + nutrient_pulse * exp(-(t - pulse_peak)^2 / (2 * pulse_width^2))

All sizes are illustrative assumptions, not measurements.
"""

import numpy as np

from src.parameters import make_parameters

SCENARIOS = {
    #                   coral_loss  seaweed_arrival  herbivore_loss  nutrient_pulse
    "No event":         (0.00,      0.00,            0.00,           0.0),
    "Mild storm":       (0.20,      0.02,            0.10,           0.0),
    "Strong storm":     (0.60,      0.05,            0.30,           0.0),
    "Heatwave":         (0.40,      0.05,            0.10,           0.0),
    "Pollution event":  (0.10,      0.03,            0.10,           1.5),
    "Overfishing":      (0.00,      0.00,            0.60,           0.0),
    "Compound event":   (0.60,      0.08,            0.60,           1.5),
}
PULSE_PEAK = 1.0    # years after the event
PULSE_WIDTH = 0.5   # years


def apply_shock(x, scenario):
    coral_loss, seaweed_arrival, herbivore_loss, _ = SCENARIOS[scenario]
    C, M, H = x
    M_new = min(M + seaweed_arrival, 1 - (1 - coral_loss) * C)
    return np.array([(1 - coral_loss) * C, M_new, (1 - herbivore_loss) * H])


def conditions_over_time(scenario, **base):
    """Return p(t): parameters at time t, including any nutrient pulse."""
    pulse = SCENARIOS[scenario][3]
    p0 = make_parameters(**base)
    if pulse == 0:
        return lambda t: p0

    def p(t):
        extra = pulse * np.exp(-((t - PULSE_PEAK) ** 2) / (2 * PULSE_WIDTH ** 2))
        return {**p0, "N": p0["N"] + extra}
    return p
