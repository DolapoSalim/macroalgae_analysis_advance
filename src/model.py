"""Coral–macroalgae–herbivore model.

State x = (C, M, H):
    C  coral cover        (proportion of reef, 0–1)
    M  macroalgal cover   (proportion of reef, 0–1)
    H  herbivore biomass  (standardised, ~1 = well stocked)

Nutrients N are held constant (a parameter) in this first version.

Space rule: coral and macroalgae can only grow into free space
B = 1 - C - M (bare rock / turf). This keeps C + M <= 1 automatically,
with no need to clip values during simulation.
"""

import numpy as np


def free_space(C, M):
    """Proportion of the reef not occupied by coral or macroalgae."""
    return 1.0 - C - M


def grazing(M, H, p):
    """Saturating grazing: herbivores get 'swamped' when algae is abundant.

    At low M, grazing ~ (G_max / K_G) * H * M  (each bit of algae gets eaten fast).
    At high M, grazing -> G_max * H            (herbivores are eating flat out).
    """
    return p["G_max"] * H * M / (p["K_G"] + M)


def process_rates(x, p):
    """Every individual process in the model, as a dictionary.

    Useful for seeing *why* a variable is increasing or decreasing.
    """
    C, M, H = x
    B = free_space(C, M)
    return {
        # coral
        "coral_growth": p["r_C"] * C * B,
        "coral_recruitment": p["r_rec"] * B,
        "coral_overgrown_by_algae": p["a_CM"] * M * C,
        "coral_mortality": p["mu_C"] * C,
        # macroalgae
        "algae_growth": (p["r_M"] + p["beta_N"] * p["N"]) * M * B,
        "algae_overgrowing_coral": p["a_CM"] * M * C,
        "algae_grazed": grazing(M, H, p),
        "algae_loss": p["mu_M"] * M,
        # herbivores
        "herbivore_growth": p["r_H"] * H * (1 - H / p["K_H"]),
        "herbivore_food_boost": p["q_M"] * M * H,
        "herbivore_mortality": p["mu_H"] * H,
        "herbivore_fishing": p["F_H"] * H,
    }


def rhs(t, x, p):
    """Right-hand side dx/dt = f(x). Signature (t, x, p) matches scipy's solve_ivp."""
    r = process_rates(x, p)
    dC = (r["coral_growth"] + r["coral_recruitment"]
          - r["coral_overgrown_by_algae"] - r["coral_mortality"])
    dM = (r["algae_growth"] + r["algae_overgrowing_coral"]
          - r["algae_grazed"] - r["algae_loss"])
    dH = (r["herbivore_growth"] + r["herbivore_food_boost"]
          - r["herbivore_mortality"] - r["herbivore_fishing"])
    return np.array([dC, dM, dH])
