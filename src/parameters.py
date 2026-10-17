"""Baseline parameter values for the coral–macroalgae–herbivore model.

Time unit: years. Cover variables (C, M) are proportions of the reef (0–1).
Herbivore biomass H is standardised so that 1.0 is a well-stocked reef.

These are illustrative starting values (mostly from Section 6 of the guide),
chosen to give interpretable simulations. They are NOT field estimates.
"""

BASELINE = {
    # --- Coral (C) ---
    "r_C": 0.20,    # coral lateral growth into free space, per year
    "r_rec": 0.05,  # coral larval recruitment onto free space, per year
    "mu_C": 0.05,   # coral background mortality, per year
    "a_CM": 0.35,   # rate at which macroalgae overgrow coral, per year

    # --- Macroalgae (M) ---
    "r_M": 0.45,    # macroalgal growth into free space, per year
    "beta_N": 0.20, # extra macroalgal growth per unit nutrient
    "mu_M": 0.12,   # macroalgal background loss (dislodgement, senescence), per year

    # --- Grazing (herbivores eating macroalgae) ---
    "G_max": 0.40,  # maximum grazing rate per unit herbivore, per year
    "K_G": 0.30,    # algal cover at which grazing is half its maximum

    # --- Herbivores (H) ---
    "r_H": 0.15,    # herbivore intrinsic growth rate, per year
    "K_H": 1.00,    # herbivore carrying capacity (standardised)
    "q_M": 0.10,    # herbivore growth boost from abundant food (not in guide; chosen here)
    "mu_H": 0.08,   # herbivore natural mortality, per year
    "F_H": 0.00,    # fishing mortality on herbivores, per year (0 = no fishing)

    # --- Nutrients (held constant for now) ---
    "N": 0.50,      # nutrient level (0.5 = low/moderate, 2.0 = heavily enriched)
}


def make_parameters(**overrides):
    """Return a copy of the baseline parameters, with any overrides applied.

    Example: make_parameters(G_max=0.1, N=1.2)
    """
    unknown = set(overrides) - set(BASELINE)
    if unknown:
        raise KeyError(f"Unknown parameter(s): {sorted(unknown)}")
    return {**BASELINE, **overrides}
