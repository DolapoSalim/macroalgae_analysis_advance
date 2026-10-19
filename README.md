# Coral reef tipping points

[![Explore the interactive website](https://img.shields.io/badge/Explore_the_interactive_website-0A7A87?style=for-the-badge)](https://dolaposalim.github.io/macroalgae_analysis_advance/)

### → [dolaposalim.github.io/macroalgae_analysis_advance](https://dolaposalim.github.io/macroalgae_analysis_advance/)

[![Tests](https://github.com/DolapoSalim/macroalgae_analysis_advance/actions/workflows/tests.yml/badge.svg)](https://github.com/DolapoSalim/macroalgae_analysis_advance/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23471354.svg)](https://doi.org/10.5281/zenodo.23471354)

Run the reef model in your browser, see the tipping point, and read the findings for everyday readers, scientists and decision-makers.

---

A small coral–seaweed–herbivore model, built from scratch with **synthetic data**, to test how scientists study reef resilience: tipping points, recovery zones, early-warning signals, parameter fitting and management scenarios. Because the data are made up, the right answers are known, so each method can be checked.

> The model and its numbers are illustrative. It is a practice reef for testing methods, not a forecast for any real reef.

## Contents

- `src/` – the model and analysis tools
- `steps/` – one script per analysis step (run in order, 01–14)
- `tests/` – automated checks
- `data/`, `figures/`, `reports/logs/` – generated data, figures and run outputs
- `website/` – the explainer site

## Run it

```bash
pip install -r requirements.txt
python -m pytest                          # checks (about 10 s)
python -m steps.step03_grazing_sweep      # example: find the tipping point
```

## Cite and license

Use the **"Cite this repository"** button on GitHub, or cite the archived release: Olatoye, D. (2026). *Coral reef tipping points: a synthetic coral–macroalgae–herbivore resilience model* (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23471355

Code: MIT. Website text and figures: CC BY 4.0.

Author: Dolapo Olatoye
