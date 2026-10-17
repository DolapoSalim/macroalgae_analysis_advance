"""Step 8 — A 'safety certificate' (Lyapunov function) for the coral home.

Grazing = 0.19 (two-homes zone). Target = the coral home.

Three versions of the score V = (x - home)^T P (x - home):
  1. simple guess       P = identity (all three variables weighted equally)   [guide §12]
  2. from the model     P solves J^T P + P J = -I                             [guide §13.5]
  3. tuned              same, but with the weights Q searched for the biggest zone

For each: find the biggest zone where the score always goes down, then compare its
size with the REAL safe zone found by simulation (as in Step 5).

Run from the project root:
    python -m steps.step08_safety_certificate
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

from src.equilibria import find_equilibria
from src.lyapunov import P_from_model, V, largest_safe_level, sample_zone
from src.parameters import make_parameters
from src.stability import classify, numerical_jacobian
from src.stochastic import settle, which_home_many

G = 0.19
H_MAX = 1.5            # the map covers herbivores 0..1.5 (as in Step 5)
OUT = Path("figures")


def random_reefs(n, rng):
    """Reefs spread evenly over the map: C + M <= 1, 0 <= H <= H_MAX."""
    CM = rng.dirichlet([1, 1, 1], n)[:, :2].T
    return np.vstack([CM, rng.uniform(0, H_MAX, n)])


def main():
    rng = np.random.default_rng(8)
    p = make_parameters(G_max=G)
    pts = [x for x in find_equilibria(p) if x[2] > 1e-6]
    homes = [x for x in pts if classify(x, p)[1].startswith("STABLE")]   # [coral, seaweed]
    saddle = [x for x in pts if classify(x, p)[1].startswith("UNSTABLE")][0]
    home = homes[0]
    J = numerical_jacobian(home, p)

    # ---- the three scores
    candidates = {"1. simple guess": np.eye(3),
                  "2. from the model": P_from_model(J, np.eye(3))}
    print("Searching for the best weights for version 3...")
    best = None
    for q2 in np.logspace(0, 3, 7):
        for q3 in np.logspace(-1, 2, 7):
            P = P_from_model(J, np.diag([1.0, q2, q3]))
            c = largest_safe_level(home, P, p, rng)
            volume = c ** 1.5 / np.sqrt(np.linalg.det(P))   # ellipsoid size, up to a constant
            if best is None or volume > best[0]:
                best = (volume, P, (1.0, float(q2), float(q3)))
    candidates["3. tuned"] = best[1]
    print(f"  best weights Q = diag{tuple(round(q, 2) for q in best[2])}")

    # ---- size of the REAL safe zone (simulation, as in Step 5)
    reefs = random_reefs(20_000, rng)
    real_share = np.mean(which_home_many(settle(reefs, p), homes) == 0)

    # ---- size of each certified zone
    many = random_reefs(2_000_000, rng)
    levels = {}
    print(f"\nREAL safe zone (by simulation): {real_share:.1%} of the map\n")
    print(f"{'score':20s} {'biggest safe level c':>22s} {'share of map':>13s} {'share of REAL safe zone':>24s}")
    for name, P in candidates.items():
        c = largest_safe_level(home, P, p, rng)
        levels[name] = c
        share = np.mean(V(many, home, P) <= c) if c > 0 else 0.0
        print(f"{name:20s} {c:22.2e} {share:13.4%} {share / real_share:24.3%}")

    # ---- check the certificate: do reefs inside the tuned zone really recover?
    P3, c3 = candidates["3. tuned"], levels["3. tuned"]
    inside = sample_zone(home, P3, c3, 3000, rng)
    recovered = np.mean(which_home_many(settle(inside, p), homes) == 0)
    print(f"\nCheck: {inside.shape[1]} reefs inside the tuned zone -> {recovered:.1%} recovered to coral")

    # ---- Figure: flat slice at the home's herbivore level
    H0 = home[2]
    n = 101
    Cg, Mg = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    feasible = (Cg + Mg) <= 1
    grid = np.vstack([Cg[feasible], Mg[feasible], np.full(feasible.sum(), H0)])
    basin = np.full(Cg.shape, np.nan)
    basin[feasible] = (which_home_many(settle(grid, p), homes) == 0).astype(float)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, zoom in zip(axes, [False, True]):
        ax.contourf(Cg, Mg, basin, levels=[-0.5, 0.5, 1.5], colors=["#f4a6a6", "#a6c8f4"])
        for (name, P), colour in zip(list(candidates.items())[1:], ["#ff7f0e", "#2ca02c"]):
            grid_all = np.vstack([Cg.ravel(), Mg.ravel(), np.full(Cg.size, H0)])
            Vg = V(grid_all, home, P).reshape(Cg.shape)
            ax.contour(Cg, Mg, np.where(feasible, Vg, np.nan), levels=[levels[name]],
                       colors=colour, linewidths=2)
            ax.plot([], [], color=colour, lw=2, label=f"certified zone: {name[3:]}")
        ax.plot(*home[:2], "ko", ms=8, label="coral home")
        ax.plot(*saddle[:2], "o", mfc="white", mec="black", mew=2, ms=8, label="mountain pass")
        ax.set_xlabel("Coral cover (C)")
        ax.set_ylabel("Seaweed cover (M)")
        if zoom:
            ax.set_xlim(0.6, 0.95)
            ax.set_ylim(0, 0.15)
            ax.set_title("Zoomed in around the coral home")
        else:
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_title(f"Slice at herbivores = {H0:.2f}  (blue = really safe, red = tips)")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "step08_certificate_slice.png", dpi=120)

    # ---- interactive 3D: real safe zone (cloud) vs certified zones (solid shapes)
    safe_pts = reefs[:, which_home_many(settle(reefs, p), homes) == 0][:, :4000]
    fig3 = go.Figure()
    fig3.add_trace(go.Scatter3d(x=safe_pts[0], y=safe_pts[1], z=safe_pts[2], mode="markers",
                                marker=dict(size=2, color="#1f77b4", opacity=0.25),
                                name="REALLY safe reefs (simulation)"))
    for (name, P), colour in zip(list(candidates.items())[1:], ["#ff7f0e", "#2ca02c"]):
        shell = sample_zone(home, P, levels[name], 20_000, rng)
        fig3.add_trace(go.Mesh3d(x=shell[0], y=shell[1], z=shell[2], alphahull=0,
                                 color=colour, opacity=0.6, name=f"certified: {name[3:]}",
                                 showlegend=True))
    for x, label, colour in [(home, "coral home", "black"), (saddle, "mountain pass", "white")]:
        fig3.add_trace(go.Scatter3d(x=[x[0]], y=[x[1]], z=[x[2]], mode="markers+text",
                                    text=[label], name=label,
                                    marker=dict(size=6, color=colour, line=dict(color="black", width=3))))
    fig3.update_layout(title="Certified (guaranteed) safe zones vs the real safe zone - drag to rotate",
                       scene=dict(xaxis_title="Coral (C)", yaxis_title="Seaweed (M)",
                                  zaxis_title="Herbivores (H)"))
    fig3.write_html(OUT / "step08_certificate_3d.html", include_plotlyjs="cdn")
    print(f"\nFigures saved in {OUT}/")


if __name__ == "__main__":
    main()
