"""Step 5 — Map the basins: start reefs everywhere, colour each by the home it reaches.

Grazing is fixed at 0.19 (the two-homes zone from Step 3).

Outputs:
  figures/step05_basins_3d.html   interactive 3D map (open in a browser, drag to rotate)
  figures/step05_basins_3d.png    same map as a still picture
  figures/step05_basin_slice.png  a flat slice of the map, herbivores fixed at 0.5

Run from the project root:
    python -m steps.step05_basin_map
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from matplotlib.colors import ListedColormap

from src.basins import UNRESOLVED, grid_of_reefs, map_basins, seaweed_tolerance
from src.equilibria import find_equilibria
from src.parameters import make_parameters
from src.stability import classify

G = 0.19
HERBIVORE_LEVELS_3D = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
COLOURS = {0: "#1f77b4", 1: "#d62728", UNRESOLVED: "#999999"}  # coral home, seaweed home, unclear
OUT = Path("figures")


def main():
    p = make_parameters(G_max=G)
    points = [x for x in find_equilibria(p) if x[2] > 1e-6]
    homes = [x for x in points if classify(x, p)[1].startswith("STABLE")]   # [coral, seaweed]
    hilltop = [x for x in points if classify(x, p)[1].startswith("UNSTABLE")][0]

    # ---- 3D map -------------------------------------------------------------
    starts = grid_of_reefs(step=0.05, herbivore_levels=HERBIVORE_LEVELS_3D)
    print(f"Simulating {len(starts)} starting reefs for the 3D map...")
    ends = map_basins(starts, p, homes)

    print("\nShare of starting reefs that end in the CORAL home, by starting herbivores:")
    for H in HERBIVORE_LEVELS_3D:
        sel = np.isclose(starts[:, 2], H)
        print(f"  herbivores {H:.1f}:  {np.mean(ends[sel] == 0):5.0%} coral home, "
              f"{np.mean(ends[sel] == 1):5.0%} seaweed home, "
              f"{np.mean(ends[sel] == UNRESOLVED):4.0%} unclear")
    print(f"  ALL:              {np.mean(ends == 0):5.0%} coral home, "
          f"{np.mean(ends == 1):5.0%} seaweed home")

    # ---- 3D "tolerance wall": max seaweed a reef can start with and still recover
    coral_axis = np.round(np.arange(0.0, 0.951, 0.05), 3)
    herb_axis = np.round(np.arange(0.1, 1.51, 0.1), 3)
    print(f"\nFinding the tolerance wall on a {len(coral_axis)} x {len(herb_axis)} grid...")
    wall = np.array([[seaweed_tolerance(C, H, p, homes) for C in coral_axis] for H in herb_axis])
    for H, row in zip(herb_axis, wall):
        if round(H * 10) % 2 == 1:  # print H = 0.1, 0.3, ...
            print(f"  herbivores {H:.1f}: max seaweed for recovery at coral 20% / 50% / 80% = "
                  f"{row[4]:.2f} / {row[10]:.2f} / {row[16]:.2f}")

    # interactive 3D (plotly): the wall, the homes, and (hidden by default) every simulated reef
    fig = go.Figure()
    fig.add_trace(go.Surface(
        x=coral_axis, y=herb_axis, z=wall, colorscale="RdYlBu_r", opacity=0.9,
        colorbar=dict(title="max seaweed<br>for recovery"), name="tolerance wall"))
    names = {0: "ends in CORAL home", 1: "ends in SEAWEED home", UNRESOLVED: "unclear"}
    for k, name in names.items():
        sel = ends == k
        if sel.any():
            fig.add_trace(go.Scatter3d(
                x=starts[sel, 0], y=starts[sel, 2], z=starts[sel, 1], mode="markers",
                marker=dict(size=3, color=COLOURS[k], opacity=0.55), name=name,
                visible="legendonly"))
    for x, label, colour in [(homes[0], "coral home (valley)", "black"),
                             (homes[1], "seaweed home (valley)", "black"),
                             (hilltop, "balance point (mountain pass)", "white")]:
        fig.add_trace(go.Scatter3d(
            x=[x[0]], y=[x[2]], z=[x[1]], mode="markers+text", text=[label],
            marker=dict(size=7, color=colour, line=dict(color="black", width=3)), name=label))
    fig.update_layout(
        title=(f"Tolerance wall (grazing = {G}): reefs BELOW the wall recover to coral, "
               "reefs ABOVE it end as seaweed reefs"),
        scene=dict(xaxis_title="Coral cover (C)", yaxis_title="Herbivores (H)",
                   zaxis_title="Seaweed cover (M)"),
        legend=dict(x=0, y=1))
    OUT.mkdir(exist_ok=True)
    fig.write_html(OUT / "step05_basins_3d.html", include_plotlyjs="cdn")

    # still 3D picture (matplotlib)
    fig3 = plt.figure(figsize=(8, 6.5))
    ax = fig3.add_subplot(projection="3d")
    CC, HH = np.meshgrid(coral_axis, herb_axis)
    ax.plot_surface(CC, HH, wall, cmap="RdYlBu_r", edgecolor="grey", linewidth=0.3, alpha=0.9)
    for x, face in [(homes[0], "black"), (homes[1], "black"), (hilltop, "white")]:
        ax.scatter(x[0], x[2], x[1], s=120, color=face, edgecolor="black", linewidth=2,
                   depthshade=False)
    ax.set_xlabel("Coral cover (C)")
    ax.set_ylabel("Herbivores (H)")
    ax.set_zlabel("Max seaweed for recovery")
    ax.set_title(f"Tolerance wall (grazing = {G})\n"
                 "below the wall -> coral home, above it -> seaweed home")
    ax.view_init(elev=25, azim=-130)
    fig3.tight_layout()
    fig3.savefig(OUT / "step05_basins_3d.png", dpi=120)

    # ---- flat slice: herbivores start at 0.5 --------------------------------
    H_slice = 0.5
    step = 0.02
    slice_starts = grid_of_reefs(step=step, herbivore_levels=[H_slice])
    print(f"\nSimulating {len(slice_starts)} starting reefs for the flat slice (herbivores = {H_slice})...")
    slice_ends = map_basins(slice_starts, p, homes)

    n = int(round(1 / step)) + 1
    image = np.full((n, n), np.nan)
    for (C, M, _), k in zip(slice_starts, slice_ends):
        image[int(round(M / step)), int(round(C / step))] = {0: 0, 1: 1, UNRESOLVED: 2}[k]

    fig2, ax2 = plt.subplots(figsize=(6.5, 5.5))
    cmap = ListedColormap([COLOURS[0], COLOURS[1], COLOURS[UNRESOLVED]])
    ax2.imshow(image, origin="lower", extent=(-step / 2, 1 + step / 2, -step / 2, 1 + step / 2),
               cmap=cmap, vmin=0, vmax=2, alpha=0.8)
    ax2.plot([0, 1], [1, 0], color="black", lw=1)
    ax2.text(0.55, 0.52, "impossible\n(coral + seaweed > 100%)", fontsize=9, color="grey")
    ax2.plot([], [], "s", color=COLOURS[0], label="ends in CORAL home")
    ax2.plot([], [], "s", color=COLOURS[1], label="ends in SEAWEED home")
    ax2.set_xlabel("Starting coral cover (C)")
    ax2.set_ylabel("Starting seaweed cover (M)")
    ax2.set_title(f"Flat slice: starting herbivores = {H_slice}, grazing = {G}")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.legend(loc="upper right", fontsize=9)
    fig2.tight_layout()
    fig2.savefig(OUT / "step05_basin_slice.png", dpi=120)

    print(f"  slice: {np.mean(slice_ends == 0):.0%} coral home, {np.mean(slice_ends == 1):.0%} seaweed home")
    print(f"\nFigures saved in {OUT}/ (open step05_basins_3d.html in a browser to rotate)")


if __name__ == "__main__":
    main()
