#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from pathlib import Path

# -----------------------------
# USER SETTINGS
# -----------------------------
VMIN = 1e-6
VMAX = 1.0

# Directories to compare
DATA_DIR_A = Path("/project/cstr/Jaiman/sp26/Project1/Test3/Metrics/Twist/")
DATA_DIR_B = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Twist/")

OUT_DIR = Path("/project/cstr/Jaiman/sp26/Project1/Test5/API/Outputs/TwistComparison/")

TIME_STEPS = range(10, 81, 10)

DPI = 200
FIGSIZE = (6, 5)
CMAP_DIFF = "RdBu_r"
CMAP_BASE = "viridis"

# -----------------------------
# LOAD GRID (FROM A)
# -----------------------------
first_file = DATA_DIR_A / "TWIST_MAP_10"
data0 = np.fromfile(first_file, sep=" ").reshape(-1, 3)

x = data0[:, 0]
y = data0[:, 1]

x_unique = np.unique(x)
y_unique = np.unique(y)

NX = len(x_unique)
NY = len(y_unique)

extent = [
    x_unique.min(), x_unique.max(),
    y_unique.min(), y_unique.max()
]

print(f"Detected grid: {NX} × {NY}")

# -----------------------------
# HELPER FUNCTION
# -----------------------------
def load_twist(fname):
    data = np.fromfile(fname, sep=" ").reshape(-1, 3)
    x, y, z = data.T

    z[z <= 0.0] = VMIN

    Z = np.full((NY, NX), np.nan)

    ix = np.searchsorted(x_unique, x)
    iy = np.searchsorted(y_unique, y)

    valid = (
        (ix >= 0) & (ix < NX) &
        (iy >= 0) & (iy < NY)
    )

    Z[iy[valid], ix[valid]] = z[valid]
    return Z

# -----------------------------
# MAIN LOOP
# -----------------------------
for step in TIME_STEPS:

    fA = DATA_DIR_A / f"TWIST_MAP_{step}"
    fB = DATA_DIR_B / f"TWIST_MAP_{step}"

    if not (fA.exists() and fB.exists()):
        print(f"Skipping step {step}: missing file")
        continue

    ZA = load_twist(fA)
    ZB = load_twist(fB)

    # Difference
    Zdiff = ZA - ZB
    max_abs = np.nanmax(np.abs(Zdiff))

    # -------------------------
    # DIFFERENCE PLOT
    # -------------------------
    fig, ax = plt.subplots(figsize=FIGSIZE)

    im = ax.imshow(
        Zdiff,
        origin="lower",
        extent=extent,
        cmap=CMAP_DIFF,
        norm=TwoSlopeNorm(vcenter=0.0, vmin=-max_abs, vmax=max_abs),
        aspect="auto"
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Δ Twist (A − B)")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Twist Difference — t = {step}")

    plt.tight_layout()
    plt.savefig(OUT_DIR / f"twist_diff_{step:03d}.png", dpi=DPI)
    plt.close(fig)

    # -------------------------
    # OVERLAY PLOT
    # -------------------------
    fig, ax = plt.subplots(figsize=FIGSIZE)

    base = ax.imshow(
        ZB,
        origin="lower",
        extent=extent,
        norm=LogNorm(vmin=VMIN, vmax=VMAX),
        cmap=CMAP_BASE,
        aspect="auto"
    )

    cs = ax.contour(
        ZA,
        levels=8,
        colors="white",
        linewidths=0.8,
        extent=extent
    )

    ax.clabel(cs, fontsize=8, inline=True)

    cbar = fig.colorbar(base, ax=ax)
    cbar.set_label("Twist (B)")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Overlay: A (contours) on B — t = {step}")

    plt.tight_layout()
    plt.savefig(OUT_DIR / f"twist_overlay_{step:03d}.png", dpi=DPI)
    plt.close(fig)

    print(f"Saved comparison plots for t = {step}")

print("Comparison complete.")
