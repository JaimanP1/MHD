#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

# -----------------------------
# USER SETTINGS
# -----------------------------
VMIN = 1e-6
VMAX = 1.0


DATA_DIR = Path("/project/cstr/Jaiman/sp26/Project1/Test3/Metrics/Twist/")        # directory containing TWIST_MAP_* files
OUT_DIR  = Path("/project/cstr/Jaiman/sp26/Project1/Test3/API/Outputs/Twist/")   # output directory
OUT_DIR.mkdir(exist_ok=True)

DPI = 200
FIGSIZE = (6, 5)
CMAP = "viridis"

TIME_STEPS = range(10, 81, 10)

# -----------------------------
# LOAD GRID FROM FIRST FILE
# -----------------------------
first_file = DATA_DIR / "TWIST_MAP_10"
data0 = np.fromfile(first_file, sep=' ').reshape(-1, 3)

x = data0[:, 0]
y = data0[:, 1]

x_unique = np.unique(x)
y_unique = np.unique(y)

NX = len(x_unique)
NY = len(y_unique)

print(f"Detected grid: {NX} × {NY}")

extent = [
    x_unique.min(), x_unique.max(),
    y_unique.min(), y_unique.max()
]

# -----------------------------
# MAIN LOOP
# -----------------------------
for step in TIME_STEPS:

    fname = DATA_DIR / f"TWIST_MAP_{step}"
    if not fname.exists():
        print(f"Skipping missing file: {fname}")
        continue

    data = np.fromfile(fname, sep=' ').reshape(-1, 3)

    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    z[z <= 0.0] = VMIN

    # Build Z safely
    Z = np.full((NY, NX), np.nan)

    # Map (x,y) → grid indices
    ix = np.searchsorted(x_unique, x)
    iy = np.searchsorted(y_unique, y)

    valid = (
        (ix >= 0) & (ix < NX) &
        (iy >= 0) & (iy < NY)
    )

    Z[iy[valid], ix[valid]] = z[valid]

    fig, ax = plt.subplots(figsize=FIGSIZE)

    im = ax.imshow(
        Z,
        origin="lower",
        extent=extent,
        norm=LogNorm(vmin=VMIN, vmax=VMAX),
        cmap=CMAP,
        aspect="auto"
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Value")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"TWIST MAP — t = {step}")

    plt.tight_layout()
    plt.savefig(OUT_DIR / f"twist_map_{step:03d}.png", dpi=DPI)
    plt.close(fig)

    print(f"Saved twist_map_{step:03d}.png")

print("Done.")
