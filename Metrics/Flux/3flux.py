import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from pathlib import Path

Twist_dir     = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Twist/")
Bottom_bz_dir = Path("/project/cstr/Jaiman/sp26/Project1/Test4/VAPOR2/Merge/")
OUT_DIR       = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Flux2/")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Load TWIST files and extract timestep indices (10, 20, 30, ...)
Twist_files = sorted(Twist_dir.glob("TWIST_MAP_*"))
timesteps   = [int(f.name.split('_')[-1]) for f in Twist_files]

# Load only the BOTTOM_BZ and TIME files matching those timesteps
Bz_files    = [Bottom_bz_dir / f"BOTTOM_BZ.{t:03d}"  for t in timesteps]
Atime_files = [Bottom_bz_dir / f"TIME.{t:03d}"        for t in timesteps]

# Verify all expected files exist
for f in Bz_files + Atime_files:
    assert f.exists(), f"Missing file: {f}"

num_files = len(Twist_files)
print(f"Number of timesteps: {num_files}")
print(f"Timesteps: {timesteps}")

GRID_N  = 361
N_CELLS = GRID_N * GRID_N  # 130,321 rows after cleaning

dx = dy = 1 / 250
dS = dx * dy

Tw_cri = 1.0

def drop_every_nth_row(data, n):
    """Remove every nth row (1-indexed: n, 2n, 3n, ...) from a loaded array."""
    mask = (np.arange(len(data)) + 1) % n != 0
    return data[mask]

Tw_flux_values = []
Atime_values   = []

for i in range(num_files):
    # --- TWIST_MAP: 362x361 + 360 rows → drop every 362nd, then trim trailing 360 ---
    Twist_raw   = np.genfromtxt(Twist_files[i], invalid_raise=False)
    Twist_clean = drop_every_nth_row(Twist_raw, 362)
    Twist_data  = Twist_clean[:N_CELLS]
    assert len(Twist_data) == N_CELLS, \
        f"TWIST_MAP[{i}]: expected {N_CELLS} rows after cleaning, got {len(Twist_data)}"

    # --- BOTTOM_BZ: 362x361 rows → drop every 362nd blank row ---
    Bz_raw  = np.genfromtxt(Bz_files[i], invalid_raise=False)
    Bz_data = drop_every_nth_row(Bz_raw, 362)
    assert len(Bz_data) == N_CELLS, \
        f"BOTTOM_BZ[{i}]: expected {N_CELLS} rows after cleaning, got {len(Bz_data)}"

    Atime_data = np.loadtxt(Atime_files[i])

    Tw = Twist_data[:, 2]
    Bz = np.abs(Bz_data[:, 2])

    indices = Tw > Tw_cri
    Tw_flux = np.sum(Bz[indices] * dS)

    Tw_flux_values.append(Tw_flux)
    Atime_values.append(Atime_data)

# Save output
Output_file_Tw = OUT_DIR / "Abs_Twist_flux_vs_Time.txt"
np.savetxt(Output_file_Tw, np.column_stack((Atime_values, Tw_flux_values)), fmt="%.6f")
print(f"Saved: {Output_file_Tw}")

# --- Plot ---
data = np.genfromtxt(Output_file_Tw)
time = data[:, 0]
flux = data[:, 1]

plt.figure(figsize=(7, 5))
plt.ylim(0, max(flux) * 1.05)
plt.tick_params(axis='both', which='major', labelsize=11)
plt.plot(time, flux, marker='o', linestyle='-')
plt.xlabel("Alfven Time")
plt.ylabel("Twist Flux")
plt.tight_layout()
plt.savefig(OUT_DIR / "Toroidal_flux.png", dpi=300)
plt.show()
