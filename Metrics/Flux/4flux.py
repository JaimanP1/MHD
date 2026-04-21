import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

Twist_dir     = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Twist/")
Bottom_bz_dir = Path("/project/cstr/Jaiman/sp26/Project1/Test4/VAPOR2/Merge/")
OUT_DIR       = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Flux2/")
OUT_DIR.mkdir(parents=True, exist_ok=True)

Twist_files = sorted(Twist_dir.glob("TWIST_MAP_*"))
timesteps   = [int(f.name.split('_')[-1]) for f in Twist_files]

Bz_files    = [Bottom_bz_dir / f"BOTTOM_BZ.{t:03d}" for t in timesteps]
Atime_files = [Bottom_bz_dir / f"TIME.{t:03d}"       for t in timesteps]

for f in Bz_files + Atime_files:
    assert f.exists(), f"Missing file: {f}"

num_files = len(Twist_files)
print(f"Number of timesteps: {num_files}")
print(f"Timesteps: {timesteps}")

GRID_N  = 361
N_CELLS = GRID_N * GRID_N

dx = dy = 1 / 250
dS = dx * dy
Tw_cri = 1.0

def load_clean(filepath):
    """Load file and remove blank/NaN rows."""
    data = np.genfromtxt(filepath, invalid_raise=False)
    return data[~np.isnan(data).any(axis=1)]

Tw_flux_values = []
Atime_values   = []

for i in range(num_files):
    Twist_data = load_clean(Twist_files[i])[:N_CELLS]
    Bz_data    = load_clean(Bz_files[i])[:N_CELLS]

    assert len(Twist_data) == N_CELLS, \
        f"TWIST_MAP[{i}]: expected {N_CELLS} rows, got {len(Twist_data)}"
    assert len(Bz_data) == N_CELLS, \
        f"BOTTOM_BZ[{i}]: expected {N_CELLS} rows, got {len(Bz_data)}"

    Tw = Twist_data[:, 2]
    Bz = np.abs(Bz_data[:, 2])

    Tw_flux = np.sum(Bz[Tw > Tw_cri] * dS)

    Tw_flux_values.append(Tw_flux)
    Atime_values.append(np.loadtxt(Atime_files[i]))

# Save output
Output_file_Tw = OUT_DIR / "Abs_Twist_flux_vs_Time.txt"
np.savetxt(Output_file_Tw, np.column_stack((Atime_values, Tw_flux_values)), fmt="%.6f")
print(f"Saved: {Output_file_Tw}")

# Plot
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
