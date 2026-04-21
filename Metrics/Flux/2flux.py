import numpy as np 
import matplotlib.pyplot as plt 
import glob 
import os 


Twist_dir   = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Twist/")
#RUN_DIR   = Path("/project/cstr/Jaiman/sp26/Project1/Test3/Metrics/Twist/")
Bottom_bz_dir = Path("/project/cstr/Jaiman/sp26/Project1/Test4/VAPOR2/Merge/")
#VAPOR_DIR = Path("/project/cstr/Jaiman/fa25/Project1/Test1/VAPOR/Merge/")
OUT_DIR   = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Flux2/")
#OUT_DIR   = Path("/project/cstr/Jaiman/sp26/Project1/Test3/Metrics/Flux/")
Twist_files = sorted(glob.glob(Twist_dir + "TWIST_MAP.*"))
Bz_files    = sorted(glob.glob(Bottom_bz_dir + "BOTTOM_BZ.*"))
#Jz_files    = sorted(glob.glob(file_path + "BOTTOM_CZ.*"))
Atime_files = sorted(glob.glob(file_path + "TIME.*"))

if len(Twist_files) == len(Bz_files) == len(Atime_files):
    num_files = len(Twist_files)
    print(f"Number of timesteps: {num_files}")

GRID_N  = 361
N_CELLS = GRID_N * GRID_N  # 130,321 expected rows after cleaning

dx = dy = 1 / 250
dS = dx * dy

Tw_cri = 1.0

Tw_flux_values = []
#Jz_int_values  = []
Atime_values   = []

def drop_every_nth_row(data, n):
    """Remove every nth row (0-indexed: n-1, 2n-1, ...) from a loaded array."""
    mask = (np.arange(len(data)) + 1) % n != 0
    return data[mask]

for i in range(num_files):
    # --- TWIST_MAP: 362x361 + 360 rows ---
    # Step 1: read with genfromtxt so blank rows become NaN rows
    Twist_raw  = np.genfromtxt(Twist_files[i], invalid_raise=False)
    # Step 2: drop every 362nd row (the blank separators)
    Twist_clean = drop_every_nth_row(Twist_raw, 362)
    # Step 3: drop the final 360 trailing rows
    Twist_data  = Twist_clean[:N_CELLS]
    assert len(Twist_data) == N_CELLS, \
        f"TWIST_MAP[{i}]: expected {N_CELLS} rows after cleaning, got {len(Twist_data)}"

    # --- BOTTOM_BZ / BOTTOM_CZ: 362x361 rows ---
    # Step 1: read with genfromtxt so blank rows become NaN rows
    Bz_raw  = np.genfromtxt(Bz_files[i], invalid_raise=False)
#    Jz_raw  = np.genfromtxt(Jz_files[i], invalid_raise=False)
    # Step 2: drop every 362nd row (the blank separators)
    Bz_data = drop_every_nth_row(Bz_raw, 362)
#    Jz_data = drop_every_nth_row(Jz_raw, 362)
    assert len(Bz_data) == N_CELLS, \
        f"BOTTOM_BZ[{i}]: expected {N_CELLS} rows after cleaning, got {len(Bz_data)}"

    Atime_data = np.loadtxt(Atime_files[i])

    x, y, Tw = Twist_data[:, 0], Twist_data[:, 1], Twist_data[:, 2]
    Bz = Bz_data[:, 2]
#    Jz = Jz_data[:, 2]

    indices = Tw > Tw_cri

    Bz = np.abs(Bz)
#    Jz = np.abs(Jz)

    Tw_flux = np.sum(Bz[indices] * dS)
#    Jz_int  = np.sum(Jz[indices] * dS)

    Tw_flux_values.append(Tw_flux)
#    Jz_int_values.append(Jz_int)
    Atime_values.append(Atime_data)

# Save outputs
Output_file_Tw = os.path.join(file_path, "Abs_Twist_flux_vs_Time.txt")
#Output_file_Jz = os.path.join(file_path, "Abs_Current_vs_Time.txt")

np.savetxt(Output_file_Tw, np.column_stack((Atime_values, Tw_flux_values)), fmt="%.6f")
#np.savetxt(Output_file_Jz, np.column_stack((Atime_values, Jz_int_values)),  fmt="%.6f")

print("Successfully computed twist flux and current intensity for all timesteps.")

# --- Plot ---
Twist_file   = file_path + "Abs_Twist_flux_vs_Time.txt"
#Current_file = file_path + "Abs_Current_vs_Time.txt"

Twist   = np.genfromtxt(Twist_file)
#Current = np.genfromtxt(Current_file)

time = Twist[:, 0]
flux = Twist[:, 1]
#cur  = Current[:, 1]

plt.figure(figsize=(7, 5))
plt.ylim(0, max(flux) * 1.05)
plt.tick_params(axis='both', which='major', labelsize=11)
plt.plot(time, flux, marker='o', linestyle='-')
plt.tight_layout()
plt.savefig(file_path + "Toroidal_flux.png", dpi=300)
plt.show()
