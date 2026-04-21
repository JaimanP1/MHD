#!/usr/bin/env python3
import numpy as np
from pathlib import Path

# -----------------------------
# USER SETTINGS
# -----------------------------
RUN_DIR   = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Twist/")
#RUN_DIR   = Path("/project/cstr/Jaiman/sp26/Project1/Test3/Metrics/Twist/")
VAPOR_DIR = Path("/project/cstr/Jaiman/sp26/Project1/Test4/VAPOR2/Merge/")
#VAPOR_DIR = Path("/project/cstr/Jaiman/fa25/Project1/Test1/VAPOR/Merge/")
OUT_DIR   = Path("/project/cstr/Jaiman/sp26/Project1/Test4/Metrics/Flux/")
#OUT_DIR   = Path("/project/cstr/Jaiman/sp26/Project1/Test3/Metrics/Flux/")

TIME_STEPS = range(10, 81, 10)
TW_THRESHOLD = 0.1
DX = None
DY = None

# Known grid dimensions — set explicitly to avoid float-imprecision issues
# with np.unique across the two files.
# Clean grid (Bz): NX_CLEAN x NY_CLEAN = 361 x 361 = 130321 points
# Padded grid (Twist): same NX but NY_PAD rows per x-block
NX        = 361   # number of x-values (shared by both grids)
NY_CLEAN  = 361   # y-rows per x-block in the clean (Bz) grid
# NY_PAD is inferred below from the actual file size

# -----------------------------
# LOAD FILE
# -----------------------------
def load_map(fname):
    data = np.fromfile(fname, sep=" ").reshape(-1, 3)
    return data[:, 0], data[:, 1], data[:, 2]

# -----------------------------
# STRIP BOUNDARY ROWS
#
# The padded grid has NX blocks, each with NY_PAD y-rows.
# We drop either the first or last (NY_PAD - NY_CLEAN) rows per block.
# strip_end=False  -> drop from start of each block (boundary at top)
# strip_end=True   -> drop from end  of each block (boundary at bottom)
# -----------------------------
def strip_by_block(vals, nx, ny_pad, ny_clean, strip_end=False):
    n_extra = ny_pad - ny_clean
    assert len(vals) == nx * ny_pad, (
        f"Array length {len(vals)} != nx*ny_pad = {nx*ny_pad}"
    )
    idx = []
    for i in range(nx):
        base = i * ny_pad
        if strip_end:
            idx.extend(range(base, base + ny_clean))
        else:
            idx.extend(range(base + n_extra, base + ny_pad))
    return vals[np.array(idx)]

# -----------------------------
# STRIP WITH REMAINDER HANDLING
#
# 361 * 363 = 131043 but actual count is 131042 — one x-block has one
# fewer boundary row. We pad a dummy row to make the math exact, strip,
# then remove the dummy at the end.
# -----------------------------
def strip_padded(x_pad, y_pad, val_pad, nx, ny_pad, ny_clean, strip_end=False):
    total_needed = nx * ny_pad
    shortfall    = total_needed - len(val_pad)

    if shortfall > 0:
        x_pad   = np.concatenate([x_pad,   np.full(shortfall, x_pad[-1])])
        y_pad   = np.concatenate([y_pad,   np.full(shortfall, y_pad[-1])])
        val_pad = np.concatenate([val_pad, np.full(shortfall, np.nan)])

    x_out   = strip_by_block(x_pad,   nx, ny_pad, ny_clean, strip_end)
    y_out   = strip_by_block(y_pad,   nx, ny_pad, ny_clean, strip_end)
    val_out = strip_by_block(val_pad, nx, ny_pad, ny_clean, strip_end)

    good = ~np.isnan(val_out)
    return x_out[good], y_out[good], val_out[good]

# -----------------------------
# INFER GRID DIMENSIONS FROM FILE SIZE
# -----------------------------
print("Detecting grid structure...")

x_tw0, y_tw0, _ = load_map(RUN_DIR   / "TWIST_MAP_10")
x_bz0, y_bz0, _ = load_map(VAPOR_DIR / "BOTTOM_BZ.010")

n_tw = len(x_tw0)
n_bz = len(x_bz0)

print(f"Twist raw points : {n_tw}")
print(f"Bz    raw points : {n_bz}")

assert n_bz == NX * NY_CLEAN, (
    f"Bz has {n_bz} points but NX*NY_CLEAN = {NX*NY_CLEAN}. "
    f"Update NX or NY_CLEAN at the top of the script."
)
assert n_tw > n_bz, "Expected Twist to be the padded grid."

NY_PAD    = int(np.round(n_tw / NX))   # 363
remainder = n_tw - NX * NY_PAD         # -1 (one block is 1 row shorter)

print(f"\nClean grid : {NX} x {NY_CLEAN} = {NX*NY_CLEAN}")
print(f"Padded rows per block : {NY_PAD}  (raw total {n_tw}, remainder {remainder})")
print(f"Extra rows per block  : {NY_PAD - NY_CLEAN}")

# -----------------------------
# AUTO-DETECT strip direction
# -----------------------------
print("\nAuto-detecting boundary-row position (start vs end of block)...")

_, y_s, _ = strip_padded(x_tw0, y_tw0, np.zeros(n_tw), NX, NY_PAD, NY_CLEAN, strip_end=False)
_, y_e, _ = strip_padded(x_tw0, y_tw0, np.zeros(n_tw), NX, NY_PAD, NY_CLEAN, strip_end=True)

# Compare the first block's y-values against the Bz first block
y_bz_b0 = np.sort(y_bz0[:NY_CLEAN])
diff_start = np.mean(np.abs(np.sort(y_s[:NY_CLEAN]) - y_bz_b0))
diff_end   = np.mean(np.abs(np.sort(y_e[:NY_CLEAN]) - y_bz_b0))

print(f"  mean |delta_y| strip-from-start : {diff_start:.3e}")
print(f"  mean |delta_y| strip-from-end   : {diff_end:.3e}")

STRIP_END = diff_end < diff_start
print(f"  -> Stripping from {'end' if STRIP_END else 'start'} of each x-block")

# Final verification
_, _, tw_check = strip_padded(x_tw0, y_tw0, np.zeros(n_tw), NX, NY_PAD, NY_CLEAN, STRIP_END)
print(f"  Stripped Twist size: {len(tw_check)}  (target {n_bz})")
assert len(tw_check) == n_bz, (
    f"Strip alignment failed: got {len(tw_check)}, expected {n_bz}. "
    "Check NX / NY_CLEAN / NY_PAD."
)

# -----------------------------
# GRID SPACING from clean Bz grid
# -----------------------------
xu_bz = np.unique(x_bz0)
yu_bz = np.unique(y_bz0)
dx    = DX if DX is not None else float(np.min(np.diff(np.sort(xu_bz))))
dy    = DY if DY is not None else float(np.min(np.diff(np.sort(yu_bz))))
dA    = dx * dy
print(f"\ndx = {dx:.6e},  dy = {dy:.6e},  dA = {dA:.6e}")

# -----------------------------
# MAIN LOOP
# -----------------------------
results = []

for step in TIME_STEPS:
    twist_file = RUN_DIR   / f"TWIST_MAP_{step}"
    bz_file    = VAPOR_DIR / f"BOTTOM_BZ.{step:03d}"

    if not (twist_file.exists() and bz_file.exists()):
        print(f"\nSkipping t={step}: missing file(s)")
        continue

    print(f"\n--- t = {step} ---")

    x_tw, y_tw, Tw = load_map(twist_file)
    _,    _,    Bz = load_map(bz_file)

    _, _, Tw = strip_padded(x_tw, y_tw, Tw, NX, NY_PAD, NY_CLEAN, STRIP_END)

    if len(Tw) != len(Bz):
        print(f"  WARNING: sizes still differ ({len(Tw)} vs {len(Bz)}) — skipping.")
        continue

    print(f"  Aligned size : {len(Tw)}")

    mask        = Tw > TW_THRESHOLD
    F_tw_signed = float(np.sum(Bz[mask])) * dA
    F_tw_abs    = float(np.sum(np.abs(Bz[mask]))) * dA
    results.append((step, F_tw_signed, F_tw_abs))

    print(f"  F_Tw (signed) = {F_tw_signed:.4e}")
    print(f"  |F_Tw|        = {F_tw_abs:.4e}")

# -----------------------------
# SAVE RESULTS
# -----------------------------
if results:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_file = OUT_DIR / f"flux_tw{TW_THRESHOLD}.txt"
    np.savetxt(out_file, results, fmt="%.6e",
               header="t  F_Tw_signed  F_Tw_abs")
    print(f"\nSaved -> {out_file}")
else:
    print("\nNo results to save.")
