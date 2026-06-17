"""
batch_run.py — runs all init conditions sequentially.
Usage: python batch_run.py
       python batch_run.py rings lattice vortex   (subset only)
"""

import subprocess
import sys
import time

ALL_INITS = [
    "rings",
    "lattice",
    "two_clusters",
    # "vert_lines",
    # "horiz_rows",
    # "two_rows",
    # "four_blocks",
    # "hex_lattice",
    "shear_flow",
    "four_lines",
    "out_velocity",
    "in_velocity",
    # "single_flock",
    "vortex",
    "anti_vortex",
]

# run subset if given as args, else run all
inits = sys.argv[1:] if len(sys.argv) > 1 else ALL_INITS

print(f"=== Batch run: {len(inits)} init condition(s) ===\n")

for idx, name in enumerate(inits):
    print(f"[{idx+1}/{len(inits)}] Starting: {name}")
    t0 = time.time()
    result = subprocess.run(["./main", name])
    elapsed = time.time() - t0
    status = "OK" if result.returncode == 0 else f"FAILED (code {result.returncode})"
    print(f"[{idx+1}/{len(inits)}] {name} done in {elapsed:.1f}s — {status}\n")

print("=== All runs complete ===")