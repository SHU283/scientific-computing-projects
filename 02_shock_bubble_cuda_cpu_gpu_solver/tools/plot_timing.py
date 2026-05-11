#!/usr/bin/env python3
import os
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Read timing.csv

base = Path(__file__).resolve().parents[1]
others = f"{base}/data/others"
outdir = f"{base}/project_figs"
os.makedirs(outdir, exist_ok=True)

rows = {}
with open(f"{others}/timing.csv", "r", newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        rows[row["platform"]] = {
            k: float(row[k]) for k in ["bc_s", "dt_s", "step_s", "total_s"]
        }

cpu = rows["cpu"]
gpu = rows["gpu"]

labels = ["BC", "dt", "step"]
cpu_vals = np.array([cpu["bc_s"], cpu["dt_s"], cpu["step_s"]])
gpu_vals = np.array([gpu["bc_s"], gpu["dt_s"], gpu["step_s"]])

x = np.arange(2)

plt.figure(figsize=(7.5, 4.5))
bottom = np.zeros(2)

for k, lab in zip(["bc_s", "dt_s", "step_s"], labels):
    vals = np.array([cpu[k], gpu[k]])
    plt.bar(x, vals, bottom=bottom, label=lab)
    bottom += vals

plt.xticks(x, ["CPU", "GPU"])
plt.ylabel("Time (s)")
plt.title("Runtime breakdown (BC / dt / step)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{outdir}/timing_breakdown.png", dpi=250)
plt.close()

speedup_wall = cpu["total_s"] / gpu["total_s"]
speedup_step = cpu["step_s"] / gpu["step_s"]

with open(f"{others}/timing_summary.txt", "w") as f:
    f.write(f"Wall speedup = {speedup_wall:.3f}x\n")
    f.write(f"Step speedup = {speedup_step:.3f}x\n")

print("Wrote timing_breakdown.png to", outdir)
print("Wrote timing_summary.txt to", others)
