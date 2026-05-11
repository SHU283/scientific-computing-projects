#!/usr/bin/env python3
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

BASE   = Path(__file__).resolve().parents[1]
OUTDIR = os.path.join(BASE, "project_figs")
os.makedirs(OUTDIR, exist_ok=True)

# Timing data (seconds)
cpu_bc   = 0.0076
cpu_dt   = 0.47
cpu_step = 11.84
cpu_total = cpu_bc + cpu_dt + cpu_step

gpu_bc   = 0.19
gpu_dt   = 0.23
gpu_step = 2.06
gpu_total = gpu_bc + gpu_dt + gpu_step

# Figure 1: total runtime
labels = ["CPU", "GPU"]
times  = [cpu_total, gpu_total]

plt.figure(figsize=(5, 4))
plt.bar(labels, times)
plt.ylabel("Total wall-clock time (s)")
plt.title("Total runtime comparison (CPU vs GPU)")

for i, t in enumerate(times):
    plt.text(i, t, f"{t:.2f}s", ha="center", va="bottom")

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "runtime_total_cpu_vs_gpu.png"), dpi=300)
plt.close()

# Figure 2: runtime breakdown
bc   = [cpu_bc,   gpu_bc]
dt   = [cpu_dt,   gpu_dt]
step = [cpu_step, gpu_step]
x = np.arange(len(labels))

plt.figure(figsize=(6, 4))
plt.bar(x, bc, label="BC")
plt.bar(x, dt, bottom=bc, label="dt")
plt.bar(x, step, bottom=np.array(bc) + np.array(dt), label="step")

plt.xticks(x, labels)
plt.ylabel("Time (s)")
plt.title("Runtime breakdown (BC / dt / step)")
plt.legend()

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "runtime_breakdown_cpu_gpu.png"), dpi=300)
plt.close()

print("Performance figures written to", OUTDIR)
