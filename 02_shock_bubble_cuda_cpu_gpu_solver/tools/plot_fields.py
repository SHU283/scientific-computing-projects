#!/usr/bin/env python3
import os
from pathlib import Path
from plot_utils import load_csv, save_imshow

# Read rho_500.csv/p_500.csv in data_2nd/data_gpu_2nd

step = 500
base = Path(__file__).resolve().parents[1]
outdir = f"{base}/project_figs"
os.makedirs(outdir, exist_ok=True)

rho_cpu = load_csv(f"{base}/data/data_2nd/rho_{step}.csv")
rho_gpu = load_csv(f"{base}/data/data_gpu_2nd/rho_{step}.csv")

p_cpu   = load_csv(f"{base}/data/data_2nd/p_{step}.csv")
p_gpu   = load_csv(f"{base}/data/data_gpu_2nd/p_{step}.csv")

save_imshow(rho_cpu, f"CPU density (step {step})", f"{outdir}/rho_cpu_step{step}.png")
save_imshow(rho_gpu, f"GPU density (step {step})", f"{outdir}/rho_gpu_step{step}.png")

save_imshow(p_cpu,   f"CPU pressure (step {step})", f"{outdir}/p_cpu_step{step}.png")
save_imshow(p_gpu,   f"GPU pressure (step {step})", f"{outdir}/p_gpu_step{step}.png")

print("Wrote 4 field plots into", outdir)
