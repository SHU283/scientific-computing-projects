#!/usr/bin/env python3
import os
from pathlib import Path
from plot_utils import load_csv, save_imshow

# Read rho_diff_500.csv/p_diff_500.csv in others

step = 500
base = Path(__file__).resolve().parents[1]
others = f"{base}/data/others"
outdir = f"{base}/project_figs"
os.makedirs(outdir, exist_ok=True)

rho_diff = load_csv(f"{others}/rho_diff_{step}.csv")
p_diff   = load_csv(f"{others}/p_diff_{step}.csv")

save_imshow(rho_diff, f"|rho_CPU - rho_GPU| (step {step})",
            f"{outdir}/rho_diff_step{step}.png", log=True)

save_imshow(p_diff, f"|p_CPU - p_GPU| (step {step})",
            f"{outdir}/p_diff_step{step}.png", log=True)

print("Wrote 2 diff plots into", outdir)
