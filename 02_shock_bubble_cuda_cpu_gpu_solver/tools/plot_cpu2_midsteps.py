#!/usr/bin/env python3
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Read rho_150.csv/p_150.csv/rho_300.csv/p_300.csv in data_2nd

BASE    = Path(__file__).resolve().parents[1]
DATADIR = os.path.join(BASE, "data")
OUTDIR  = os.path.join(BASE, "project_figs")

CPU2 = os.path.join(DATADIR, "data_2nd")
STEPS = [150, 300]
DPI = 450

def load_csv(folder, name, step):
    path = os.path.join(folder, f"{name}_{step}.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")
    return np.loadtxt(path, delimiter=",")

def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close()

def heatmap(Z, title, cbar_label):
    plt.figure(figsize=(10, 4))
    plt.imshow(Z, origin="lower", aspect="auto")
    plt.title(title)
    plt.xlabel("x index")
    plt.ylabel("y index")
    cb = plt.colorbar()
    cb.set_label(cbar_label)

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("Output directory:", OUTDIR)

    for step in STEPS:
        print(f"Plotting CPU 2nd order, step = {step}")

        rho = load_csv(CPU2, "rho", step)
        p   = load_csv(CPU2, "p", step)

        heatmap(rho, f"CPU 2nd density rho (step {step})", "density ρ")
        savefig(os.path.join(OUTDIR, f"rho_cpu2_step{step}.png"))

        heatmap(p, f"CPU 2nd pressure p (step {step})", "pressure p")
        savefig(os.path.join(OUTDIR, f"p_cpu2_step{step}.png"))

    print("Done. 4 figures generated.")

if __name__ == "__main__":
    main()
