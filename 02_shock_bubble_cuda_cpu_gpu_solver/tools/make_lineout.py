#!/usr/bin/env python3
import csv, os
from pathlib import Path

# Read rho_530.csv/p_530.csv in data_1st/data_2nd/data_gpu_1st

BASE    = Path(__file__).resolve().parents[1]
DATADIR = os.path.join(BASE, "data")
OTHERS  = os.path.join(DATADIR, "others")

CPU1 = os.path.join(DATADIR, "data_1st")
CPU2 = os.path.join(DATADIR, "data_2nd")
GPU1 = os.path.join(DATADIR, "data_gpu_1st")

def load_csv(path):
    with open(path, newline='') as f:
        return [list(map(float, row)) for row in csv.reader(f)]

def save_dat(path, x, a, b, c):
    with open(path, "w") as f:
        for i in range(len(x)):
            f.write(f"{x[i]} {a[i]} {b[i]} {c[i]}\n")

step = 530

rho1 = load_csv(os.path.join(CPU1, f"rho_{step}.csv"))
p1   = load_csv(os.path.join(CPU1, f"p_{step}.csv"))

rho2 = load_csv(os.path.join(CPU2, f"rho_{step}.csv"))
p2   = load_csv(os.path.join(CPU2, f"p_{step}.csv"))

rhog = load_csv(os.path.join(GPU1, f"rho_{step}.csv"))
pg   = load_csv(os.path.join(GPU1, f"p_{step}.csv"))

ny = len(rho1)
nx = len(rho1[0])
jmid = ny // 2
x = list(range(nx))

os.makedirs(OTHERS, exist_ok=True)
save_dat(os.path.join(OTHERS, "line_rho_midY_step530.dat"), x, rho1[jmid], rho2[jmid], rhog[jmid])
save_dat(os.path.join(OTHERS, "line_p_midY_step530.dat"),   x, p1[jmid],   p2[jmid],   pg[jmid])

print("Saved:")
print(" ", os.path.join(OTHERS, "line_rho_midY_step530.dat"))
print(" ", os.path.join(OTHERS, "line_p_midY_step530.dat"))
print("midY index =", jmid, "ny =", ny)
