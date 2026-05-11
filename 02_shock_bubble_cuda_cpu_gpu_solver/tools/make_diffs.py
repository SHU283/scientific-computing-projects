#!/usr/bin/env python3
import csv, math, os
from pathlib import Path

# Read rho_530.csv/p_530.csv

BASE   = Path(__file__).resolve().parents[1]
DATADIR = os.path.join(BASE, "data")
OTHERS  = os.path.join(DATADIR, "others")

CPU1 = os.path.join(DATADIR, "data_1st")
GPU1 = os.path.join(DATADIR, "data_gpu_1st")

def load_csv(path):
    with open(path, newline='') as f:
        return [list(map(float, row)) for row in csv.reader(f)]

def save_csv(path, A):
    with open(path, "w", newline='') as f:
        w = csv.writer(f)
        w.writerows(A)

def l2(A):
    s = 0.0
    n = 0
    for r in A:
        for x in r:
            s += x * x
            n += 1
    return math.sqrt(s / max(n, 1))

step = 530

cpu_rho = load_csv(os.path.join(CPU1, f"rho_{step}.csv"))
gpu_rho = load_csv(os.path.join(GPU1, f"rho_{step}.csv"))
cpu_p   = load_csv(os.path.join(CPU1, f"p_{step}.csv"))
gpu_p   = load_csv(os.path.join(GPU1, f"p_{step}.csv"))

nr = len(cpu_rho)
nc = len(cpu_rho[0])

rho_diff = [[abs(cpu_rho[j][i] - gpu_rho[j][i]) for i in range(nc)] for j in range(nr)]
p_diff   = [[abs(cpu_p[j][i]   - gpu_p[j][i])   for i in range(nc)] for j in range(nr)]

rho_max = max(max(r) for r in rho_diff)
p_max   = max(max(r) for r in p_diff)

print(f"[DIFF] rho: max={rho_max:.3e}  L2={l2(rho_diff):.3e}")
print(f"[DIFF] p  : max={p_max:.3e}  L2={l2(p_diff):.3e}")

os.makedirs(OTHERS, exist_ok=True)
save_csv(os.path.join(OTHERS, f"rho_absdiff_{step}.csv"), rho_diff)
save_csv(os.path.join(OTHERS, f"p_absdiff_{step}.csv"), p_diff)
