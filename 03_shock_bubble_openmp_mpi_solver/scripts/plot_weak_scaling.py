import re
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

BASE = "output/scaling/weak"
OUTDIR = "figs"
os.makedirs(OUTDIR, exist_ok=True)

RE_OMP = re.compile(r"\[TIMER\]\s+per-step adv\s*=\s*([0-9.eE+-]+)\s*s")
RE_MPI = re.compile(r"\[TIMER MPI\]\s+per-step adv\s*=\s*([0-9.eE+-]+)\s*s")

def parse_logs(pattern_glob, regex):
    data = []
    for path in sorted(glob.glob(pattern_glob)):
        mN = re.search(r"_weak_(\d+)\.log$", path)
        if not mN:
            continue
        N = int(mN.group(1))
        with open(path, "r", errors="ignore") as f:
            txt = f.read()
        m = regex.search(txt)
        if not m:
            raise RuntimeError(f"Cannot find per-step adv in {path}")
        t = float(m.group(1))
        data.append((N, t))
    data.sort(key=lambda x: x[0])
    return np.array([x[0] for x in data]), np.array([x[1] for x in data])

N_omp, T_omp = parse_logs(os.path.join(BASE, "run_omp_weak_*.log"), RE_OMP)
N_mpi, T_mpi = parse_logs(os.path.join(BASE, "run_mpi_weak_*.log"), RE_MPI)

def efficiency(N, T):
    idx = np.where(N == 1)[0]
    if len(idx) != 1:
        raise RuntimeError("Need exactly one N=1 result to compute efficiency.")
    T1 = T[idx[0]]
    return T1 / T

E_omp = efficiency(N_omp, T_omp)
E_mpi = efficiency(N_mpi, T_mpi)

plt.figure()
plt.plot(N_omp, T_omp, marker="o", label="OpenMP")
plt.plot(N_mpi, T_mpi, marker="o", label="MPI")
plt.xscale("log", base=2)
plt.xlabel("N (threads or MPI ranks)")
plt.ylabel("Runtime per advance step (s)")
plt.title("Weak scaling: runtime per step vs N")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "weak_runtime_per_step.png"), dpi=200)
plt.close()

plt.figure()
plt.plot(N_omp, E_omp, marker="o", label="OpenMP")
plt.plot(N_mpi, E_mpi, marker="o", label="MPI")
plt.xscale("log", base=2)
plt.xlabel("N (threads or MPI ranks)")
plt.ylabel("Weak scaling efficiency  T(1)/T(N)")
plt.title("Weak scaling: efficiency vs N")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "weak_efficiency.png"), dpi=200)
plt.close()

print("Wrote:")
print("  figs/weak_runtime_per_step.png")
print("  figs/weak_efficiency.png")
