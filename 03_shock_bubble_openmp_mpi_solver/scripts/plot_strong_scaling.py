import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "output" / "scaling" / "strong"
FIGS = ROOT / "figs"
FIGS.mkdir(parents=True, exist_ok=True)

CORES = np.array([1, 2, 4, 8, 16], dtype=float)

def parse_per_step_adv(log_path: Path) -> float:
    if not log_path.exists():
        raise FileNotFoundError(f"Missing log: {log_path}")
    txt = log_path.read_text(errors="ignore")
    m = re.search(r"per-step adv\s*=\s*([0-9.eE+-]+)\s*s", txt)
    if not m:
        m = re.search(r"per-step\s+adv\s*=\s*([0-9.eE+-]+)\s*s", txt)
    if not m:
        raise RuntimeError(f"Cannot find 'per-step adv = ... s' in {log_path}")
    return float(m.group(1))

def speedup_eff(cores, times):
    t1 = times[0]
    sp = t1 / times
    eff = sp / cores
    return sp, eff

def plot_line(x, y, title, xlabel, ylabel, outpng):
    plt.figure()
    plt.plot(x, y, marker="o")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)
    plt.xticks(x)
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()

omp_times = np.array([parse_per_step_adv(OUT / f"run_omp_{int(c)}.log") for c in CORES])
omp_sp, omp_eff = speedup_eff(CORES, omp_times)

plot_line(CORES, omp_sp,  "OpenMP strong scaling: speedup",    "Threads",   "Speedup",
          FIGS / "omp_speedup.png")
plot_line(CORES, omp_eff, "OpenMP strong scaling: efficiency", "Threads",   "Efficiency",
          FIGS / "omp_efficiency.png")

mpi_times = np.array([parse_per_step_adv(OUT / f"run_mpi_{int(c)}.log") for c in CORES])
mpi_sp, mpi_eff = speedup_eff(CORES, mpi_times)

plot_line(CORES, mpi_sp,  "MPI strong scaling: speedup",       "Processes", "Speedup",
          FIGS / "mpi_speedup.png")
plot_line(CORES, mpi_eff, "MPI strong scaling: efficiency",    "Processes", "Efficiency",
          FIGS / "mpi_efficiency.png")

print("Done. Wrote:")
print("  ", FIGS / "omp_speedup.png")
print("  ", FIGS / "omp_efficiency.png")
print("  ", FIGS / "mpi_speedup.png")
print("  ", FIGS / "mpi_efficiency.png")
