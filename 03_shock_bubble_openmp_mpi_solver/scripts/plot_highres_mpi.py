import re
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FIGS = ROOT / "figs"
LOG_DIR = ROOT / "output" / "scaling" / "highres_mpi"

FIGS.mkdir(parents=True, exist_ok=True)

# Match: run_mpi_highres_16.log
re_n = re.compile(r"run_mpi_highres_(\d+)\.log$")
# Match: [TIMER MPI] per-step adv = 0.0231575 s
re_per_step = re.compile(r"\[TIMER MPI\]\s*per-step adv\s*=\s*([0-9.eE+-]+)\s*s")

data = []  # list of (N, per_step_adv_s)

for f in sorted(LOG_DIR.glob("run_mpi_highres_*.log")):
    m = re_n.search(f.name)
    if not m:
        continue
    N = int(m.group(1))

    txt = f.read_text(errors="ignore")
    m2 = re_per_step.search(txt)
    if not m2:
        raise RuntimeError(f"Cannot find '[TIMER MPI] per-step adv' in {f}")
    per_step = float(m2.group(1))
    data.append((N, per_step))

if not data:
    raise RuntimeError(f"No logs found in {LOG_DIR}")

data.sort(key=lambda x: x[0])
Ns = [n for n, _ in data]
Ts = [t for _, t in data]

# Strong scaling metrics (using per-step runtime as proxy)
T1 = Ts[Ns.index(1)] if 1 in Ns else Ts[0]
speedup = [T1 / t for t in Ts]
eff = [s / n for s, n in zip(speedup, Ns)]

# Plot 1: runtime per step vs N
plt.figure()
plt.plot(Ns, Ts, marker="o")
plt.xscale("log", base=2)
plt.xlabel("MPI ranks N")
plt.ylabel("Runtime per advance step (s)")
plt.title("High-resolution MPI strong scaling: runtime per step vs N")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(FIGS / "highres_mpi_runtime_per_step.png", dpi=200)
plt.close()

# Plot 2: speedup vs N
plt.figure()
plt.plot(Ns, speedup, marker="o")
plt.xscale("log", base=2)
plt.xlabel("MPI ranks N")
plt.ylabel("Speedup (T1 / TN)")
plt.title("High-resolution MPI strong scaling: speedup vs N")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(FIGS / "highres_mpi_speedup.png", dpi=200)
plt.close()

# Plot 3: efficiency vs N 
plt.figure()
plt.plot(Ns, eff, marker="o")
plt.xscale("log", base=2)
plt.xlabel("MPI ranks N")
plt.ylabel("Parallel efficiency (speedup / N)")
plt.title("High-resolution MPI strong scaling: efficiency vs N")
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(FIGS / "highres_mpi_efficiency.png", dpi=200)
plt.close()

print("Parsed high-res logs:")
for n, t in data:
    print(f"  N={n:>2d}  per-step adv={t:.6g} s")
print("Saved figures to:", FIGS)