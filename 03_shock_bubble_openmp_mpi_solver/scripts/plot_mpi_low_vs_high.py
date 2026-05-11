import re
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FIGS = ROOT / "figs"

LOW_DIR  = ROOT / "output" / "scaling" / "strong"
HIGH_DIR = ROOT / "output" / "scaling" / "highres_mpi"

FIGS.mkdir(parents=True, exist_ok=True)

re_n_low  = re.compile(r"run_mpi_(\d+)\.log$")
re_n_high = re.compile(r"run_mpi_highres_(\d+)\.log$")
re_per_step = re.compile(r"\[TIMER MPI\]\s*per-step adv\s*=\s*([0-9.eE+-]+)\s*s")

def parse_logs(log_dir, re_n):
    data = []
    for f in sorted(log_dir.glob("run_mpi*.log")):
        m = re_n.search(f.name)
        if not m:
            continue
        N = int(m.group(1))
        txt = f.read_text(errors="ignore")
        m2 = re_per_step.search(txt)
        if not m2:
            raise RuntimeError(f"Cannot find per-step adv in {f}")
        t = float(m2.group(1))
        data.append((N, t))
    data.sort(key=lambda x: x[0])
    return data

low_data  = parse_logs(LOW_DIR,  re_n_low)
high_data = parse_logs(HIGH_DIR, re_n_high)

Ns_low  = [n for n, _ in low_data]
Ts_low  = [t for _, t in low_data]
Ns_high = [n for n, _ in high_data]
Ts_high = [t for _, t in high_data]

plt.figure()
plt.plot(Ns_low,  Ts_low,  marker="o", label="Low resolution (500×197)")
plt.plot(Ns_high, Ts_high, marker="o", label="High resolution")

plt.xscale("log", base=2)
plt.xlabel("MPI ranks N")
plt.ylabel("Runtime per advance step (s)")
plt.title("MPI strong scaling: low-res vs high-res")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(FIGS / "mpi_runtime_low_vs_high.png", dpi=200)
plt.close()

print("Low-res data:")
for n, t in low_data:
    print(f"  N={n:>2d}  per-step adv={t:.6g} s")

print("High-res data:")
for n, t in high_data:
    print(f"  N={n:>2d}  per-step adv={t:.6g} s")

print("Saved figure to:", FIGS / "mpi_runtime_low_vs_high.png")
