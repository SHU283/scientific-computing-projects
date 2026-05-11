import re, glob, os
import numpy as np
import matplotlib.pyplot as plt

def parse_adv_time(path, is_mpi=False):
    txt = open(path, "r", encoding="utf-8", errors="ignore").read()
    if is_mpi:
        m = re.search(r"\[TIMER MPI\]\s+adv\s+total\s+=\s+([0-9.eE+-]+)", txt)
    else:
        m = re.search(r"\[TIMER\]\s+adv\s+total\s+=\s+([0-9.eE+-]+)", txt)
    return float(m.group(1)) if m else None

def parse_steps(path, is_mpi=False):
    txt = open(path, "r", encoding="utf-8", errors="ignore").read()
    if is_mpi:
        m = re.search(r"\[TIMER MPI\]\s+steps\s+=\s+(\d+)", txt)
    else:
        m = re.search(r"\[TIMER\]\s+steps\s+=\s+(\d+)", txt)
    return int(m.group(1)) if m else None

def load_omp():
    rows=[]
    for f in sorted(glob.glob("output/old_runs/timing_omp_t*.log")):
        m=re.search(r"_t(\d+)\.log$", f)
        if not m:
            continue
        t=int(m.group(1))
        adv=parse_adv_time(f, is_mpi=False)
        steps=parse_steps(f, is_mpi=False)
        if adv is None:
            print("missing adv in", f)
            continue
        rows.append((t, adv, steps))
    rows.sort()
    return rows

def load_mpi():
    rows=[]
    for f in sorted(glob.glob("output/old_runs/timing_mpi_np*.log")):
        m=re.search(r"_np(\d+)\.log$", f)
        if not m:
            continue
        p=int(m.group(1))
        adv=parse_adv_time(f, is_mpi=True)
        steps=parse_steps(f, is_mpi=True)
        if adv is None:
            print("missing adv in", f)
            continue
        rows.append((p, adv, steps))
    rows.sort()
    return rows

def print_table(name, rows):
    if not rows:
        print(f"\n{name}: no rows found.")
        return
    base_time = rows[0][1]
    print(f"\n{name} (baseline = {rows[0][0]}):")
    print("n, adv_time(s), speedup, efficiency, steps")
    for n, t, steps in rows:
        sp = base_time / t
        eff = sp / n
        print(f"{n:>2}, {t:>10.4f}, {sp:>7.3f}, {eff:>10.3f}, {steps}")

def plot_speedup_eff(name, rows, out_prefix):
    if not rows:
        return
    n = np.array([r[0] for r in rows], dtype=float)
    t = np.array([r[1] for r in rows], dtype=float)
    base = t[0]
    speedup = base / t
    eff = speedup / n

    plt.figure()
    plt.plot(n, speedup, marker="o")
    plt.xlabel("n")
    plt.ylabel("Speedup")
    plt.title(f"{name} speedup (adv time)")
    plt.tight_layout()
    plt.savefig(out_prefix + "_speedup.png", dpi=200)
    plt.close()

    plt.figure()
    plt.plot(n, eff, marker="o")
    plt.xlabel("n")
    plt.ylabel("Efficiency")
    plt.title(f"{name} efficiency (adv time)")
    plt.tight_layout()
    plt.savefig(out_prefix + "_efficiency.png", dpi=200)
    plt.close()

def main():
    os.makedirs("figs", exist_ok=True)

    omp_rows = load_omp()
    mpi_rows = load_mpi()

    print_table("OpenMP", omp_rows)
    print_table("MPI", mpi_rows)

    plot_speedup_eff("OpenMP", omp_rows, "figs/omp_oldruns")
    plot_speedup_eff("MPI", mpi_rows, "figs/mpi_oldruns")

    print("\nWrote: figs/omp_oldruns_speedup.png, figs/omp_oldruns_efficiency.png, figs/mpi_oldruns_speedup.png, figs/mpi_oldruns_efficiency.png")

if __name__ == "__main__":
    main()
