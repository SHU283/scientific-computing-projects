import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs("figs", exist_ok=True)

def plot_field(csv_path, out_png, title):
    a = np.loadtxt(csv_path, delimiter=",")

    plt.figure(figsize=(7,4))

    plt.imshow(
        a,
        origin="lower",
        aspect="auto",
        cmap="viridis",
        interpolation="bilinear"
    )

    plt.colorbar()
    plt.title(title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


todo = [
    # old snapshots
    ("output/old_runs/rho_0.csv",    "figs/rho_0.png",    "rho (t=0)"),
    ("output/old_runs/rho_400.csv",  "figs/rho_400.png",  "rho (t=400)"),
    ("output/old_runs/rho_800.csv",  "figs/rho_800.png",  "rho (t=800)"),
    ("output/old_runs/rho_1200.csv", "figs/rho_1200.png", "rho (t=1200)"),
    ("output/old_runs/rho_1600.csv", "figs/rho_1600.png", "rho (t=1600)"),
    ("output/old_runs/rho_2000.csv", "figs/rho_2000.png", "rho (t=2000)"),

    ("output/old_runs/p_0.csv",    "figs/p_0.png",    "p (t=0)"),
    ("output/old_runs/p_400.csv",  "figs/p_400.png",  "p (t=400)"),
    ("output/old_runs/p_800.csv",  "figs/p_800.png",  "p (t=800)"),
    ("output/old_runs/p_1200.csv", "figs/p_1200.png", "p (t=1200)"),
    ("output/old_runs/p_1600.csv", "figs/p_1600.png", "p (t=1600)"),
    ("output/old_runs/p_2000.csv", "figs/p_2000.png", "p (t=2000)"),

    # MPI final fields
    ("output/scaling/mpi_np1_rho.csv", "figs/mpi_np1_rho_final.png",
     "MPI np=1 rho (final)"),
    ("output/scaling/mpi_np1_p.csv",   "figs/mpi_np1_p_final.png",
     "MPI np=1 p (final)"),

    ("output/scaling/mpi_np2_rho.csv", "figs/mpi_np2_rho_final.png",
     "MPI np=2 rho (final)"),
    ("output/scaling/mpi_np2_p.csv",   "figs/mpi_np2_p_final.png",
     "MPI np=2 p (final)"),

    ("output/scaling/mpi_np4_rho.csv", "figs/mpi_np4_rho_final.png",
     "MPI np=4 rho (final)"),
    ("output/scaling/mpi_np4_p.csv",   "figs/mpi_np4_p_final.png",
     "MPI np=4 p (final)"),
]

for c, p, t in todo:
    if os.path.exists(c):
        plot_field(c, p, t)
        print("Wrote", p)
    else:
        print("missing:", c)

print("done -> figs/")
