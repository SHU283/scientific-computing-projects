import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs("figs", exist_ok=True)

p1_path = "output/omp/p_final_ompN1.csv"
p4_path = "output/omp/p_final_ompN4.csv"

a = np.loadtxt(p1_path, delimiter=",")
b = np.loadtxt(p4_path, delimiter=",")

vmin = min(a.min(), b.min())
vmax = max(a.max(), b.max())

fig = plt.figure(figsize=(12, 3.9), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.05])

ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1], sharex=ax0, sharey=ax0)
cax = fig.add_subplot(gs[0, 2])

im0 = ax0.imshow(a, origin="lower", aspect="auto",
                 cmap="viridis", interpolation="bilinear",
                 vmin=vmin, vmax=vmax)
ax0.set_title(r"OpenMP $p$ (N=1)")
ax0.set_xlabel("i")
ax0.set_ylabel("j")

im1 = ax1.imshow(b, origin="lower", aspect="auto",
                 cmap="viridis", interpolation="bilinear",
                 vmin=vmin, vmax=vmax)
ax1.set_title(r"OpenMP $p$ (N=4)")
ax1.set_xlabel("i")
ax1.set_ylabel("j")

cb = fig.colorbar(im1, cax=cax)
cb.set_label(r"Pressure $p$ (Pa)")

out = "figs/omp_np1_np4_p_final_side_by_side.png"
plt.savefig(out, dpi=300)
plt.close()
print("Wrote", out)

print("L_inf |Δp| =", np.max(np.abs(a - b)))
