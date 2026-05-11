import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs("figs", exist_ok=True)

a = np.loadtxt("output/omp/p_final_ompN1.csv", delimiter=",")
b = np.loadtxt("output/omp/p_final_ompN4.csv", delimiter=",")

d = np.abs(a - b)

plt.figure(figsize=(7,4))
plt.imshow(d, origin="lower", aspect="auto")
plt.colorbar()
plt.title("OpenMP |p(N=1) - p(N=4)| at final time")
plt.tight_layout()
plt.savefig("figs/omp_np1_np4_p_absdiff.png", dpi=300)
print("Wrote figs/omp_np1_np4_p_absdiff.png")
