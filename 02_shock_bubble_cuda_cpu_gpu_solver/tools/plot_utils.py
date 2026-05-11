#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def load_csv(path):
    return np.loadtxt(path, delimiter=",")

def save_imshow(data, title, fname, cmap="viridis", log=False, dpi=250):
    plt.figure(figsize=(9, 3.2))
    if log:
        data = np.log10(data + 1e-16)
        title += " (log10)"
    im = plt.imshow(data, origin="lower", aspect="auto", cmap=cmap)
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.title(title)
    plt.xlabel("i")
    plt.ylabel("j")
    plt.tight_layout()
    plt.savefig(fname, dpi=dpi)
    plt.close()