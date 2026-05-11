#!/usr/bin/env python3
import os, glob
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Read rho_xxx.csv/p_xxx.csv/rho_diff_500.csv/p_diff_500.csv

# Paths
BASE   = Path(__file__).resolve().parents[1]
DATADIR = os.path.join(BASE, "data")
OTHERS = os.path.join(DATADIR, "others")
OUTDIR = os.path.join(BASE, "project_figs")

CPU1 = os.path.join(DATADIR, "data_1st")
CPU2 = os.path.join(DATADIR, "data_2nd")
GPU1 = os.path.join(DATADIR, "data_gpu_1st")
GPU2 = os.path.join(DATADIR, "data_gpu_2nd")
DIFF2 = OTHERS   # p_diff_500.csv / rho_diff_500.csv now live in data/others

STEP0 = 0
STEP1_LAST = 530
STEP2_LAST = 540
STEP_LOGDIFF = 500

DPI = 450
LOG_EPS = 1e-30
MASK_TOL = 0.0

def ensure_dir(d):
    os.makedirs(d, exist_ok=True)

def load_csv(folder, name, step):
    path = os.path.join(folder, f"{name}_{step}.csv")
    if not os.path.exists(path):
        return None
    return np.loadtxt(path, delimiter=",")

def load_diff_csv(folder, name, step):
    path = os.path.join(folder, f"{name}_diff_{step}.csv")
    if not os.path.exists(path):
        return None
    return np.loadtxt(path, delimiter=",")

def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close()

def heatmap(Z, title, cbar_label, xlabel="x index", ylabel="y index"):
    plt.figure(figsize=(10, 4))
    plt.imshow(Z, origin="lower", aspect="auto")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    cb = plt.colorbar()
    cb.set_label(cbar_label)

def diff_mask(A, B, tol, title):
    D = np.abs(A - B)
    M = (D > tol).astype(np.int32)
    plt.figure(figsize=(10, 4))
    plt.imshow(M, origin="lower", aspect="auto", vmin=0, vmax=1)
    plt.title(title + f" (tol={tol:g})")
    plt.xlabel("x index")
    plt.ylabel("y index")
    cb = plt.colorbar()
    cb.set_label("0=same, 1=different")
    return M

def triple_panel(A, B, D, titleA, titleB, titleD, cbarA, cbarB, cbarD, outpath):
    plt.figure(figsize=(10, 10))

    plt.subplot(3, 1, 1)
    plt.imshow(A, origin="lower", aspect="auto")
    plt.title(titleA)
    plt.colorbar(label=cbarA)

    plt.subplot(3, 1, 2)
    plt.imshow(B, origin="lower", aspect="auto")
    plt.title(titleB)
    plt.colorbar(label=cbarB)

    plt.subplot(3, 1, 3)
    plt.imshow(D, origin="lower", aspect="auto")
    plt.title(titleD)
    plt.colorbar(label=cbarD)

    savefig(outpath)

def midline_plot(arrs, labels, title, ylabel, outpath):
    plt.figure(figsize=(12, 3))
    for A, lab in zip(arrs, labels):
        mid = A[A.shape[0] // 2, :]
        plt.plot(mid, label=lab, linewidth=2)
    plt.title(title)
    plt.xlabel("x (cell index)")
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True, alpha=0.25)
    savefig(outpath)

def list_steps(folder, name):
    patt = os.path.join(folder, f"{name}_*.csv")
    steps = []
    for p in glob.glob(patt):
        b = os.path.basename(p)
        try:
            s = int(b.replace(f"{name}_", "").replace(".csv", ""))
            steps.append(s)
        except:
            pass
    return sorted(steps)

def max_abs_error_curve(folderA, folderB, name, outpath, title):
    stepsA = set(list_steps(folderA, name))
    stepsB = set(list_steps(folderB, name))
    steps = sorted(list(stepsA & stepsB))
    if not steps:
        print(f"[skip] no common steps for {name} between {folderA} and {folderB}")
        return

    errs = []
    for s in steps:
        A = load_csv(folderA, name, s)
        B = load_csv(folderB, name, s)
        errs.append(float(np.max(np.abs(A - B))))

    plt.figure(figsize=(10, 4))
    plt.plot(steps, errs, linewidth=2)
    plt.title(title)
    plt.xlabel("step")
    plt.ylabel(f"max |Δ{name}|")
    plt.grid(True, alpha=0.25)
    savefig(outpath)

def main():
    ensure_dir(OUTDIR)
    print("Writing figures to:", OUTDIR)

    # 1) step 0 snapshots
    for fld, cbar in [("rho", "density ρ"), ("p", "pressure p (Pa)")]:
        for tag, folder in [("cpu1", CPU1), ("cpu2", CPU2), ("gpu1", GPU1), ("gpu2", GPU2)]:
            A = load_csv(folder, fld, STEP0)
            if A is None:
                print(f"[skip] {tag} {fld} step0 missing")
                continue
            heatmap(A, f"{tag.upper()} {fld} (step {STEP0})", cbar)
            savefig(os.path.join(OUTDIR, f"{fld}_{tag}_step{STEP0}.png"))

    # optional HD aliases matching list
    for fld, folder, outname in [
        ("rho", CPU1, "rho_0_hd.png"),
        ("p", CPU1, "p_0_hd.png"),
    ]:
        A = load_csv(folder, fld, STEP0)
        if A is not None:
            heatmap(A, f"{fld} step 0", fld)
            savefig(os.path.join(OUTDIR, outname))

    # 2) final snapshots
    for fld, cbar in [("rho", "density ρ"), ("p", "pressure p (Pa)")]:
        for tag, folder, step in [("cpu1", CPU1, STEP1_LAST), ("gpu1", GPU1, STEP1_LAST)]:
            A = load_csv(folder, fld, step)
            if A is None:
                print(f"[skip] {tag} {fld} step{step} missing")
                continue
            heatmap(A, f"{tag.upper()} {fld} (step {step})", cbar)
            savefig(os.path.join(OUTDIR, f"{fld}_{tag}_step{step}.png"))

        for tag, folder, step in [("cpu2", CPU2, STEP2_LAST), ("gpu2", GPU2, STEP2_LAST)]:
            A = load_csv(folder, fld, step)
            if A is None:
                print(f"[skip] {tag} {fld} step{step} missing")
                continue
            heatmap(A, f"{tag.upper()} {fld} (step {step})", cbar)
            savefig(os.path.join(OUTDIR, f"{fld}_{tag}_step{step}.png"))

    # 3) midline profiles at step 500
    for fld, ylab in [("rho", "rho"), ("p", "p")]:
        arrs, labels = [], []
        for lab, folder in [("CPU 1st", CPU1), ("CPU 2nd", CPU2), ("GPU 1st", GPU1), ("GPU 2nd", GPU2)]:
            A = load_csv(folder, fld, STEP_LOGDIFF)
            if A is not None:
                arrs.append(A)
                labels.append(lab)
        if len(arrs) >= 2:
            midline_plot(
                arrs, labels,
                f"Midline profiles (y=mid), {fld}, step={STEP_LOGDIFF}",
                ylab,
                os.path.join(OUTDIR, f"{fld}_midline_step{STEP_LOGDIFF}_all.png")
            )

    # 4) CPU1 vs GPU1 at step 500
    for fld, cbar in [("rho", "ρ"), ("p", "p (Pa)")]:
        A = load_csv(CPU1, fld, STEP_LOGDIFF)
        B = load_csv(GPU1, fld, STEP_LOGDIFF)
        if A is None or B is None:
            print(f"[skip] CPU1 vs GPU1 {fld} step500 missing")
            continue

        D = np.abs(A - B)
        heatmap(D, f"|{fld}_CPU1 - {fld}_GPU1| (step {STEP_LOGDIFF})", f"|Δ{fld}|")
        savefig(os.path.join(OUTDIR, f"{fld}_absdiff_cpu1_gpu1_step{STEP_LOGDIFF}.png"))

        L = np.log10(D + LOG_EPS)
        heatmap(L, f"log10(|{fld}_CPU1 - {fld}_GPU1|) (step {STEP_LOGDIFF})", f"log10(|Δ{fld}|)")
        savefig(os.path.join(OUTDIR, f"{fld}_log10_absdiff_cpu1_gpu1_step{STEP_LOGDIFF}.png"))

        diff_mask(A, B, MASK_TOL, f"{fld} diff mask CPU1 vs GPU1 step={STEP_LOGDIFF}")
        savefig(os.path.join(OUTDIR, f"{fld}_diffmask_cpu1_gpu1_step{STEP_LOGDIFF}.png"))

        triple_panel(
            A, B, D,
            f"CPU 1st {fld} (step {STEP_LOGDIFF})",
            f"GPU 1st {fld} (step {STEP_LOGDIFF})",
            f"abs diff |CPU1-GPU1| (step {STEP_LOGDIFF})",
            cbar, cbar, f"|Δ{fld}|",
            os.path.join(OUTDIR, f"{fld}_triple_cpu1_gpu1_absdiff_step{STEP_LOGDIFF}.png")
        )

    # 5) CPU2 vs GPU2 at step 500
    for fld, cbar in [("rho", "ρ"), ("p", "p (Pa)")]:
        A = load_csv(CPU2, fld, STEP_LOGDIFF)
        B = load_csv(GPU2, fld, STEP_LOGDIFF)
        Dfile = load_diff_csv(DIFF2, fld, STEP_LOGDIFF)
        if A is None or B is None:
            print(f"[skip] CPU2/GPU2 {fld} step500 missing")
            continue

        D = np.abs(A - B) if Dfile is None else np.abs(Dfile)

        heatmap(D, f"|{fld}_CPU2 - {fld}_GPU2| (step {STEP_LOGDIFF})", f"|Δ{fld}|")
        savefig(os.path.join(OUTDIR, f"{fld}_absdiff_cpu2_gpu2_step{STEP_LOGDIFF}.png"))

        L = np.log10(D + LOG_EPS)
        heatmap(L, f"log10(|{fld}_CPU2 - {fld}_GPU2|) (step {STEP_LOGDIFF})", f"log10(|Δ{fld}|)")
        savefig(os.path.join(OUTDIR, f"{fld}_log10_absdiff_cpu2_gpu2_step{STEP_LOGDIFF}.png"))

        diff_mask(A, B, MASK_TOL, f"{fld} diff mask CPU2 vs GPU2 step={STEP_LOGDIFF}")
        savefig(os.path.join(OUTDIR, f"{fld}_diffmask_cpu2_gpu2_step{STEP_LOGDIFF}.png"))

        triple_panel(
            A, B, D,
            f"CPU 2nd {fld} (step {STEP_LOGDIFF})",
            f"GPU 2nd {fld} (step {STEP_LOGDIFF})",
            f"abs diff |CPU2-GPU2| (step {STEP_LOGDIFF})",
            cbar, cbar, f"|Δ{fld}|",
            os.path.join(OUTDIR, f"{fld}_triple_cpu2_gpu2_absdiff_step{STEP_LOGDIFF}.png")
        )

    # 6) Error-vs-step curves
    max_abs_error_curve(
        CPU1, GPU1, "rho",
        os.path.join(OUTDIR, "rho_max_abs_error_vs_step_cpu1_gpu1.png"),
        "Max abs error vs step (rho, CPU 1st vs GPU 1st)"
    )
    max_abs_error_curve(
        CPU1, GPU1, "p",
        os.path.join(OUTDIR, "p_max_abs_error_vs_step_cpu1_gpu1.png"),
        "Max abs error vs step (p, CPU 1st vs GPU 1st)"
    )
    max_abs_error_curve(
        CPU2, GPU2, "rho",
        os.path.join(OUTDIR, "rho_max_abs_error_vs_step_cpu2_gpu2.png"),
        "Max abs error vs step (rho, CPU 2nd vs GPU 2nd)"
    )
    max_abs_error_curve(
        CPU2, GPU2, "p",
        os.path.join(OUTDIR, "p_max_abs_error_vs_step_cpu2_gpu2.png"),
        "Max abs error vs step (p, CPU 2nd vs GPU 2nd)"
    )

    print("Done.")

if __name__ == "__main__":
    main()
