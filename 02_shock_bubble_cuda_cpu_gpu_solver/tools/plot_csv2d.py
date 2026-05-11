#!/usr/bin/env python3
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

# Read any two-dimensional matrix CSV

def main():
    ap = argparse.ArgumentParser(description="Plot 2D CSV (matrix) as heatmap (hi-res).")
    ap.add_argument("input", help="path to input .csv (2D matrix)")
    ap.add_argument("-o", "--output", default=None, help="output .png path")
    ap.add_argument("--title", default=None, help="figure title")
    ap.add_argument("--cblabel", default=None, help="colorbar label")
    ap.add_argument("--cmap", default="viridis", help="matplotlib colormap (default viridis)")
    ap.add_argument("--dpi", type=int, default=300, help="output dpi (default 300)")
    ap.add_argument("--vmin", type=float, default=None, help="colorbar min")
    ap.add_argument("--vmax", type=float, default=None, help="colorbar max")
    ap.add_argument("--clip-percent", type=float, default=None,
                    help="optional robust clipping by percent (e.g. 1 means clip to [1,99])")
    ap.add_argument("--no-axes", action="store_true", help="hide axis ticks/labels")
    args = ap.parse_args()

    data = np.loadtxt(args.input, delimiter=",")

    # robust clipping to avoid one shock cell ruining the color scale
    vmin, vmax = args.vmin, args.vmax
    if args.clip_percent is not None:
        p = float(args.clip_percent)
        lo, hi = np.percentile(data, [p, 100.0 - p])
        if vmin is None: vmin = lo
        if vmax is None: vmax = hi

    # output name
    if args.output is None:
        base = os.path.splitext(os.path.basename(args.input))[0]
        args.output = base + ".png"

    fig_w = 10.5
    fig_h = 4.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(
        data,
        origin="lower",
        aspect="auto",
        cmap=args.cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )

    if args.title:
        ax.set_title(args.title)

    if args.no_axes:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.set_xlabel("x index")
        ax.set_ylabel("y index")

    cb = fig.colorbar(im, ax=ax)
    if args.cblabel:
        cb.set_label(args.cblabel)

    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)
    plt.close(fig)
    print(f"[OK] wrote {args.output}")

if __name__ == "__main__":
    main()

