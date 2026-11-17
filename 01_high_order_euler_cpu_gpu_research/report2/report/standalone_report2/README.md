# Formal Report 2 source

This folder is a self-contained LaTeX project and can be uploaded directly to Overleaf.

## Build

Set `main.tex` as the main document and compile with pdfLaTeX. Compile twice so the table of contents and cross-references are resolved.

## Report scope

The report covers the Report 1 bridge, integrated numerical methodology,
convergence validation, backend/physical verification with the discontinuous-flow
evolution and slice analysis, the main accuracy--performance analysis, integrated
discussion, and conclusions.
The chapter order follows the Report 2 marking priorities without displaying
artificial `Part I` or `Part II` headings in the report itself:

- **Project development material (20%)**: a short Report 1 bridge, one integrated Methodology chapter containing the new HLLC/WENO-3/ADER-3 theory and compact implementation listings, followed by two complete 1D convergence studies and a finest-interval 2D confirmation;
- **Computational evidence**: backend consistency, CFL sensitivity, the exact Sod reference,
  accepted high-resolution fields, discontinuous evolution and density slices, smooth and
  discontinuous performance,
  matched HLL--HLLC comparisons, phase-resolved bottleneck analysis, OpenMP thread scaling,
  and the official CUDA-only scaling sequence through `5000 x 5000`.

## Structure

- `main.tex`: title page, packages, and section order.
- `sections/`: editable report chapters.
- `figures/`: only the confirmed figures used by the report.

## Results

Numerical values and figures use the HLLC experiments in
`../../results/hllc/`. The accuracy-performance chapter includes the
matched HLL/Rusanov comparison from `../../analysis/hll_vs_hllc/`.

The [results summary](../../docs/results.md) describes the experiments and
the WENO-2 OpenMP correction used in the performance tables.
