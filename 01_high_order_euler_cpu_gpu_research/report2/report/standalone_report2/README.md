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

## Data provenance

The numerical values and figures in this report come from the canonical HLLC
results stored in this package under `../../results/hllc/`.

No HLL/Rusanov performance data are used in the validation or backend-consistency
chapters. HLL appears only as a strictly matched comparator in the accuracy--performance
chapter; HLLC remains authoritative.

The full formal-report audit dated 20 July 2026 is recorded in
`../../docs/FORMAL_REPORT_HLLC_AUDIT_2026-07-20.md`; the earlier data-version
audit is retained in `../../docs/DATA_VERSION_AUDIT_2026-07-19.md`. The July
WENO-2 OpenMP correction changes performance measurements, not the numerical
fields: a fresh local serial/OpenMP replay of all three Stage B WENO-2 cases
produced byte-identical CSV files. Stage C performance values are taken from
the post-refresh canonical tables included in this package.
