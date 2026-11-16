Stage C

Contents:
- Runtime benchmark CSV files
- Runtime scaling studies
- Accuracy-vs-runtime merged tables
- Efficiency / speedup figures

Note:
- These files are the main performance-study archive.
- The Lovelace A30 no-I/O files are the final timing data used by the report
  for the all-backend 2D Gaussian runtime-scaling, speedup, and
  accuracy-runtime checks. These runs exclude CSV field-output time by using
  --output=none so that the timing focuses on numerical update and backend
  execution cost.
- The discontinuous Chapter 6 performance table uses the
  runtime_discontinuous_2d_lovelace_a30_n600_noio.csv file, also measured with
  --output=none.
- The Lovelace OpenMP thread-scaling and discontinuous runtime checks are
  retained as supporting performance data.
- The all-backend OpenMP timing runs use OMP_NUM_THREADS=16 and
  OMP_DYNAMIC=FALSE; the separate OpenMP thread-scaling study intentionally
  varies OMP_NUM_THREADS.

