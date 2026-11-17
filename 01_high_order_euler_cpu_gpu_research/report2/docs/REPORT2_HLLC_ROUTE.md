# Report 2 HLLC Route

The Report 2 solver in `code/hllc_solver/` uses HLLC interface fluxes with
WENO-2, WENO-3, and ADER-3. This note describes the flux and CFL selection.

## Active Methods

- WENO-2 uses HLLC for the Euler interface flux.
- WENO-3 uses HLLC for the Euler interface flux.
- ADER-3 uses HLLC for the time-averaged Euler interface flux.

Godunov and MUSCL-Hancock remain as historical Report 1 methods. Their CPU
implementations are not the active Report 2 development target. HLL/Rusanov may
remain as a degenerate-denominator safety fallback inside HLLC functions, but
HLLC is the official Report 2 flux choice.

## HLLC Implementation Notes

The HLLC solver uses left and right acoustic wave-speed estimates,

```text
S_L = min(u_L - c_L, u_R - c_R),
S_R = max(u_L + c_L, u_R + c_R),
```

with the corresponding normal velocity in the y-direction. The contact speed is
computed from the Rankine-Hugoniot relation,

```text
S_M = (p_R - p_L + rho_L u_L (S_L - u_L)
       - rho_R u_R (S_R - u_R))
      / (rho_L (S_L - u_L) - rho_R (S_R - u_R)).
```

The one-dimensional star states use `(rho*, rho* S_M, E*)`. In two dimensions,
the normal momentum is replaced by `rho* S_M` and the tangential velocity is
preserved, so the x-flux star state uses `rho* v_K` while the y-flux star state
uses `rho* u_K`.

Flux selection follows the standard HLLC wave regions:

- return the left physical flux if `S_L >= 0`;
- return the right physical flux if `S_R <= 0`;
- otherwise return the left or right star flux according to the sign of `S_M`.

## Current CFL Decision

The current official discontinuous CFL values are:

```text
WENO-2: 0.20
WENO-3: 0.20
ADER-3: 0.14  (official selected value)
```

ADER-3 uses a lower discontinuous CFL because the positivity/retry logic is
sensitive to small floating-point and reduction-order differences near strong
interfaces. The important observation is that the CFL effect is not monotone.
For a discontinuous problem, a smaller CFL does not simply mean better backend
agreement; it changes the sequence of intermediate states and can place HLLC
wave-speed selection or fallback decisions closer to a branch threshold.

The current diagnostic record is:

- `quadrant_riemann_2d`: the existing `100 x 100` probe recorded the following
  serial/CUDA maximum field differences:

  | CFL | max serial/CUDA difference | Note |
  |---:|---:|---|
  | 0.15 | `4.93e-12` | roundoff-level agreement |
  | 0.12 | `3.36e-3` | different fallback counts; visible difference |
  | 0.10 | `1.60e-12` | roundoff-level agreement |
  | 0.08 | `1.60e-12` | roundoff-level agreement |
  | 0.06 | `2.78e-5` | small but visible above roundoff |

- `shock_bubble_2d`: a single-case `CFL=0.12` confirmation gave roundoff-level
  serial/CUDA differences, with maximum component differences below `6e-8`.
  The earlier shock-bubble CFL probe also recorded `CFL=0.15`, but that setting
  was not attractive for backend consistency: serial and CUDA took different
  step counts (`1124` vs `1116`), different rejected-step counts (`11` vs `13`),
  and different fallback counts. That probe did not store exported field
  differences for `CFL=0.15`, so it is treated as diagnostic timing/path
  evidence rather than a final field-difference result. The refreshed full
  backend-consistency check at `CFL=0.10` still showed an ADER-3 CUDA difference
  of about `4e-4`, even though WENO-2 and WENO-3 remained at roundoff level.
- The two-case intermediate probe then tested `CFL=0.13` and `0.14`.
  `CFL=0.13` kept quadrant Riemann at roundoff level but gave a large
  shock-bubble serial/CUDA difference of about `1.58e-2`. `CFL=0.14` kept
  quadrant Riemann at roundoff level and reduced the shock-bubble difference to
  about `2.02e-7`, with matching step/retry/fallback counts. The value `0.14`
  is therefore the official ADER-3 discontinuous CFL for Report 2.

## Lovelace Entry Points

Backend consistency:

```powershell
Set-Location <package-root>\workflows
.\run_report2_hllc_backend_consistency_from_windows.ps1 -User <user> -HostName <host>
```

Phase-profile CUDA smoke and resume:

```powershell
Set-Location <package-root>\workflows
.\smoke_report2_hllc_phase_profile_from_windows.ps1 -User <user> -HostName <host>
```

These wrappers preserve the original collection environment and therefore
still contain explicit workspace and remote-account roots internally. Update
those variables before rerunning from another checkout or account.
