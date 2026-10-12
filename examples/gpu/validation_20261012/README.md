# Serial and GPU boundary validation, October 12, 2026

Particles leaving the radial domain now terminate independently with finite
saved states. Serial fixed-step and adaptive tracing roll back the rejected
step and its output. The device solver uses `isoutside(u,p,t)`. Interior NaNs
remain numerical failures. Absorption stays fixed at 200 km altitude, even
when a field mesh extends down to 100 km. Initially absorbed particles stop
at the initial time without querying fields or aborting their batch.

`trace_forward` saves the last accepted integration node with
`ReturnCode.Terminated`. `trace_forward_bounded` saves exact sphere intersections
with `:inner`, `:outer`, or `:time_limit`; initially absorbed particles retain
their initial state and have `:inner` at the initial time.

The environment uses Julia 1.12.6, four CPU threads, TestParticle 0.24.2, and
an NVIDIA GeForce RTX 4070 SUPER. Registered BorisPushers 0.1.1 does not contain
the new device boundary checks, so Project.toml and Manifest.toml select its
upstream subdirectory at commit `7da93a7612d3c36a6f39377cb87ee15bd0d4b5a5`.
Its declared version is still 0.1.1. MHD readers now use complete stored axis
metadata; legacy files retain their historical grid convention.

## Actual MHD checks

Two read-only fields were tested: the legacy 120 x 72 x 144 grid starting at
200 km and the pe_hall 400 x 61 x 121 grid starting at 100 km. Both extend to
4 Mars radii. Rm = 3390 km. Vector components follow the simulation Cartesian
axes; no new MSO/MSE transformation was applied. E and B are static, in V/m
and T. The pe_hall source header identifies its axes as GSM.

Each field used 20 deterministic O2+ initial states along the six Cartesian
axis directions: six inward, six outward, six interior, one inward from the
absorption boundary and one initially at 99 km. Boundary velocities were
100 km/s. The fixed timestep was 0.005 s, with a maximum time of 1 s. Ordinary
serial and GPU tracing were compared with every-step and 0.02 s saving.
Exact sphere stopping was compared between sequential single-particle CPU
calls and a GPU batch. All saved states were finite and statuses agreed.

| Comparison | Maximum position difference (m) | Maximum velocity difference (m/s) | Maximum end-time difference (s) |
| --- | ---: | ---: | ---: |
| Ordinary serial versus GPU | 1.463e-8 | 1.463e-6 | 0 |
| Exact spherical CPU versus GPU | 1.464e-8 | 1.463e-6 | 1.910e-14 |

Each exact-sphere batch returned eight `:inner`, six `:outer` and six
`:time_limit` states. Moving inward particles ended at 200 km to within
4.55e-13 km. Additional initial heights of 99, 50, 0 and -1 km returned finite
unchanged initial states at t=0 on both serial and GPU paths.

The actual-field script passed 1047 checks. Full regression, including CUDA
spherical tests, passed 868 checks. The deliberately invalid-field cases test
failure handling; they do not indicate NaNs in these actual field inputs.

The CSV files beside this document record each particle's comparisons.
Large field inputs and solver snapshots are not included in the repository.

## Reproduce

From the project directory, with the appropriate VTS inputs available:

```powershell
julia --threads=4 --project=. test/runtests.jl
$env:MARSTP_TEST_CUDA = 'true'
julia --threads=4 --project=. test/runtests.jl
julia --threads=4 --project=. examples/gpu/actual_mhd_domain_validation.jl legacy.vts pe_hall_100km.vts
```

The actual-field script writes to a new output directory and preserves inputs.
