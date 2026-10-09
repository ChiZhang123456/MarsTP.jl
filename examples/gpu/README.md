# CPU/GPU tracing comparison

The project uses TestParticle 0.24.1 and its native BorisPushers ensemble kernel,
CUDA and KernelAbstractions. CUDA is loaded by the example, not on CPU-only
MarsTP imports. Julia 1.12.6 and an NVIDIA RTX 4070 SUPER were used locally.

Run from the repository root:

```powershell
julia --threads=4 --project=. examples/gpu/benchmark_cpu_gpu.jl
```

The script checks zero-field free motion, proton gyration direction, magnetic
energy conservation, timestep convergence and CPU/GPU agreement. It compares
100, 1000 and 10000 particles, each with 1000 fixed Boris steps, in analytic
uniform fields and the existing MarsTP MHD fields. All calculations use Float64.
Each result is the median of three warmed complete solves. Allocation, GPU
upload/download and host solution construction are included. First compilation
and VTK file loading are excluded. Raw trials, input SHA256, software versions,
hardware, seed and numerical diagnostics are saved in a new timestamped
`outputs/gpu_benchmark_*` directory. No full trajectory files are saved.

The MHD example traces O2+ for 1 s at dt=0.001 s, saves every 0.1 s, and starts
at Cartesian (1.5,0.1,0.4) Mars radii with seeded Gaussian velocities of component
standard deviation 10 km/s. Radius is 3390 km. Fields are static, with the total
electric field, SI units and the existing spherical interpolation/basis
conversion. The axes match the source MHD data; the benchmark does not claim a
new MSO/MSE transformation. No collisions, gravity or feedback are included.

## Use GPU through MarsTP

```julia
using MarsTP, CUDA, StaticArrays
fields = load_mhd_fields() # cache input, reused across batches
states = [SA[1.5Rm, 0.1Rm, 0.4Rm, 1e4, 0., 0.]]
cfg = ForwardTraceConfig(species="O2+", solver=:boris, dt=.001, tspan=(0.,1.))
solutions = trace_forward(states; config=cfg, fields,
    backend=CUDA.CUDABackend(), saveat=.1)
trajectory = only(solutions[1].u)
```

Omit `backend` for the original serial CPU path; use `KernelAbstractions.CPU()`
for the threaded ensemble kernel. CPU/GPU use identical timesteps, precision
and output cadence. MarsTP preserves the vector of single-member ensembles
returned by `trace_forward`. TestParticle itself changed its single-particle
`solve(prob, Boris(); ...)` result in 0.24 to an `ODESolution`, so standalone
older examples that directly use TestParticle must read `sol.t` and `sol.u`
instead of `only(sol.u).t` and `only(sol.u).u`.

GPU covers fixed-step forward tracing. Adaptive tracing and the
source-accumulating backward VDF remain on CPU. `trace_forward` uses the
upstream time-limited kernel; `trace_forward_bounded` adds a MarsTP kernel that
stops each particle independently at spherical absorption/escape boundaries.
Nonfinite states remain failures, with no field extrapolation or replacement.
The short benchmark below measures `trace_forward`; it does not measure the
new boundary kernel or the complete Monte Carlo or backward VDF pipeline.

## Per-particle absorption and escape

```julia
using MarsTP, CUDA
cfg = ForwardTraceConfig(species="O2+",solver=:boris,dt=.1,tspan=(0.,500.))
trajectories = trace_forward_bounded(initial_states;config=cfg,
    backend=CUDA.CUDABackend(),inner_radius_m=Rm+200e3,
    outer_radius_m=Router,save_every=10)
statuses = getproperty.(trajectories,:status)
```

Each returned trajectory has `t`, synchronized Cartesian `u` (m, m/s),
`species`, `retcode`, and `status`: `:inner`, `:outer`, `:time_limit`, or
`:numerical_failure`. Absorption is at 200 km altitude by default; escape is
at 4 Mars radii. Source sampling at 500 km and n*norm(U) source weights are
unchanged and are supplied separately. Field coverage must include both spheres.

Every Boris drift is checked for its first sphere intersection, even when
both endpoints lie outside the inner sphere. The final state lies at that
intersection and the time includes its fractional step. Endpoint velocity
uses the same clipped-step synchronization as the CPU shell tracer. Check
step-size convergence for curved trajectories and termination diagnostics.
Boundary checks precede field evaluation at a trial endpoint, preventing
out-of-domain interpolation. A one-micrometre query guard places boundary
field evaluation on the covered side; the saved endpoint is unshifted.

`save_every=0` saves only initial and final states; a positive integer saves
every N steps plus the terminal state. Run bounded particle batches to control
GPU memory. Detector residence and postprocessed work diagnostics need sufficiently
dense histories. Pass `work_itp` to accumulate all integration-step work on the
device, independently of output thinning; the returned `work` includes signed
total, convection and Hall work in eV, their positive/negative contributions,
kinetic-energy change and closure residuals. General callbacks and on-device
detector accumulators are not implemented.

## Production Monte Carlo GPU path

The shell driver now accepts `tracing_backend=:cuda`. It uses the GPU boundary
kernel and device work accumulation, retains n*norm(U) injection at 500 km and
absorption at 200 km, and processes detector residence/crossings on the CPU
from every accepted integration step. All physical source weights are unchanged.

```julia
include("examples/forward_tracing/monte_carlo_forward_tracing/monte_carlo_shell.jl")
ShellMonteCarlo.run_monte_carlo("outputs/my_gpu_run",
    ShellMonteCarlo.Config(tracing_backend=:cuda,work_mode=:summary,
        batch_size=256,dt=.1,tmax=500.,per_cell=100))
```

Use `work_mode=:summary` with `save_power=false` to write the device-computed
work summaries directly to the existing JLD2 format. `:steps` or endpoint
power additionally evaluates diagnostics on the host to retain the existing
per-step/power arrays. `particles.csv` preserves `work_eV` (total), and adds
`work_conv_eV`, `work_hall_eV`, and `field_sum_residual_eV`. Existing kinetic
closure and termination columns are retained. Metadata records the backend;
the numerical kernel is copied into each run's source snapshots.

The GPU accelerates integration and work accumulation. Source sampling,
detector analysis, field uploads, trajectory transfers, compression and disk
output still contribute to end-to-end runtime; the earlier short tracing
speed ratios do not describe this full production path. Full histories are
kept for detector estimates, so choose batches according to steps and GPU RAM
(roughly 56 bytes per particle per saved state, before fields and diagnostics).

Small end-to-end actual-MHD validation (60 particles per backend):
`julia --threads=4 --project=. examples/gpu/production_smoke.jl`.

Small actual-field CPU/GPU comparison:
`julia --threads=4 --project=. examples/gpu/bounded_forward.jl`.
To include device checks in the package tests, set `MARSTP_TEST_CUDA=true`.

Initial boundary validation on 2026-10-08 (Julia 1.12.6, TestParticle 0.24.1, RTX 4070 SUPER):
550 package assertions passed, including 47 CPU boundary assertions; the same
47 device assertions passed on CUDA. Cases cover independent stopping times,
initial boundary positions, a drift through the entire inner sphere, invalid
fields, synchronized magnetic motion, and magnetic-event timestep convergence.
The actual-MHD example passed 11 CPU/GPU comparison assertions: absorption at
200 km at 0.0099995681 s, escape at 4 Rm at 0.0099999996 s, and a third particle
reaching 10 s. These are small validation cases, not a full production ensemble.

After production/work integration, 674 package assertions passed and the
152 boundary/work assertions passed on CUDA. The shell suite passed 435
assertions with both CPU-kernel and CUDA adapter comparisons against the
original CPU tracer. Four Python analysis tests passed. The 60-particle,
0.2 s actual-MHD end-to-end run passed 1274 output assertions (plus two
record-count checks), including source weights, CSV work columns and JLD2
summaries. Maximum CPU/GPU discrepancies were 1.60e-14 eV for the three work
totals, 2.33e-10 m for position, and 4.07e-12 m/s for velocity. Compilation
is included in the smoke-run elapsed times; they are not performance measurements.

The same end-to-end comparison extended to 30 s also passed all 1274 output
assertions: both backends absorbed 14 particles at 200 km and kept 46 until the
time limit. Maximum discrepancies were 2.73e-12 eV for work, 1.92e-9 m for
position and 5.03e-10 m/s for velocity. Escape and partial terminal-step work
are additionally covered by the analytic production-adapter tests. Run this
longer small check with `examples/gpu/production_smoke.jl 30`; it remains a
validation sample, not a converged production population.

The pre-update Project.toml and Manifest.toml were preserved in
`outputs/gpu_update_20261008_1349/`.

## TestParticle 0.24.1 measurements (2026-10-08)

The new run uses the same particle counts, initial-condition seed, Float64,
1000 steps per particle and saved output cadence as the 0.24.0 example.
Times below are medians of three warmed complete solves, including field
adaptation, transfers and solution construction, excluding compilation and
VTK loading. No other Julia processes were running during the measurements.

Archived [raw timing trials](benchmarks/testparticle_0.24.1/timings.csv) and
[metadata and physics checks](benchmarks/testparticle_0.24.1/results.toml)
are small benchmark records; field arrays and trajectories are not bundled.
The metadata's MarsTP commit identifies the baseline checkout with local
modifications used for the run, before the publication commit.

MarsTP MHD forward tracing:

| Particles | CPU serial (s) | CPU 4 threads (s) | GPU (s) | Serial/GPU | Threads/GPU |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.042861 | 0.026749 | 0.061594 | 0.70 | 0.43 |
| 1000 | 0.267259 | 0.065351 | 0.077209 | 3.46 | 0.85 |
| 10000 | 2.335658 | 0.371996 | 0.093749 | 24.91 | 3.97 |

Analytic uniform-field tracing:

| Particles | CPU serial (s) | CPU 4 threads (s) | GPU (s) |
| ---: | ---: | ---: | ---: |
| 100 | 0.018099 | 0.004951 | 0.005563 |
| 1000 | 0.116426 | 0.042634 | 0.026092 |
| 10000 | 1.230331 | 0.374732 | 0.032753 |

Hardware and supporting dependencies are unchanged: i7-14700KF, four Julia
threads, RTX 4070 SUPER, CUDA.jl 6.4.1, KernelAbstractions 0.9.44,
BorisPushers 0.1.1 and Julia 1.12.6. At 10000 MHD particles the GPU is 24.9
times faster than serial CPU and 3.97 times faster than four CPU threads.
Smaller MHD batches remain faster on four CPU threads. Maximum saved-state
CPU/GPU differences were 1.87e-9 m and 2.55e-11 m/s. Magnetic energy and
timestep-convergence checks also passed. These timings do not demonstrate
a performance improvement caused by the patch upgrade; small differences
between runs include ordinary system and timing variability.

## Previous TestParticle 0.24.0 measurements (2026-10-08)

The timing tables below were measured with 0.24.0. The environment was later
updated to 0.24.1, which exports `solve` and tolerates floating-point roundoff
in the spherical GPU grid boundary checks. The old timings are retained as
0.24.0 results, not relabeled as measurements of 0.24.1.
After updating to 0.24.1, all 503 MarsTP tests passed. A 100-particle MHD
CPU/GPU smoke test (O2+, 1 s, dt=0.001 s, saveat=0.1 s, Float64) passed,
with maximum saved-state differences of 9.32e-10 m and 7.95e-12 m/s.

The old raw records remain local in `outputs/gpu_benchmark_20261008_140053/`.
CPU: Intel Core i7-14700KF, four Julia threads. GPU: RTX 4070 SUPER, 12 GB,
driver 591.86. CUDA.jl 6.4.1, KernelAbstractions 0.9.44, TestParticle 0.24.0,
BorisPushers 0.1.1, Julia 1.12.6. No other Julia tests ran during this final
benchmark. An earlier exploratory run is retained separately.

MarsTP MHD forward tracing, 1000 steps per particle:

| Particles | CPU serial (s) | CPU 4 threads (s) | GPU (s) | Serial/GPU | Threads/GPU |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.054061 | 0.028763 | 0.058034 | 0.93 | 0.50 |
| 1000 | 0.284286 | 0.063661 | 0.074031 | 3.84 | 0.86 |
| 10000 | 2.464084 | 0.339063 | 0.092972 | 26.50 | 3.65 |

Analytic uniform-field tracing, also 1000 steps per particle:

| Particles | CPU serial (s) | CPU 4 threads (s) | GPU (s) |
| ---: | ---: | ---: | ---: |
| 100 | 0.009695 | 0.005095 | 0.014334 |
| 1000 | 0.099646 | 0.038148 | 0.024326 |
| 10000 | 1.205213 | 0.390866 | 0.030210 |

Small MHD batches are dominated by field preparation, transfers and output.
At 10000 particles the GPU is 26.5 times faster than the serial API and 3.65
times faster than the four-thread kernel. These are warmed end-to-end ratios
for this short forward example, not kernel-only or full VDF speed estimates.
The MHD grid has shape (120,72,144), and the maximum CPU/GPU discrepancy over
saved states was 1.87e-9 m in position and 2.55e-11 m/s in velocity.
In the uniform magnetic field, relative energy error was 3.11e-15; halving
dt from 0.02 to 0.01 s reduced the analytic endpoint position error from
32.78 to 8.19 m, consistent with second-order convergence.
