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

GPU currently covers fixed-step forward tracing. Adaptive tracing and the
source-accumulating backward VDF remain on CPU. The GPU ensemble kernel has no
per-particle boundary callback. Use a short interval known to stay inside the
field domain. Nonfinite states or solver failures raise an error in MarsTP;
there is no field extrapolation or missing-value replacement. Production
escape/impact statistics need a GPU event/termination implementation before
using this path for long traces. A speed ratio from this short example does
not establish acceleration of the complete backward VDF pipeline.

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
