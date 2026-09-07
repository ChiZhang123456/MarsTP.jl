# MarsTP.jl

Clean Mars test-particle tracing package for O2+ work with MHD fields, prepared
source rates, GITM neutral temperature, and detector velocity-distribution
analysis.

## Reproducible figure examples

See [examples/README.md](examples/README.md) for the 800 km O2+ trajectory panels
and the dayside/nightside electric-work maps, including Julia/Python source,
the two PNG figures, prerequisites, run commands, and interpretation notes.

Saved forward trajectories can be analyzed with `forward_psd_saved`, using the
same accumulator as `forward_psd`. `write_trajectory_batch` provides bounded
batch output for large ensembles. See the [Monte Carlo API and reproducible
three-probe example](examples/forward_tracing/monte_carlo_forward_tracing/README.md).

## Data

Required runtime inputs live in `data/`:

- `amps_sph.mat`
- `gitm_sph.mat`
- `O2plus_source_rates.mat`

The MHD `.vts` file is intentionally not tracked in Git. Put the converted file
at `data/mars_fields_spherical_from_dat.vts` before running tracing.

## Basic Use

Monte Carlo Maxwellian initial states and physical weights are available via
`sample_maxwellian_source`. See [weight definitions and units](examples/forward_tracing/monte_carlo_weights.md)
and the [local MHD source example](examples/forward_tracing/maxwellian_source.jl).
For a finite-volume detector, see the [3D phase-space density derivation](examples/forward_tracing/monte_carlo_forward_tracing/README.md),
including residence-time weighting, units, and the distinction from source density weights.
Use `forward_psd(solutions; detector_m, side_m, vlim, vgrid, species,
rate_weights_s, option="3D")` to compute the detector PSD from saved forward
trajectories. `vgrid` is the bin count per axis; `vlim` is in m/s by default.
Options `"Vx-Vy"`, `"Vy-Vz"`, and `"Vx-Vz"` return integrated 2D distributions.
See the guide above for units, source-rate normalization, and a complete example.

```julia
using MarsTP

cfg = BacktraceConfig(
    detector_Rm = SA[-1.5, 0.0, 1.0],
    vx_min_kms = -500,
    vx_max_kms = 500,
    vy_min_kms = -500,
    vy_max_kms = 500,
    vz_min_kms = -500,
    vz_max_kms = 500,
    dv_kms = 1,
    dvy_kms = 1,
    solver = :boris,
    dt = -0.2,
    stream_vy = true,
)

result = run_backtrace_vdf(cfg)
```

The returned `f2d_xz` is the 3D VDF integrated over `Vy` using `sum(f) * dvy`.

## Per-particle electric work

Build the electric-field interpolators once, then reuse them for saved in-memory
trajectories (Cartesian positions in m, velocities in m/s, times in s):

```julia
itp = build_field_work_interpolators()
summary = field_work(sol, itp; species = "O2+")
profile = Electric_field_work_profile(sol, itp; species = "O2+")
# profile.t, profile.conv_eV, profile.hall_eV, profile.total_eV
# profile.delta_kinetic_eV, profile.energy_residual_eV
# profile.power_conv_eV_s, profile.power_hall_eV_s

# solutions: TestParticle ensemble or vector returned by trace_forward
per_particle = particle_field_work(solutions, itp; threaded = true)
histories = particle_field_work(solutions, itp; profiles = true, threaded = true)
```

Each component is signed work `q ∫ E_component ⋅ v dt` evaluated on the **same
total-field trajectory**. Positive work adds energy; negative work removes it.
Power is in eV/s, cumulative profiles in eV; `field_work` also returns J.
The work functions do not run trajectories or save files. Single-particle calls
reject ensembles with multiple members to avoid silently using only the first.

The residual `ΔK - W_total` diagnoses numerical energy closure for nonrelativistic
trajectories subject only to electromagnetic forces. Check timestep and saved
sample cadence: quadrature cannot recover missing trajectory structure. Reverse
time integrations preserve signed time intervals. Work fractions do not predict
the result of disabling a field component, because doing so changes the trajectory.
Avoid evaluating interpolators beyond their field domain.

Run the focused tests without the MHD input file:

```sh
julia --threads=2 --project=. test/runtests.jl
```

The former name `field_work_profile` remains available as a compatibility alias for
`Electric_field_work_profile`, with identical arguments and return values.

`examples/forward_tracing/hemisphere_work.jl` demonstrates `Electric_field_work_profile` on a
reproducible 1000-particle, 800 km O2+ shell release. It uses four
threads when launched by `examples/forward_tracing/plot_hemisphere_work.py`, starts at rest, and stops
at the 200 km or 4 Mars-radius boundary (20,000 s guard). Work ends at the last
valid in-domain state. Energy closure above 0.1% triggers timestep refinement;
this is a diagnostic, not a guarantee of full trajectory convergence.

The Python script requires NumPy, Matplotlib and the user's `py_space_zc` library.
Set `HEMISPHERE_WORK_PREVIEW` to customize the output PNG path. It plots six
XZ panels with dayside and nightside particles in separate rows and convection,
Hall and total work in three columns. Each row uses its own symmetric-log color
scale. See [the example documentation](examples/README.md). Only the preview image is
saved; trajectory data pass through memory. Threaded output is not ordered, so
each record carries its original particle ID.

## Perlmutter

Use `scripts/perlmutter/run_backtrace_vdf.slurm` as the starting point. Copy or
download the MHD `.vts` file into `data/` on Perlmutter before submitting.

## MHD ionosphere thin-sheet source in backtracing

The source is a penetrable shell, **400 km by default**. The absorbing inner
boundary is always **200 km**, independent of source altitude (allowed range:
200 to 800 km). Backtracing continues below the shell and accumulates volume
production until it reaches 200 km, the outer boundary, or the time limit.

```julia
cfg = BacktraceConfig(ionosphere_altitude_km=400.0, dt=-0.05)
result = run_backtrace_vdf(cfg)
```

MHD O2+ density, ion temperature and Cartesian bulk velocity define
`F = n*norm(Ui)` in m^-2 s^-1 and a normalized drifting Maxwellian `g` in
s^3 m^-3. Each transverse crossing adds `F*g/abs(dot(v,er))` in s^3 m^-6,
where `v` is the traced particle velocity. All resolved crossings count, in
either direction. The source shell does not terminate or scatter particles.
This replaces the previous density-normalized boundary-VDF prescription.

Volume production still uses the MAT volume rate and a zero-drift
neutral-temperature Maxwellian. The saved `f2d_volume` and `f2d_ionosphere`
are integrated over Vy in s^2 m^-5; `f2d_xz` is their sum.
`ionosphere_crossings` counts sheet encounters per Vx/Vz cell across sampled Vy.
`model="thin_shell_source_v1"` distinguishes these files from older results.

`load_ionosphere_source`, `ionosphere_properties` and `ionosphere_distribution`
expose the MHD shell moments. The last helper returns both normalized `g` and
`f=n*g`; the transport source uses `flux*g`, not `flux*f`.
`include_ionosphere=false` disables only the sheet source, retaining 200 km
termination. Detectors below the source sheet are supported. Detectors exactly
on the ideal sheet and unresolved grazing crossings raise errors; no arbitrary
radial-speed floor is applied. A finite-thickness model would be required to
regularize these cases.

See the [full derivation and unit audit](examples/backward_tracing/backtracing_derivation.md)
for the transport equation, delta-function change of variables, volume
quadrature, crossing treatment, unit table, background assumptions and limits.
Run `julia --project=. test/runtests.jl` for analytic tests and
`julia --project=. scripts/smoke_ionosphere.jl` for the local MHD smoke test.

## Trajectory visualization

[src/visualization](src/visualization/README.md) provides standalone Python plots of Mars, BS/MPB and two- or three-dimensional trajectories, with species and particle-ID selection from forward or backward trajectory files.
