# MarsTP.jl

Clean Mars test-particle tracing package for O2+ work with MHD fields, prepared
source rates, GITM neutral temperature, and detector velocity-distribution
analysis.

## Reproducible figure examples

See [examples/README.md](examples/README.md) for the 800 km O2+ trajectory panels
and the dayside/nightside electric-work maps, including Julia/Python source,
the two PNG figures, prerequisites, run commands, and interpretation notes.

## Data

Required runtime inputs live in `data/`:

- `amps_sph.mat`
- `gitm_sph.mat`
- `O2plus_source_rates.mat`

The MHD `.vts` file is intentionally not tracked in Git. Put the converted file
at `data/mars_fields_spherical_from_dat.vts` before running tracing.

## Basic Use

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

`examples/hemisphere_work.jl` demonstrates `Electric_field_work_profile` on a
reproducible 1000-particle, 800 km O2+ shell release. It uses four
threads when launched by `examples/plot_hemisphere_work.py`, starts at rest, and stops
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

## MHD ionosphere boundary in backtracing

`BacktraceConfig` includes `ionosphere_altitude_km` (200 to 800 km, default
200 km) and `include_ionosphere=true`. The detector must be strictly above the
selected shell and below the outer field boundary:

```julia
cfg = BacktraceConfig(ionosphere_altitude_km=400.0, dt=-0.05)
result = run_backtrace_vdf(cfg)
```

At the first backward entry into this shell, the boundary VDF uses MHD O2+
density `n_O^2^p [m^-3]`, ion temperature `T_O^2^p [K]`, and Cartesian bulk
velocity `U_O^2^p [m/s]` interpolated at that altitude and angular position:
`f = n/(pi^(3/2)*vth^3)*exp(-sum(abs2,v-Ui)/vth^2)`, where
`vth=sqrt(2*kB*Ti/m)`. No radial bulk-velocity clipping or hemisphere
renormalization is applied. Crossing the source shell ends the trajectory.
Raising the shell samples MHD moments at the new height, rather than rescaling
a 200 km source map. Field-grid boundaries remain unchanged.

The scalar flux diagnostic is `n*norm(Ui)` in m^-2 s^-1. It is not a signed
normal flux or a one-way thermal flux through the sphere. The density-normalized
boundary VDF is added once, without a flux, timestep or surface-area multiplier.

```julia
source = load_ionosphere_source(; altitude_km=400.0)
properties = ionosphere_properties(source, SA[1.0, 0.0, 0.0])
# Position specifies direction; sampling is on the configured shell.
# properties.n, properties.Ti, properties.Ui, properties.flux
boundary = ionosphere_distribution(source, SA[1.0, 0.0, 0.0], SA[1000., 0., 0.])
# boundary.f: s^3 m^-6
```

The original volume production model (MAT production density and a zero-drift
neutral-temperature Maxwellian) remains outside the shell. `f2d_volume` and
`f2d_ionosphere` separately store the Vy-integrated contributions; their sum is
`f2d_xz`, all in s^2 m^-5. `include_ionosphere=false` disables the boundary VDF
and uses the 200 km inner termination surface.

Volume quadrature is trapezoidal at accepted step endpoints, including the
fractional boundary segment. Boris staggered velocities are synchronized before
VDF evaluation. Segment-sphere intersections stop before out-of-domain field
evaluation. Invalid MHD moments or source samples raise an error rather than
being replaced with zero. Nonfinite trajectories have a separate status.
Status-count indices: 1 unused, 2 time limit, 3 inner/ionosphere shell,
4 nonfinite trajectory, 5 outer boundary.

`test/backtrace_ionosphere.jl` checks the Maxwellian peak, nonradial flux,
allowed heights, both Boris solvers, straight-line crossings, fractional-step
quadrature, outer exits and time limits. `scripts/smoke_ionosphere.jl` checks the
local MHD input and a single-particle VDF. Step convergence must still be checked
for each scientific detector/velocity-grid configuration.
