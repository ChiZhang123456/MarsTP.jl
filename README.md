# MarsTP.jl

Clean Mars test-particle tracing package for O2+ work with MHD fields, prepared
source rates, GITM neutral temperature, and detector velocity-distribution
analysis.

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
profile = field_work_profile(sol, itp; species = "O2+")
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

## Perlmutter

Use `scripts/perlmutter/run_backtrace_vdf.slurm` as the starting point. Copy or
download the MHD `.vts` file into `data/` on Perlmutter before submitting.
