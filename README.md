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

## Perlmutter

Use `scripts/perlmutter/run_backtrace_vdf.slurm` as the starting point. Copy or
download the MHD `.vts` file into `data/` on Perlmutter before submitting.
