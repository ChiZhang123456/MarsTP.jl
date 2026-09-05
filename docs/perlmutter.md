# Perlmutter Notes

Before submitting a job:

1. Clone `MarsTP.jl`.
2. Copy `data/mars_fields_spherical_from_dat.vts` into the repository `data/`
   directory. This file is intentionally not tracked by Git.
3. Instantiate the Julia environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Submit a detector backtracing run:

```bash
sbatch scripts/perlmutter/run_backtrace_vdf.slurm
```

The default grid is `Vx,Vy,Vz=-500:1:500 km/s`, with `Vy` streamed one slice at
a time and integrated into `f(Vx,Vz)` using `sum(f) * dvy`.
