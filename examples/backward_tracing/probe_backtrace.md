# Backtracing 5000 O₂⁺ particles from a fixed detector

![XZ trajectories](images/probe_backtrace_xz_5000.png)

## Code and reproduction

- `random_probe_backtrace.jl`: trajectory integration, boundary termination and time-step comparison.
- `plot_probe_backtrace.py`: runs Julia, saves paths and diagnostics, and plots XZ.
- `images/probe_backtrace_xz_5000.png`: example figure.

Run from the repository root:

```sh
python examples/backward_tracing/plot_probe_backtrace.py
# Replot saved paths:
python examples/backward_tracing/plot_probe_backtrace.py --data outputs/<run_id>/trajectories.jsonl
```

Python requires NumPy, Matplotlib, h5py and Arial. Mars and boundary plotting use `src/visualization` and its bundled image. Julia must be on PATH and use Project.toml/Manifest.toml. The input `data/mars_fields_spherical_from_dat.vts` is not distributed with the code. The original run used Julia 1.12.6 and TestParticle 0.23.3.

`PARTICLE_COUNT` overrides the default 5000; in replot mode it must match the record count. Diagnostics and initial velocities are saved in a new `outputs/probe_backtrace_<timestamp>/` directory. New integrations also save `trajectories.jsonl`, with particle ID, initial velocity, negative times, three-dimensional positions and convergence diagnostics. Replot mode reads the original file. Images are written to `examples/backward_tracing/images/` and replace matching filenames.

## Initial conditions and model

| Parameter | Value |
| --- | --- |
| Detector | (0,0,2 Rm), Mars-centered coordinates |
| Mars radius | 3390 km |
| Species | Nonrelativistic O₂⁺; mass and positive charge from TestParticle SpeciesDict |
| Initial speed | Uniform from 10 to 200 km/s |
| Initial direction | Uniform angle in XZ, initial Vy=0 |
| Random generator | Xoshiro(20260905) |
| Integration | Boris from t=0 toward negative time, dt=−0.05 s |
| Comparison step | −0.1 s for all 5000 particles |
| Saved cadence | 0.5 s, plus boundary intersections |
| Inner boundary | 200 km altitude |
| Outer boundary | 4 Rm from the center |
| Maximum lookback | 500 s |
| Fields | Static total E and B, with MarsTP spherical-grid interpolation and SI conversion |

Velocity and charge retain their physical signs; negative time implements backtracing. All three position and velocity components evolve. Collisions, gravity, chemical weights and feedback are omitted. A line-segment intersection locates the final spherical-boundary point without extrapolating fields outside the domain.

## Figure

The XZ projection uses equal axes and all 5000 paths. Turbo indicates initial speed, and the red star marks the detector. Arial is used for nonmathematical text. BS and MPB reference curves are drawn as dashed and dotted lines, with a radius-1 Mars image.

Coordinates follow the input Cartesian axes. The MSO/MSE designation has not been independently established from metadata. Reference boundaries assume +X sunward. Since the XZ projection includes every Y, a path crossing the displayed disk need not intersect the planet in three dimensions.

## Recorded checks

In the original run, 4989 particles reached the outer boundary and 11 reached the inner boundary; none reached 500 s. Count, finite positions, initial speeds and Vy, decreasing time and boundary intersections were checked.

Reducing dt from 0.1 to 0.05 s preserved all termination categories. The largest same-time separation was 0.001034 Rm (3.51 km), the largest endpoint difference was 0.000865 Rm (2.93 km), and the largest termination-time difference was about 0.0100 s. These are results of this comparison, not general error bounds. The figure was visually checked.
