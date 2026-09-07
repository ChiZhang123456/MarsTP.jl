# Local power and cumulative work along backtraced paths

![Path power](images/probe_path_power_xz_5000.png)

Columns show convective, Hall and total electric-field contributions. The first row is local power (eV/s), and the second is cumulative work (eV). The example reuses 5000 O₂⁺ initial velocities: uniform speed from 10 to 200 km/s, uniform direction in XZ, initial Vy=0 and seed 20260905. The detector is at (0,0,2 Rm), with Rm=3390 km. Motion is three-dimensional; only its XZ projection is shown.

## Local power

Charge q is in C, electric field E in V/m and velocity v in m/s. The following power is in J/s; divide by 1.602176634e-19 J/eV for eV/s:

$$
P_j=q\mathbf E_j\cdot\mathbf v,\qquad j=\mathrm{conv,Hall,total}.
$$

Positive values (red) indicate energy gain in physical forward time; negative values (blue) indicate loss. The magnetic force does no work. Negative integration steps do not reverse the charge, physical velocity or power sign. For a backward step dt<0, physical forward-time work is `-q*dot(E_mid,v_mid)*dt`, converted from J to eV.

Work is accumulated at every −0.05 s integration step. Display segments span 0.5 s, with the actual duration used for the final partial segment. Their color is the segment work divided by duration. Velocity is not estimated from sparse saved positions.

Coolwarm and SymLogNorm use a linear interval of ±1 eV/s and logarithmic scaling outside it, shared across all three panels over ±1000 eV/s. Overlapping transparent paths are not spatial averages. Mars and boundaries use `src/visualization`.

All field terms are evaluated along the same trajectory driven by the total field. Local power differs from the [net-energy map](energy_gain.md): a particle may gain and later lose energy along one path.

## Cumulative work

Let the detector time be 0, the past endpoint be t_b<0, and t be a time along the path, all in s. With q in C, E in V/m and v in m/s, C_j and ΔK_j below are in J (converted to eV for plotting):

$$
C_j(t)=\int_{t_b}^{t}q\mathbf E_j\cdot\mathbf v\,dt',\qquad
C_j(t_b)=0,\qquad C_j(0)=\Delta K_j.
$$

Work accumulates from the past endpoint toward the detector. That endpoint is a boundary or time limit, not necessarily the particle's birthplace. The code saves B(t), the forward-time work accumulated while integrating backward from the detector, then uses `C(t)=DeltaK-B(t)`. Each segment is colored by the mean of its endpoint C values.

The second row uses coolwarm and SymLogNorm, linear within ±1 eV. Its three panels share a color range, independently of the first row. Signed values are retained. Exact per-particle totals are in `particles.csv`.

## Run

Use the project Julia environment and Python with NumPy, Matplotlib, h5py and Arial:

```sh
julia --startup-file=no --compiled-modules=existing --threads=1 --project=. examples/backward_tracing/probe_path_power.jl outputs/path_power_new_run 5000
python examples/backward_tracing/plot_path_power.py outputs/path_power_new_run
```

Replace 5000 with 3 for a small run. The original position-only paths do not contain velocity, so this example reintegrates their initial conditions. Boundaries are 200 km and 4 Rm, with a 500 s limit. The 400 km source surface does not terminate these unweighted paths. Only PNG figures are exported.

The original `segments.csv` contains 480,389 display segments and remains local. Records contain particle ID, XZ endpoints, duration, three mean-power terms and cumulative backward-scan work. Particle totals (`path_power_particles.csv`, generated run output) and plot checks (`path_power_qa.json`, generated run output) are generated with the run.

## Recorded checks and interpretation

Initial velocities matched the original example. Eleven particles reached the inner boundary and 4989 the outer boundary, without time limits or nonfinite states. The maximum energy-closure residual was 0.2063 eV. Relative to `max(abs(work),abs(delta K),1 eV)`, the 99th percentile and maximum residuals were 0.00170% and 0.1672%. Total power differed from convective plus Hall power by at most 2.49e-5 eV/s. The original 210 tests covered observer segment sums, both Boris solvers and signed power.

These paths illustrate mechanisms without source-PSD weighting. The 0.5 s averaging can hide shorter sign changes; energy closure does not establish convergence of local power structure. Cumulative endpoints matched particle totals within 3.64e-11 eV for all 5000 particles. Initial conditions, termination and total work were preserved when the cumulative row was added; its shared range is ±100000 eV.
