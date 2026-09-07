# Mars boundaries and trajectory visualization

Python visualization with NumPy, Matplotlib, h5py and Python 3.11+. No `py_space_zc` installation is required. Run commands from the repository root.

## Trajectory plots

`plot_trajectory.py` reads saved trajectories and uses `mars.py` to draw Mars, the bow shock (BS) and magnetic pileup boundary (MPB). Output is PNG.

```powershell
# Forward JLD2 batches: select O2+ and draw three projections.
python src/visualization/plot_trajectory.py outputs/my_run --species O2+ --count 5000 --output outputs/forward.png

# Backward JSONL: positions in Rm and times in s.
python src/visualization/plot_trajectory.py outputs/my_backtrace/trajectories.jsonl --position-unit Rm --species O2+ --planes XZ XY YZ 3D --output outputs/backward.png

# Select species and computational-particle IDs.
python src/visualization/plot_trajectory.py outputs/my_run --species O2+ --particle-ids 12 45 91 --planes 3D --output outputs/selected.png
```

`--species` selects the ion species; `--particle-ids` selects saved particle numbers. Species labels are read from per-particle fields, file fields or adjacent `metadata.toml`/`metadata.json`. Mixed-species data can label individual records. For old data without labels, `--assume-species O2+` explicitly declares the species. It does not override existing labels or transform particles into another species.

`--count` limits displayed trajectories (default 5000), using uniform reservoir sampling with a fixed seed. `--max-points` limits displayed points per trajectory (default 300), retaining both endpoints. A file and particle ID together identify a path; pass a single file to narrow the selection.

### File formats

- MarsTP JLD2 batches: `p<ID>/state` is N×7 through h5py, with t,x,y,z,vx,vy,vz in s, m and m/s. HDF5 files with the same layout are supported. Increasing and decreasing time sequences are retained for forward and backward paths.
- JSONL: each object requires `id`, N×3 `points`, and N `times_s` values. Optional fields include `species`, `coordinate_system` and `position_unit`. Position units are `m`, `km` or `Rm`, optionally specified with `--position-unit`.
- PSD-only or endpoint-only files cannot reconstruct full paths; save trajectories during integration.

The default Mars radius is 3,390,000 m, configurable with `--rm-m`. Plot coordinates use Rm. For JLD2, only selected trajectories and display points are read into memory.

### Python interface

```python
import sys
sys.path.insert(0, "src")
from visualization import plot_trajectory

fig, axes, records = plot_trajectory(
    "outputs/my_run", species="O2+", count=100,
    planes=("XZ", "3D"), output="outputs/trajectory.png")
```

## Mars, BS and MPB

In two dimensions, `plot_mars` defaults to the bundled [Mars image](mars_globe_true_color.png). No image path is required. This image was copied unchanged from the user-provided local `py_space_zc.maven` resource. Use `texture=False` for a solid disk or `texture_path=...` for a custom image.

Three-dimensional axes display a sphere. The bundled disk photograph is not a global longitude/latitude texture and is not wrapped onto the sphere.

`bs_mpb` uses the same conic parameters as `py_space_zc.maven.bs_mpb`. Let x,ρ,r,x₀,L be in Rm, θ in rad and eccentricity ε dimensionless:

$$
r=\frac{L}{1+\epsilon\cos\theta},\qquad
x=x_0+r\cos\theta,\qquad \rho=r\sin\theta.
$$

| Boundary | x₀ (Rm) | L (Rm) | ε | Region |
| --- | ---: | ---: | ---: | --- |
| BS | 0.600 | 2.081 | 1.026 | All |
| Dayside MPB | 0.640 | 1.080 | 0.770 | x≥0 |
| Nightside MPB | 1.600 | 0.528 | 1.009 | x<0 |

Only positive-radius branches are used. Rotation about X generates the three-dimensional surfaces; XY and XZ show symmetric sections. For YZ, `--x-slice 0` displays the boundary section at X=0 Rm. No YZ boundary is drawn by default.

The model assumes +X sunward. For other coordinate conventions, transform positions first or use `--no-boundaries`.

```python
import matplotlib.pyplot as plt
from visualization import plot_mars_context, bs_mpb, plot_mars

fig, ax = plt.subplots()
plot_mars_context(ax, plane="XZ", xmin=-4.2)
ax.set(xlim=(-4.2, 4.2), ylim=(-4.2, 4.2), aspect="equal")

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
bs_mpb(ax, plane="3D")
plot_mars(ax)
ax.set_box_aspect((1, 1, 1))
```
