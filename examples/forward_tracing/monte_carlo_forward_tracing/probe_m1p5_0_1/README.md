# O₂⁺ detector at (−1.5, 0, 1) Rm

![Velocity-integrated VDF](probe_psd_projections.png)

The cubic detector is centered at **(−1.5, 0, 1) Rm**, with side length **0.2 Rm=678 km**. Residence times are extracted from saved trajectories of the 500 km source simulation.

The three-dimensional PSD uses **5 km/s bins from −500 to 500 km/s** on each axis and has units **s³ m⁻⁶**. Integrating the omitted axis gives `fxy=sum_z(f3d)*5000 m/s` and `fxz=sum_y(f3d)*5000 m/s`, in **s² m⁻⁵**. Plot limits are ±300 km/s, with turbo and logarithmic normalization, without smoothing.

| Quantity | Value |
| --- | ---: |
| Density (m⁻³) | 117721.665 |
| Density (cm⁻³) | 0.117721665 |
| Independent particles intersecting the detector | 2,011 |
| Residence-weight effective sample size | 353.706 |
| Residence segments | 124,374 |
| Face crossings | 4,042 |
| Nonzero three-dimensional bins | 3,285 |

Local `probe_psd_sparse.npz` stores nonzero three-dimensional values, zero-based indices and integrated projections for the 200³ grid. Unlisted bins represent zero sample contributions. `analysis_summary.json` stores detailed statistics. These generated files are not included with the GitHub example.

See [the main guide](../README.md) for calculation steps. To plot this detector after analysis:

```powershell
python examples/forward_tracing/monte_carlo_forward_tracing/plot_library_psd.py outputs/my_probes/probe_m1p5_0_1
```

Coordinates follow the native MHD Cartesian axes. The maximum flight age is 500 s; steady-state convergence has not been established. Color scales are normalized separately between detectors. Fine velocity structure requires sufficient independent weighted samples.
