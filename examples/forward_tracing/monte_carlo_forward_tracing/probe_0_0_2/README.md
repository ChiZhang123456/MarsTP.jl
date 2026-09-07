# O₂⁺ detector at (0, 0, 2) Rm

![Velocity-integrated VDF](probe_psd_projections.png)

The cubic detector is centered at **(0, 0, 2) Rm**, with side length **0.2 Rm=678 km**. Residence times are extracted from saved trajectories of the 500 km source simulation.

The three-dimensional PSD uses **5 km/s bins from −500 to 500 km/s** on each axis and has units **s³ m⁻⁶**. Integrating the omitted axis gives `fxy=sum_z(f3d)*5000 m/s` and `fxz=sum_y(f3d)*5000 m/s`, in **s² m⁻⁵**. Plot limits are ±300 km/s, with turbo and logarithmic normalization, without smoothing.

| Quantity | Value |
| --- | ---: |
| Density (m⁻³) | 1002.16762 |
| Density (cm⁻³) | 0.00100216762 |
| Independent particles intersecting the detector | 115 |
| Residence-weight effective sample size | 9.708 |
| Residence segments | 3,328 |
| Face crossings | 230 |
| Nonzero three-dimensional bins | 487 |

Local `probe_psd_sparse.npz` stores nonzero three-dimensional values, zero-based indices and integrated projections for the 200³ grid. Unlisted bins represent zero sample contributions. `analysis_summary.json` stores detailed statistics. These generated files are not included with the GitHub example.

See [the main guide](../README.md) for calculation steps. To plot this detector after analysis:

```powershell
python examples/forward_tracing/monte_carlo_forward_tracing/plot_library_psd.py outputs/my_probes/probe_0_0_2
```

Coordinates follow the native MHD Cartesian axes. The maximum flight age is 500 s; steady-state convergence has not been established. Color scales are normalized separately between detectors. Fine velocity structure requires sufficient independent weighted samples.
