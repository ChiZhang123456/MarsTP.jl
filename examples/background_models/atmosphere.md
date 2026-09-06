# GITM and AMPS atmosphere example

![Atmosphere profiles and maps](images/atmosphere_gitm200km_amps500km.png)

Run from the repository root with Python, NumPy, SciPy and Matplotlib installed:

```sh
python examples/background_models/plot_atmosphere.py
```

Arial must be available. Inputs are the tracked `data/gitm_sph.mat` and
`data/amps_sph.mat`. The PNG (300 dpi) is saved in images/; source arrays and QA metadata are saved beside the script. Rerunning replaces these example outputs. Only the current PNG figure is tracked; no PDF or SVG is exported.

Panels, in row order:

1. Number density versus altitude at longitude = latitude = 0 degrees.
2. GITM neutral temperature versus altitude at the same location.
3. GITM CO2 density at 200 km.
4. GITM O density at 200 km.
5. GITM neutral temperature at 200 km.
6. AMPS hot O density at 500 km.

All maps use turbo with black labeled contours and independent colorbars.
Density maps use logarithmic normalization; temperature uses linear normalization.
The altitude range shown in profiles is 100 to 600 km. Mars radius is 3390 km.
Density inputs are interpreted in SI as consumed by MarsTP and displayed in
cm^-3. Temperatures are in K. Coordinates retain the input longitude/latitude
definitions; no additional geographic or MSO transformation is imposed.

GITM input covers only 100 to 220 km. Above 220 km this **plotting example** uses
the existing constant-temperature, constant-gravity exponential prescription,
but explicitly removes its artificial density cutoff at 20 scale heights:

```text
T0 = T(219995 m)
H = kB * T0 / (mass * 3.71 m/s^2)
n(h) = n(219995 m) * exp(-(h - 220000 m) / H)
T(h) = T0
```

Species masses follow MarsTP: 15.999 and 44.01 times the proton mass.
Dashed profile segments identify this extension. The three GITM maps are native 200 km slices, checked directly against the input arrays. The d panel uses nO, not nCO2 or nO_hot. This example does not change
the tracing library's `neutral_properties` implementation.

AMPS provides hot O density only, with no temperature. Its 500 km slice is
native in altitude. Original 0/360-degree longitude endpoint differences are
preserved. No missing values are replaced and no smoothing is applied.

QA checks validate dimensions, coordinate relations, monotonic axes, finite
positive inputs, interpolation at a grid point, and positive finite map values.
`metadata.json` records input hashes, dependency versions, ranges and zero counts.
`atmosphere_source_data.npz` contains profile densities in m^-3, temperatures in K,
map densities in cm^-3, map temperatures in K, and coordinate arrays.
