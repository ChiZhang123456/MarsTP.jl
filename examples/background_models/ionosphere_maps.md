# MHD O₂⁺ flux and ion temperature maps

![O2+ flux and temperature at 200 and 400 km](images/flux_Ti_200_400km.png)

## Definitions

Rows correspond to altitudes of 200 and 400 km. Columns show the scalar flux `n*norm(Ui)` (cm⁻² s⁻¹) and ion temperature Ti (K). Rm=3390 km. Latitude is 90 degrees minus colatitude and longitude is atan2(y,x), using the MHD Cartesian axes rather than geographic coordinates.

The flux uses the full bulk speed. It is neither the signed normal flux through a sphere nor a thermal flux. Density, temperature and the three Cartesian velocity components are interpolated first; flux is then calculated and converted from m⁻² s⁻¹ to cm⁻² s⁻¹ by division by 10⁴.

Both altitudes share one logarithmic color scale per column. The figure uses turbo, Arial, a 183×145 mm layout, and 350 dpi PNG output.

## Code and reproduction

`sample_ionosphere_maps.jl` samples the two surfaces using MarsTP ionosphere interpolation. `plot_ionosphere_maps.py` runs the sampler and draws the figure.

```sh
python examples/background_models/plot_ionosphere_maps.py
# Replot an existing directory containing source.csv:
python examples/background_models/plot_ionosphere_maps.py outputs/ionosphere_maps_<timestamp>
```

Run from the repository root. A new run saves CSV samples, logs and metadata in a new output directory. Figures are written to `examples/background_models/images/` and may replace images with the same name.

Use the project Julia environment and Python with NumPy, Matplotlib and Arial. The original run used Julia 1.12.6 and TestParticle 0.23.3. No external Mars plotting package is required.

The input `data/mars_fields_spherical_from_dat.vts` must provide `n_O^2^p [m^-3]`, `T_O^2^p [K]` and `U_O^2^p [m/s]` with the MarsTP grid conventions. Sampling uses `load_ionosphere_source` and `ionosphere_properties`, not MAT outflow rates.

## Sampling and invalid pole values

Longitude spans −180 to 180 degrees and latitude spans −90 to 90 degrees, both at 2-degree spacing (91×181 points per altitude). Linear interpolation is applied separately to Cartesian velocity components. At the 200 km boundary, a 1e-8 m inward offset avoids floating-point queries outside the grid.

The 32,942 records were checked for count, finiteness, positive temperature, nonnegative flux and consistency with `flux=n*sqrt(ux²+uy²+uz²)`.

The south pole contains temperatures near 1e-10 K and zero speeds. Its entire row is displayed in gray and excluded from color limits. The 91 zero-flux samples at 400 km occur in that row. Original values remain in `source.csv` and metadata. This display mask does not repair pole values used in transport calculations; the original grid construction needs independent verification. Color limits cover all unmasked values without percentile clipping.

## Related backtracing model

`BacktraceConfig(ionosphere_altitude_km=400.0)` defines a traversable source surface while the absorbing inner boundary remains at 200 km. Each transverse crossing contributes `F*g/abs(v dot er)`, with `F=n*norm(Ui)` and normalized drifting Maxwellian g. Volume sources continue to accumulate after crossing. See [the derivation and units](../backward_tracing/backtracing_derivation.md).
