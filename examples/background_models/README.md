# Background models: MHD, GITM and AMPS

This directory provides maps of electromagnetic fields, ion moments and neutral atmosphere parameters. Physical inputs remain in the repository `data/` directory.

## MHD O₂⁺ flux and temperature

![MHD flux and temperature](images/flux_Ti_200_400km.png)

Maps at 200 and 400 km show `n*norm(Ui)` on the left and ion temperature on the right. Spherical angles are defined relative to the model axes. Anomalous south-pole samples are masked for display. See [map definitions](ionosphere_maps.md).

```sh
python examples/background_models/plot_ionosphere_maps.py
```

The plotting script calls [sample_ionosphere_maps.jl](sample_ionosphere_maps.jl). To reuse saved samples:

```sh
python examples/background_models/plot_ionosphere_maps.py outputs/ionosphere_maps_<timestamp>
```

## GITM and AMPS

![GITM and AMPS](images/atmosphere_gitm200km_amps500km.png)

The figure includes GITM neutral density and temperature profiles and 200 km maps, plus AMPS hot oxygen profiles and a 500 km map. AMPS does not provide temperature. See [atmospheric assumptions](atmosphere.md).

```sh
python examples/background_models/plot_atmosphere.py
```

Figures are saved in `images/`. Generated `atmosphere_source_data.npz` and `metadata.json` remain local.

## Inputs and dependencies

- MHD: `data/mars_fields_spherical_from_dat.vts`, fields and ion moments; not distributed with the repository.
- GITM: `data/gitm_sph.mat`, neutral atmosphere parameters.
- AMPS: `data/amps_sph.mat`, hot oxygen density.
- Python: NumPy, Matplotlib, SciPy for atmospheric plots, and Arial.
- Julia: the repository Project.toml/Manifest.toml environment, with `julia` on PATH.

Run from the repository root. Scripts do not install dependencies automatically. Return to [examples](../README.md).
