# Data Files

This directory should contain only prepared runtime inputs.

| File | Description |
|---|---|
| `amps_sph.mat` | AMPS hot atomic oxygen spherical grid. |
| `gitm_sph.mat` | GITM neutral atmosphere spherical grid. |
| `O2plus_source_rates.mat` | Prepared O2+ production and outflow rates. |

The tracing code also expects the converted MHD field file at:

```text
data/mars_fields_spherical_from_dat.vts
```

This `.vts` file is not uploaded to GitHub because it is larger than GitHub's
normal file-size limit. Generate it from `3d__mhd_1_n00010000.dat` or copy it to
Perlmutter manually.

Important fields in `O2plus_source_rates.mat`:

| Variable | Unit |
|---|---|
| `r_m`, `theta_rad`, `phi_rad` | m, rad, rad |
| `O2plus_production_density_m3_inv_s_inv` | m^-3 s^-1 |
| `O2plus_production_rate_s_inv` | s^-1 |
| `O2plus_outflow_rate_s_inv` | s^-1 |
| `cell_volume_m3` | m^3 |
| `surface_area_m2` | m^2 |
