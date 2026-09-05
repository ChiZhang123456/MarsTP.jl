# MHD Conversion

The MHD field file used by tracing is:

```text
data/mars_fields_spherical_from_dat.vts
```

It is generated from:

```text
3d__mhd_1_n00010000.dat
```

The `.dat` and `.vts` files are intentionally not stored in GitHub. Keep them as
local or Perlmutter runtime data.

The converted `.vts` must contain:

- `B_Field [T]`
- `E_Total [V/m]`
- `E_conv [V/m]`
- `E_hall [V/m]`
- MHD ion moments for H+, O+, O2+, and CO2+
