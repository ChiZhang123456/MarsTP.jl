# Direction-averaged differential energy flux

`detector_omni_def` computes the differential energy flux (DEF) averaged over all directions. It retains sr⁻¹ and does not assume an isotropic velocity distribution. The nonrelativistic energy is `E = m*norm(v)^2/(2e)` in eV, where v is the particle velocity in the saved coordinate frame, in m/s. No bulk-velocity subtraction is applied.

For steady forward source rates Q_i in particles/s:

$$
\mathrm{DEF}_k = \frac{\sum_i Q_i\int_{\mathbf{x}_i\in V_{\mathrm{det}},\,E_i\in k} E_i v_i\,dt}
{4\pi V_{\mathrm{det}}\Delta E_k}.
$$

The result has units **eV/(m² s eV sr)**. The factor 4π defines the directional average, not an isotropy assumption. This uses the particle speed along the trajectory, not the bulk speed used to assign source rates. Source rates already include Monte Carlo normalization and must not be divided by sample count or total flight age again.

## Directly from forward trajectories

```julia
using MarsTP
edges = 10.0 .^ range(-1, 4; length=51) # 50 energy bins, eV
omni = detector_omni_def(solutions;
    detector_m=[-1.5,0,1]*Rm, side_m=0.2Rm,
    rate_weights_s=Q, energy_edges_eV=edges, species="O2+")
omni.def                  # eV/(m^2 s eV sr)
omni.energy_centers_eV
omni.dn_dE_m3_eV          # m^-3 eV^-1
```

Trajectories use synchronized positions in m, velocities in m/s and increasing time in s. The detector cube is axis-aligned. Each linear saved segment is clipped to the cube and split at its energy-edge crossings, including two crossings when speed reverses. The Ev integral uses adaptive quadrature with relative tolerance 10⁻⁹ on each resulting interval. Unresolved quadrature raises an error. This calculation needs neither a velocity cube nor pitch-angle/gyrophase bins. Check saved-cadence convergence because unresolved orbit curvature cannot be recovered.

## Together with a forward velocity PSD

```julia
r = forward_psd(solutions;
    detector_m=[-1.5,0,1]*Rm, side_m=0.2Rm,
    vlim=500e3, vgrid=100, rate_weights_s=Q,
    energy_edges_eV=edges, species="O2+", option="Vx-Vz")
r.omni_def.def

# Read Q directly from saved particle files without loading the full ensemble.
r = forward_psd_saved("outputs/my_run";
    detector_m=[-1.5,0,1]*Rm, side_m=0.2Rm,
    vlim=500e3, vgrid=100, energy_edges_eV=edges)
```

`ForwardPSDAccumulator` also accepts `energy_edges_eV`; `finish_forward_psd` returns `omni_def`. The DEF uses every saved in-detector velocity, independently of Cartesian PSD `vlim`, `vgrid`, projection, or sparse/dense storage. Without `energy_edges_eV`, `omni_def` is `nothing` and no spectrum is accumulated.

## Backward detector output

```julia
c = BacktraceConfig(energy_edges_eV=edges, stream_vy=true,
    output_file="outputs/backtrace_with_def.jld2",
    checkpoint_file="outputs/backtrace_with_def_checkpoint.jld2")
r = run_backtrace_vdf(c)
r.omni_def.def
```

At each detector velocity node, the backward solver supplies f in s³ m⁻⁶. The shared analysis code accumulates `f*dV*E*v/(4pi*DeltaE)` before integrating away vy. It uses the same rectangular velocity-node weights as the existing backward VDF: `dV=dvx*dvy*dvz`, including endpoint nodes with full weights. This works with `stream_vy=true` and `include_work=false`, without retaining a full 3D PSD. The JLD2 output and checkpoints include `omni_def`; checkpoint values cover only the completed vy slices. A failed backward velocity sample causes an error when DEF is enabled, rather than being treated as zero.

For an existing full 3D PSD, ordered (vx,vy,vz):

```julia
omni = detector_omni_def(f3d, (vx_m_s,vy_m_s,vz_m_s);
    velocity_cell_widths_m_s=(dvx_m_s,dvy_m_s,dvz_m_s),
    energy_edges_eV=edges, species="O2+")
```

Widths may be scalars or per-axis vectors. Each node's density is assigned to its node energy. Converge the velocity grid relative to the energy-bin widths. The computed velocity domain may not contain the entire distribution: missing directions and velocity tails are not extrapolated or rescaled. An already projected 2D PSD is insufficient to reconstruct energy, so it is not accepted.

## Output and interpretation

`def`, `energy_edges_eV`, `energy_centers_eV` (arithmetic midpoints), `energy_widths_eV`, `density_per_bin_m3`, and `dn_dE_m3_eV` share the energy-bin ordering. `density_in_range_m3`, `density_outside_energy_range_m3`, and `density_total_m3` report the density budget within the available input. Edges must be finite, nonnegative and strictly increasing. Bins use [lower,upper), with the final upper edge included. Values outside the requested energy range are reported in the density budget and not redistributed.

For narrow bins, DEF is approximately `E*v(E)/(4pi)*dn/dE`. The trajectory estimator integrates Ev within each bin instead of assuming this approximation. Zero-speed particles contribute density but zero DEF. Energy is in eV for both the numerator and bin width; there is no additional eV-to-joule factor in the final normalization. To express DEF per cm², divide the returned values by 10⁴.

Forward/backward spectra represent different detector measures unless their physical models agree: forward averages over a finite detector volume, whereas backward evaluates a detector point. Compare them only with consistent sources, boundaries, velocity coverage and converged spatial/velocity/time resolution.

## Validation

The existing 420 package assertions and 27 DEF assertions passed with Julia 1.12.6 and TestParticle 0.23.3. DEF checks cover an anisotropic mono-directional beam, exact residence clipping, energy-bin crossings through a velocity reversal, smooth non-collinear velocity variation, zero speed, final energy-edge inclusion, density outside energy coverage, invalid input, and equivalence of direct, streaming and VDF-node estimators.

An eight-velocity-node local-field backward smoke test compared streaming output with full-3D-PSD postprocessing and checked the final JLD2 file and checkpoint. All six comparisons passed. This validates integration and saving, not scientific convergence of a production energy spectrum. Run the unit checks with `julia --project=. test/runtests.jl`.
