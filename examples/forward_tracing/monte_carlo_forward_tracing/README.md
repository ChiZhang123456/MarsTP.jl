# O₂⁺ Monte Carlo sampling, forward tracing and detector PSD

This example samples O₂⁺ source particles, follows their trajectories through MHD fields, and reconstructs velocity distributions from residence times inside a cubic detector.

## 1. Monte Carlo source sampling

Monte Carlo represents a continuous velocity distribution with a finite set of weighted particles. Here, velocities are drawn from a broader Maxwellian to improve sampling of the high-speed tail. Importance weights recover the physical distribution.

### 1.1 Two-dimensional illustration

![Maxwellian and Monte Carlo particle weights](monte_carlo_sampling.png)

The [plotting script](plot_monte_carlo_sampling.py) uses O₂⁺ at density 5 cm⁻³ (5×10⁶ m⁻³), bulk velocity (−10,0) km/s and physical temperature 10 eV. It draws 100,000 two-dimensional velocities with vz=0. The **sampling Maxwellian** has temperature 40 eV. A source area of **1 m²** and bulk speed **10 km/s** define the fixed total patch injection rate Q_cell = n A norm(U) = 5e10 s^-1. Both signs of each sampled velocity component are retained.

All three panels use turbo. The left panel shows the analytical two-dimensional Maxwellian (s² m⁻⁵). The middle and right panels color individual samples by `density_weight` (m⁻³) and `flux_weight` (s⁻¹), respectively. Both panels use the same normalized importance weights: Wn_i = n w_i / sum(w) and Q_i = Q_cell w_i / sum(w) = (Q_cell/n) Wn_i. Their normalized colors and spatial patterns are identical; only colorbar values and units differ. No individual radial-speed factor or sign cutoff is applied in this illustration. This figure demonstrates the prescribed injection model (`flux_model=:bulk_speed`), while the shell default described below remains the bidirectional crossing model. Reconstruct the VDF with weighted velocity-bin sums, not scatter-point colors.

### 1.2 Maxwellian sampling and importance weights

Let T be the physical temperature expressed as kBT in eV, κ_E=1.602176634×10⁻¹⁹ J/eV the energy conversion factor, ε_T the thermal energy in J, and m=5.352390155808×10⁻²⁶ kg the O₂⁺ mass. The one-component velocity standard deviation σ is in m/s:

$$
\epsilon_T=\kappa_E T,\qquad \sigma=\sqrt{\epsilon_T/m}.
$$

Let d be the dimensionless velocity-space dimension (2 in the illustration and 3 for the source). Velocities v and U are in m/s, the normalized probability density g_d has units (s/m)^d, number density n is in m⁻³, and f_d has units m⁻³(s/m)^d:

$$
g_d(\mathbf v)=\frac{\exp[-|\mathbf v-\mathbf U|^2/(2\sigma^2)]}{(2\pi\sigma^2)^{d/2}},\qquad f_d=n g_d.
$$

The sampling temperature T_s is in eV and its standard deviation σ_s is in m/s. The dimensionless temperature ratio c=4 specifies how much broader the sampling Maxwellian is:

$$
T_s=cT,\qquad \sigma_s^2=c\sigma^2.
$$

For samples v_i drawn from the normalized sampling density g_s, g_d and g_s have the same units and their ratio w_i is dimensionless. With v_i, U and σ in m/s and c and d dimensionless:

$$
w_i=\frac{g_d(\mathbf v_i)}{g_s(\mathbf v_i)}
=c^{d/2}\exp\left[-\frac{|\mathbf v_i-\mathbf U|^2}{2\sigma^2}\left(1-\frac1c\right)\right].
$$

The prefactor is c in the two-dimensional illustration and c^(3/2) in the three-dimensional implementation. See [sample_maxwellian_source and maxwellian_importance_weight_3d](../../../src/tracing/monte_carlo_weight.jl).

The density-ratio principle also applies to a κ distribution: replace the physical and sampling probability densities by the intended normalized distributions, and evaluate their ratio at each sample. The sampling distribution must cover the support of the physical distribution. Efficient sampling of κ tails generally requires a suitable heavy-tailed sampling distribution. The current source function implements Maxwellian sampling; κ sampling requires an additional sampler and density evaluation. The Maxwellian exponential expression above does not apply unchanged to κ distributions.

### 1.3 Density weights and bidirectional crossing rates

The diagnostic self-normalized density weights are

$$W_i^{(n)}=n w_i/\sum_jw_j.$$

They sum exactly to n, but finite-sample moment estimates are biased. They do not determine the default crossing-rate weights.

The source is a prescribed bidirectional Maxwellian crossing ensemble at 500 km. It launches particles both into the region above the source and into the region below it. Define vr = v dot er, with radial unit vector er. For N untruncated velocity draws from gs,

$$\boxed{Q_i=\frac{nA}{N}|\mathbf v_i\cdot\hat{\mathbf r}|\frac{g(\mathbf v_i)}{g_s(\mathbf v_i)}}.$$

Q_i is in s^-1 and is positive for both signs of vr when n>0 and vr is nonzero. It is a directional crossing-rate magnitude, not signed net flux. N includes every draw. There is no normalization by sum(w), by the number of draws of either sign, or by flight duration. `sample_maxwellian_source(...; normal=er, area_m2=A)` defaults to `flux_model=:bidirectional`.

For radial bulk drift Ur and physical one-component thermal standard deviation sigma, let varphi and Phi be the standard normal density and cumulative distribution:

$$F_+=n[\sigma\varphi(U_r/\sigma)+U_r\Phi(U_r/\sigma)],$$
$$F_-=n[\sigma\varphi(U_r/\sigma)-U_r\Phi(-U_r/\sigma)].$$

F+ and F- are positive outward and inward crossing magnitudes. The sampled total rate estimates A(F+ + F-); signed net flux is F+ - F- = n Ur. At zero drift each direction has F = n sigma / sqrt(2 pi). Tangential bulk drift does not change these local rates.

This source prescription treats initial inward and outward launches as independent injection histories. A propagated particle may cross 500 km again, but keeps its original Q_i; no new particle or source weight is created at a return crossing. MHD moments alone do not establish a volume-production rate, and this imposed crossing ensemble is not a self-consistent ion-production model.

The explicit `flux_model=:reservoir` uses max(vr,0), for an outward-only reservoir source. The compatibility option `flux_model=:bulk_speed` uses n A norm(U) w_i / sum(w), for a prescribed bulk-speed injection model. Shell selections are `reservoir_maxwellian_rate` and `n_bulk_speed_maxwellian`, respectively. Neither option is the default. The older conditional-outward prescribed-rate option also remains explicit.

### 1.4 Source at 500 km, absorption at 200 km

MHD density, Cartesian bulk velocity and temperature are evaluated at the area-coordinate midpoint of each native angular cell at 500 km altitude. The absorbing sphere is independently located at 200 km altitude. Both altitudes are above a Mars radius Rm=3390 km. With source radius rs in m, colatitude theta and longitude phi in rad,

$$A_{cell}=r_s^2(\cos\theta_0-\cos\theta_1)(\phi_1-\phi_0).$$

The default shell model is `bidirectional_maxwellian_rate`. Each cell supplies 100 untruncated three-dimensional velocity draws at Ts=4T. Inward draws carry positive normal-speed rate weights and propagate below 500 km. They may turn and return above the source; crossing the source is not a termination condition. Absorption occurs only on reaching the 200 km sphere. No volume production is included outside the source.

`Config.altitude_km` remains the source-altitude setting for API compatibility; `Config.absorption_altitude_km` is the separate lower boundary. Defaults are 500 and 200. The stored field mesh must cover both radii. The code validates coverage and does not extrapolate below the field domain.

All draws within a cell launch from the same midpoint. This is spatial midpoint quadrature, not uniform random area sampling. Spatial refinement must be checked, especially for localized detector contributions. `cell_stride>1` skips cells without compensating their area and therefore describes only the sampled patches, not a full-shell rate estimate.

At c=Ts/T=4, the asymptotic importance-weight ESS fraction is (2c-1)^(3/2)/c^3 = 0.289. Thus 100 draws give roughly 29 effective density samples per cell. This ESS does not measure uncertainty of crossing rates or detector occupancy. Compare independent seeds, particles per cell, sampling temperature and spatial resolution before interpreting rare trajectories.

### 1.5 Bulk-speed flux maps

![O2+ and O+ bulk-speed flux at 200, 400 and 600 km](bulk_speed_flux_200_400_600km.png)

See [map definitions, source data and reproduction](bulk_speed_flux_maps.md). All six panels share one turbo logarithmic scale. O+ is shown for comparison; this forward-tracing example still propagates O2+.

## 2. Forward tracing

### 2.1 Equations and functions

Position x is in m, velocity v in m/s, time t in s, mass m in kg, charge q in C, electric field E in V/m and magnetic field B in T. For O₂⁺, q=+1.602176634×10⁻¹⁹ C. The nonrelativistic equations are:

$$
\frac{d\mathbf x}{dt}=\mathbf v,\qquad
\frac{d\mathbf v}{dt}=\frac qm(\mathbf E+\mathbf v\times\mathbf B).
$$

[ShellMonteCarlo.run_monte_carlo](monte_carlo_shell.jl) calls `ShellMonteCarlo.trace_particle`, which uses `MarsTP.TP.boris_velocity_update`. The fields are static total MHD E and B. Integration uses half-step velocities, while saved velocities are synchronized with saved positions.

| Parameter | Setting |
| --- | --- |
| Mars radius | Rm=3390 km |
| Source altitude | 500 km |
| Absorbing inner-boundary altitude | 200 km |
| Outer radius | 4 Rm |
| Time step and maximum flight age | 0.1 s and 500 s |
| Omitted processes | Volume production, collisions, recombination, gravity and feedback |
| Termination | inner, outer, time_limit, zero_rate; numerical failures raise errors |

Each trajectory retains its Q_i during lossless propagation.

### 2.2 Trajectory illustration

![20,000 randomly selected O2+ trajectories in XY, XZ and YZ](trajectories_5000.png)

This historical figure has not been recomputed with the bidirectional default. It shows 20,000 randomly selected O₂⁺ forward trajectories using the updated `n*norm(U_bulk)` source. Panels show XY, XZ and YZ projections, from left to right. Dashed and dotted curves mark the bow shock (BS) and magnetic pileup boundary (MPB); in the YZ panel these boundaries are cross-sections at X=0. Positions are normalized by the Mars radius, Rm=3390 km. The image retains the filename `trajectories_5000.png` for link compatibility; the displayed sample contains 20,000 trajectories.

### 2.3 Satellite orbit and particle distributions

![O2+ trajectories, a circular orbit at 2.5 Rm and corresponding energy and reduced velocity distributions](satellite_orbit_distributions.png)

This historical figure has not been recomputed with the bidirectional default. The left panel combines 20,000 O₂⁺ trajectories with a prescribed circular satellite orbit in the Y=0 plane at a Mars-centered radius of **2.5 Rm**. The right panels show the corresponding direction-averaged DEF and one-dimensional reduced distributions in vx, vy and vz versus orbit angle. See [orbit geometry, distribution definitions and interpretation](satellite_orbit_distributions.md).

## 3. Detector PSD and integrated VDF

### 3.1 Particle rate times residence time

Let dN be particle count, x position in m, v velocity in m/s and n number density in m⁻³. The three-dimensional velocity distribution f has units s³ m⁻⁶:

$$
dN=f(\mathbf x,\mathbf v)d^3x\,d^3v,\qquad n(\mathbf x)=\int f(\mathbf x,\mathbf v)d^3v.
$$

The cube side L=0.2 Rm=678,000 m and volume V_D is in m³. Velocity-bin widths Δv_x, Δv_y and Δv_z are in m/s, giving velocity-space volume Δ³v_b in m³ s⁻³:

$$
V_D=L^3,\qquad \Delta^3v_b=\Delta v_x\Delta v_y\Delta v_z.
$$

Let τ_i,D,b be residence time in s, T_i the saved flight duration in s, x_i position in m and v_i velocity in m/s. Indicator functions for the spatial cube D and velocity bin B_b are dimensionless:

$$
\tau_{i,D,b}=\int_0^{T_i}\mathbf1_D[\mathbf x_i(t)]\mathbf1_{B_b}[\mathbf v_i(t)]\,dt.
$$

For a steady source, Q_i is in s⁻¹, τ_i,D,b in s, and N_D,b is mean particle count. With V_D in m³ and Δ³v_b in m³ s⁻³, the bin-averaged PSD f̄_D,b is in s³ m⁻⁶:

$$
N_{D,b}=\sum_iQ_i\tau_{i,D,b},\qquad
\boxed{\bar f_{D,b}=\frac{\sum_iQ_i\tau_{i,D,b}}{V_D\Delta^3v_b}}.
$$

Q_i already contains sampling normalization, so there is no additional division by particle number or total integration time. For steady sources and static fields, the maximum flight age sets the truncation of the residence integral. This detector estimator does not depend on choosing a Maxwellian source.

### 3.2 Functions and binning

- [forward_psd](../../../src/tracing/detector_psd_forward.jl): trajectories in memory.
- [forward_psd_saved and foreach_saved_trajectory](../../../src/tracing/trajectory_io.jl): trajectories on disk.
- [analyze_saved_probes.jl](analyze_saved_probes.jl): multiple detectors in one pass through the files.

Position and synchronized velocity are linearly interpolated between saved states. Let α, α_a and α_b be dimensionless segment fractions; x_0 and x_1 are in m, v_0 and v_1 in m/s, and Δt and δt in s:

$$
\mathbf x(\alpha)=\mathbf x_0+\alpha(\mathbf x_1-\mathbf x_0),\quad
\mathbf v(\alpha)=\mathbf v_0+\alpha(\mathbf v_1-\mathbf v_0),\quad
\delta t=\Delta t(\alpha_b-\alpha_a).
$$

Each segment is clipped to the cube, then split at velocity-bin boundaries. Each piece adds its rate times duration. Crossings are detected even when both saved endpoints lie outside the cube. Spatial intervals are lower-inclusive and upper-exclusive; the final global velocity edge belongs to the last bin.

Each velocity axis spans −500 to 500 km/s with 5 km/s bins (200³ bins). Nonzero three-dimensional PSD values are stored sparsely.

```julia
using MarsTP
Rm = 3_390_000.0 # m
settings = (; detector_m=[0.,0.,2.] * Rm, side_m=0.2Rm,
             vlim=(-500.,500.), vgrid=200, velocity_unit=:km_s,
             species="O2+", coordinate_system="native MHD Cartesian axes")

# Synchronized in-memory times and states: s, m, m/s. Q is in s^-1.
r3 = forward_psd(solutions; settings..., rate_weights_s=Q, option="3D")
# Read saved trajectories individually.
r3 = forward_psd_saved("outputs/my_run"; settings..., option="3D", storage=:sparse)
rxy = forward_psd_saved("outputs/my_run"; settings..., option="Vx-Vy")
```

`vgrid=200` is the number of bins per axis. `velocity_unit=:km_s` specifies the input range unit; output velocity edges and centers remain in m/s. `storage=:sparse` returns a sparse three-dimensional PSD.

### 3.3 Velocity integration

Let f_ijk be the three-dimensional PSD in s³ m⁻⁶ and Δv_y=Δv_z=5000 m/s. The integrated distributions F_xy and F_xz have units s² m⁻⁵:

$$
F_{xy}(v_{x,i},v_{y,j})=\sum_kf_{ijk}\Delta v_{z,k},\qquad
F_{xz}(v_{x,i},v_{z,k})=\sum_jf_{ijk}\Delta v_{y,j}.
$$

The omitted axis is integrated over its full ±500 km/s range. Plot axes show ±300 km/s, while color values retain SI units.

Let n_D be detector density in m⁻³, f in s³ m⁻⁶, F in s² m⁻⁵, Δv in m/s, Q in s⁻¹, τ in s and V_D in m³:

$$
n_D=\sum_{ijk}f_{ijk}\Delta v^3
=\sum_{ij}F_{xy,ij}\Delta v^2
=\sum_{ik}F_{xz,ik}\Delta v^2
=\frac{\sum_iQ_i\tau_{i,D}}{V_D}.
$$

The last equality requires the velocity range to contain all residence contributions.

### 3.4 Detector illustrations

Center (1,0,2) Rm:

![Probe at 1 0 2](probe_psd_projections.png)

Colors show velocity-integrated VDFs in s² m⁻⁵ with a turbo logarithmic scale.

### 3.5 Surface flux and residence density

Let A_face be face area in m², Q_i rate in s⁻¹, and C the set of crossings for a chosen face and direction. The surface flux Γ is in m⁻² s⁻¹:

$$
\Gamma=\frac{\sum_{i\in C}Q_i}{A_{\rm face}}.
$$

For a uniform beam crossing a cube normally, let n be density in m⁻³, u speed in m/s, L side length in m, Q rate in s⁻¹ and τ residence time in s:

$$
Q=nuL^2,\qquad \tau=L/u,\qquad \frac{Q\tau}{L^3}=n.
$$

This monoenergetic beam identity checks the detector estimator. For the prescribed broad-distribution source, individual rates scale with bulk speed and importance weight, so the source is not a thermal reservoir crossing model.

> The trajectory, satellite and detector PSD figures are historical results and have not been recomputed with the bidirectional default. Only the sampling illustration reflects the corrected source weights.

### 3.6 Omnidirectional differential energy flux

`detector_omni_def` computes the 4π direction-averaged DEF in eV/(m² s eV sr), without an isotropy assumption or angular binning. Both forward detector APIs accept `energy_edges_eV` and return `omni_def`; backward runs use `BacktraceConfig(energy_edges_eV=...)`. See [definitions, interfaces and examples](detector_omni_def.md).

## 4. Running and saving

Run from the repository root using the Julia project environment. Sampling and tracing require `data/mars_fields_spherical_from_dat.vts`.

```julia
include("examples/forward_tracing/monte_carlo_forward_tracing/monte_carlo_shell.jl")
ShellMonteCarlo.run_monte_carlo("outputs/my_run",
    ShellMonteCarlo.Config(per_cell=100, dt=0.1, tmax=500., batch_size=1024,
                          altitude_km=500., absorption_altitude_km=200.,
                          flux_model="bidirectional_maxwellian_rate", compress_trajectories=false))
```

[write_trajectory_batch](../../../src/tracing/trajectory_io.jl) saves batches of paths and weights. Position, velocity and time use m, m/s and s. `rate_weights_s` is in s⁻¹ and `source_density_weights_m3` in m⁻³:

```julia
using MarsTP
write_trajectory_batch("outputs/my_run/trajectories_00001.jld2", trajectory_batch;
    particle_ids=ids, rate_weights_s=Q,
    source_density_weights_m3=density_weights, cell_ids=cell_ids,
    termination_codes=statuses, species="O2+",
    coordinate_system="native MHD Cartesian axes", compress=false)
```

State order is t,x,y,z,vx,vy,vz; particle IDs, weights and termination types are also saved. `compress=true` enables lossless compression. The output file must not already exist.

Compute and plot detector PSD:

```powershell
$example = 'examples/forward_tracing/monte_carlo_forward_tracing'
julia --project=. "$example/analyze_saved_probes.jl" outputs/my_run outputs/my_probes
python "$example/plot_library_psd.py" outputs/my_probes
```

Regenerate the sampling illustration:

```powershell
python examples/forward_tracing/monte_carlo_forward_tracing/plot_monte_carlo_sampling.py
```

Python requires NumPy, Matplotlib and h5py. Trajectory backgrounds use [src/visualization](../../../src/visualization/README.md).

| Code | Purpose |
| --- | --- |
| [plot_monte_carlo_sampling.py](plot_monte_carlo_sampling.py) | Maxwellian and particle-weight illustration |
| [monte_carlo_shell.jl](monte_carlo_shell.jl) | Source sampling, Boris tracing and batch output |
| [plot_trajectories.py](plot_trajectories.py) | Trajectory projections |
| [analyze_saved_probes.jl](analyze_saved_probes.jl) | Detector PSD from saved paths |
| [plot_library_psd.py](plot_library_psd.py) | Integrated two-dimensional VDFs |

## Historical validation of the bulk-speed update

The following results refer to the earlier bulk-speed implementation, not the corrected bidirectional default. See the reservoir validation below.

Validated with Julia 1.12.6 and TestParticle 0.23.3. The package suite passed 420 assertions and the shell suite passed 240 assertions. Checks include exact patch-rate normalization, both radial velocity signs with positive weights, independence from normal direction, zero rate at zero bulk speed, the explicit legacy estimator, and immediate absorption of inward launches.

A local-field smoke run used 6 source cells, 10 particles per cell, dt=0.1 s and maximum age 0.2 s. All 60 saved particles had positive rate weights; 28 inward launches terminated at time zero and 32 reached the time limit. Total prescribed injection was 2.984432090205178e21 s^-1 for the sampled cells, with per-cell sums and all saved rate weights verified. This short run is an implementation check, not an escape-rate or detector-PSD convergence result. The full million-particle ensemble was not rerun.

```powershell
julia --startup-file=no --compiled-modules=existing --project=. test/runtests.jl
julia --startup-file=no --compiled-modules=existing --project=. examples/forward_tracing/monte_carlo_forward_tracing/test_monte_carlo.jl
julia --startup-file=no --compiled-modules=existing --project=. examples/forward_tracing/monte_carlo_forward_tracing/monte_carlo_shell.jl outputs/new_bulk_speed_smoke 2000 0.2 0.1 10 n_bulk_speed_maxwellian
```

## Validation of the separated source and absorption boundaries

With Julia 1.12.6 and TestParticle 0.23.3, the corrected package suite passed 532 assertions and the shell suite passed 257 assertions. Tests cover bidirectional crossing-rate sums, positive inward weights, outward-only compatibility, zero-field slab density recovery, travel from 500 km to absorption at 200 km, and an inward launch that turns under an outward electric field and returns above 500 km without termination at the source.

A local MHD smoke check used 6 sampled cells, 10 draws per cell and 0.2 s of propagation. The actual stored mesh minimum is 199.646643 km. All 28 inward launches out of 60 samples carried positive rate and propagated to the time limit, with none terminated at the 500 km source. This short check is not a convergence result. The updated prescribed-rate 2D illustration has total density 5e6 m^-3 and total injection rate 5e10 s^-1. Its per-particle rate weights are exactly the density weights multiplied by A norm(U); this plot is a separate source prescription from the bidirectional shell default.

These controlled tests establish implementation behavior, not the frequency of returning trajectories in the actual MHD field. Historical trajectory and detector figures have not been recomputed with the separated boundaries.
