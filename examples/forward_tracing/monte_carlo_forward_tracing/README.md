# O₂⁺ Monte Carlo sampling, forward tracing and detector PSD

This example samples O₂⁺ source particles, follows their trajectories through MHD fields, and reconstructs velocity distributions from residence times inside a cubic detector.

## 1. Monte Carlo source sampling

Monte Carlo represents a continuous velocity distribution with a finite set of weighted particles. Here, velocities are drawn from a broader Maxwellian to improve sampling of the high-speed tail. Importance weights recover the physical distribution.

### 1.1 Two-dimensional illustration

![Maxwellian and Monte Carlo particle weights](monte_carlo_sampling.png)

The [plotting script](plot_monte_carlo_sampling.py) uses O₂⁺ at density 5 cm⁻³ (5×10⁶ m⁻³), bulk velocity (−10,0) km/s and physical temperature 10 eV. It draws 100,000 two-dimensional velocities with vz=0. The **sampling Maxwellian** has temperature 40 eV. A source area of **1 m²** with normal **−X** defines the particle-rate weights.

All three panels use turbo. The left panel shows the analytical two-dimensional Maxwellian (s² m⁻⁵). The middle and right panels color individual samples by `density_weight` (m⁻³) and `flux_weight` (s⁻¹), respectively. Gray points have zero crossing rate.

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

### 1.3 Density and particle-rate weights

Let N be the dimensionless number of draws, n the number density in m⁻³, and w_i dimensionless importance weights. The density weight W_i^(n), in m⁻³, is:

$$
W_i^{(n)}=n\frac{w_i}{\sum_{j=1}^{N}w_j},\qquad \sum_iW_i^{(n)}=n.
$$

`particle_density_weight` evaluates this expression.

Let A be source area in m², n̂ a dimensionless unit normal, v_i and v_n,i velocities in m/s, n the number density in m⁻³, and N and w_i dimensionless. The particle-rate weight Q_i is in s⁻¹, treating particle number as a count:

$$
v_{n,i}=\mathbf v_i\cdot\hat{\mathbf n},\qquad
Q_i=\frac{nA}{N}\max(v_{n,i},0)w_i.
$$

Q_i corresponds to `rate_weights_s`, called `flux_weight` in the illustration. Inward samples have zero rate; N includes all draws. The same area-crossing expression applies to other normalized velocity distributions when the appropriate importance ratio is supplied.

### 1.4 The 500 km ionospheric source

MHD density, bulk velocity and temperature are read at each source cell at 500 km altitude. The source surface is also the absorbing inner boundary, with a radial outward normal. Let r_s be the source radius in m, θ colatitude in rad and ϕ longitude in rad. Cell area A_cell is in m²:

$$
A_{\rm cell}=r_s^2(\cos\theta_0-\cos\theta_1)(\phi_1-\phi_0).
$$

Each angular cell supplies 100 three-dimensional velocity samples. Weights use the local density, cell area and radial particle velocity. The sampling Maxwellian temperature is four times the local physical temperature. No volume production is included outside the ionosphere.

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
| Source and inner-boundary altitude | 500 km |
| Outer radius | 4 Rm |
| Time step and maximum flight age | 0.1 s and 500 s |
| Omitted processes | Volume production, collisions, recombination, gravity and feedback |
| Termination | inner, outer, time_limit, zero_rate; numerical failures raise errors |

Each trajectory retains its Q_i during lossless propagation.

### 2.2 Trajectory illustration

![5000 trajectories in XZ, XY and YZ](trajectories_5000.png)

[plot_trajectories.py](plot_trajectories.py) shows 5000 paths projected onto XZ, XY and YZ.

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

Center (0,0,2) Rm:

![Probe at 0 0 2](probe_0_0_2/probe_psd_projections.png)

Center (−1.5,0,1) Rm:

![Probe at minus1.5 0 1](probe_m1p5_0_1/probe_psd_projections.png)

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

Faster particles have higher crossing rates but shorter residence times. Their product recovers density.

## 4. Running and saving

Run from the repository root using the Julia project environment. Sampling and tracing require `data/mars_fields_spherical_from_dat.vts`.

```julia
include("examples/forward_tracing/monte_carlo_forward_tracing/monte_carlo_shell.jl")
ShellMonteCarlo.run_monte_carlo("outputs/my_run",
    ShellMonteCarlo.Config(per_cell=100, dt=0.1, tmax=500., batch_size=1024,
                          flux_model="reservoir_maxwellian_rate", compress_trajectories=false))
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
