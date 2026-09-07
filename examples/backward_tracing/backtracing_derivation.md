# O₂⁺ backtracing: volume and thin-shell sources

This derivation describes `thin_shell_source_v1` in `src/tracing/detector_psd_backward.jl`. The ionosphere is a traversable source surface, normally at 400 km altitude. The absorbing inner boundary remains at 200 km.

## 1. Model and coordinates

Nonrelativistic O₂⁺ trajectories evolve in static MHD electromagnetic fields. Mass and positive charge come from TestParticle `SpeciesDict["O2+"]`. Positions and velocities use Cartesian components; the MHD spherical grid uses radius, colatitude and azimuth.

Let R_M=3.390×10⁶ m, r=|x| in m and altitude h=r−R_M in m. The inner and outer radii and source radius are in m:

$$
r_{\rm in}=R_M+200\times10^3\ \mathrm{m},\qquad r_{\rm out}=4R_M,
\qquad r_s=R_M+h_s,\qquad h_s=400\times10^3\ \mathrm{m}.
$$

`ionosphere_altitude_km` accepts 200 to 800 km without changing the inner boundary. Volume production exists throughout the integration domain, including below the source shell. Surface production exists only at r=r_s. A particle observed at 300 km can nevertheless carry a contribution from an earlier 400 km crossing.

The model omits collisions, losses, gravity, chemical feedback, self-consistent particle fields and time-varying sources. Outflow rates at different heights are not remapped to 400 km.

## 2. Definitions and units

| Quantity | Definition | SI unit |
| --- | --- | --- |
| x, r, R_M | Position or Mars-centered radius | m |
| t, Δt | Physical time and step | s |
| v | Instantaneous particle velocity | m s⁻¹ |
| U_i | MHD O₂⁺ bulk velocity | m s⁻¹ |
| E, B | Electric and magnetic fields | V m⁻¹, T |
| n_i | MHD O₂⁺ number density | m⁻³ |
| T_i, T_n | Ion and neutral temperatures | K |
| m_i, q_i | Particle mass and charge | kg, C |
| Q_V | Volume production rate | m⁻³ s⁻¹ |
| F_s | Surface production rate | m⁻² s⁻¹ |
| g_V, g_M | Normalized three-dimensional velocity densities | s³ m⁻³ |
| δ(r−r_s) | Radial Dirac delta | m⁻¹ |
| S_V, S_s | Phase-space production rates | s² m⁻⁶ |
| f, f_V, f_s | Velocity-space PSD | s³ m⁻⁶ |
| f_xz | Distribution integrated over v_y | s² m⁻⁵ |

Particle number dN is a dimensionless count. Using x in m, v in m/s, n in m⁻³ and f in s³ m⁻⁶:

$$
dN=f(\mathbf x,\mathbf v,t)d^3x\,d^3v,\qquad
n(\mathbf x,t)=\int f(\mathbf x,\mathbf v,t)d^3v.
$$

For nonrelativistic momentum p=m_i v, the velocity and momentum distributions obey f_v=m_i³f_p; their numerical values cannot be compared directly. `velocity_axes` converts km/s to m/s before integration. Display fluxes in cm⁻² s⁻¹ are obtained by dividing SI fluxes by 10⁴.

## 3. Transport along characteristics

Use the quantities and SI units above, with acceleration a in m/s². The source-driven, lossless phase-space continuity equation is:

$$
\frac{\partial f}{\partial t}
+\nabla_{\mathbf x}\cdot(\mathbf v f)
+\nabla_{\mathbf v}\cdot(\mathbf a f)=S_V+S_s,
\qquad \mathbf a=\frac{q_i}{m_i}(\mathbf E+\mathbf v\times\mathbf B).
$$

Lorentz flow has zero phase-space divergence. With x in m, v in m/s, a in m/s² and source terms in s² m⁻⁶, the characteristic equations are:

$$
\dot{\mathbf x}=\mathbf v,\qquad \dot{\mathbf v}=\mathbf a,\qquad
\frac{df}{dt}=S_V[\mathbf x(t),\mathbf v(t)]+S_s[\mathbf x(t),\mathbf v(t)].
$$

Let observation time t_0=0 and past endpoint t_*<0 be in s, with f in s³ m⁻⁶ and S_V,S_s in s² m⁻⁶:

$$
f(\mathbf x_0,\mathbf v_0,0)=f[\mathbf x(t_*),\mathbf v(t_*),t_*]
+\int_{t_*}^{0}(S_V+S_s)\,dt.
$$

The implementation sets the imposed past-endpoint background to zero. It terminates at 200 km, 4R_M or the lookback limit. A time-limited answer is a truncated source integral, not automatically a steady-state solution. Nonzero external-boundary background VDFs are not provided.

Negative time steps implement backtracing without reversing mass, charge or physical velocity. Defining positive lookback age τ=−t expresses the source integral over increasing τ. This is why positive production accumulates with |Δt|.

## 4. Volume production

Q_V(x) from MAT inputs has units m⁻³ s⁻¹. Let T_n be GITM neutral temperature in K, k_B the Boltzmann constant in J/K, m_i mass in kg and w_n speed in m/s. The normalized birth density g_V has units s³ m⁻³:

$$
g_V(\mathbf v;\mathbf x)=\frac{1}{\pi^{3/2}w_n^3}
\exp\left(-\frac{|\mathbf v|^2}{w_n^2}\right),\qquad
w_n=\sqrt{\frac{2k_BT_n(\mathbf x)}{m_i}},\qquad \int g_Vd^3v=1.
$$

The volume source uses zero drift and neutral temperature, whereas the shell uses MHD ion temperature and bulk velocity.

With Q_V in m⁻³ s⁻¹, g_V in s³ m⁻³, S_V in s² m⁻⁶, time in s and f_V in s³ m⁻⁶:

$$
S_V=Q_Vg_V,\qquad
\boxed{f_V=\int_{t_*}^{0}Q_V[\mathbf x(t)]g_V[\mathbf v(t);\mathbf x(t)]dt}.
$$

The distribution is evaluated at the instantaneous trajectory velocity. The integration variable is time; no additional speed or grid-volume factor is needed.

## 5. MHD surface production and birth velocities

The shell samples `n_O^2^p [m^-3]`, `T_O^2^p [K]` and `U_O^2^p [m/s]`. Let Ω represent angular location, n_i density in m⁻³, U_i velocity in m/s and F_s surface production in m⁻² s⁻¹:

$$
\boxed{F_s(\Omega)=n_i(r_s,\Omega)|\mathbf U_i(r_s,\Omega)|}.
$$

This is a source prescription, distinct from signed normal flux or an integrated one-sided thermal flux. It converts local moments into an assumed independent surface source. It does not imply that MHD directly provides an independent birth rate.

Let T_i be in K, m_i in kg, k_B in J/K, and v,U_i,w_i in m/s. The normalized drifting Maxwellian g_M has units s³ m⁻³:

$$
g_M(\mathbf v;\Omega)=\frac{1}{\pi^{3/2}w_i^3}
\exp\left[-\frac{|\mathbf v-\mathbf U_i|^2}{w_i^2}\right],\qquad
w_i=\sqrt{\frac{2k_BT_i}{m_i}},\qquad \int g_Md^3v=1.
$$

The density factor is already in F_s. `ionosphere_distribution` returns both `g=g_M` and `f=n_i*g_M`; the shell uses `flux*g` to avoid double-counting density.

The full three-dimensional Maxwellian includes both inward and outward velocities. Restricting emission to outward directions would require a redefined normalized distribution. This thin-shell prescription differs from the reservoir-crossing source used in the forward Monte Carlo example.

## 6. Radial velocity in the shell contribution

Let F_s be in m⁻² s⁻¹, g_M in s³ m⁻³ and the radial delta in m⁻¹. The shell phase-space source S_s is in s² m⁻⁶:

$$
S_s(\mathbf x,\mathbf v)=F_s(\Omega)g_M(\mathbf v;\Omega)\delta(r-r_s).
$$

With d³x in m³, area dA and r_s² in m², solid angle dΩ dimensionless and F_s in m⁻² s⁻¹, integrating over space gives a rate in s⁻¹:

$$
\int F_s\delta(r-r_s)d^3x=\int F_s r_s^2d\Omega=\int F_s dA.
$$

Area is already represented by the spatial measure and delta function; no extra r_s² belongs in S_s.

Let t_j be a transverse crossing time in s, r in m, v in m/s and r̂_j a dimensionless radial unit vector. With δ(t−t_j) in s⁻¹ and dr/dt in m/s:

$$
\delta[r(t)-r_s]=\sum_j\frac{\delta(t-t_j)}{|\dot r(t_j)|},\qquad
\dot r(t_j)=\mathbf v(t_j)\cdot\hat{\mathbf r}_j.
$$

Using F_s in m⁻² s⁻¹, g_M in s³ m⁻³ and radial speed in m/s, the shell PSD f_s is in s³ m⁻⁶:

$$
\boxed{f_s=\sum_j\frac{F_s(\Omega_j)g_M[\mathbf v(t_j);\Omega_j]}
{|\mathbf v(t_j)\cdot\hat{\mathbf r}_j|}}.
$$

The denominator is the instantaneous test-particle radial speed, not the MHD radial bulk speed. It describes crossing residence geometry and does not change the definition F_s=n_i|U_i|. Repeated transverse crossings add independently, shared step endpoints count once, and the source shell neither changes particle velocity nor terminates the trajectory.

## 7. Total PSD and integrated output

With zero imposed background, let f,f_V,f_s have units s³ m⁻⁶. Their source terms and velocities retain the units above:

$$
\boxed{f(\mathbf x_0,\mathbf v_0)=
\int_{t_*}^{0}Q_V[\mathbf x(t)]g_V[\mathbf v(t);\mathbf x(t)]dt
+\sum_j\frac{n_i(\Omega_j)|\mathbf U_i(\Omega_j)|g_M[\mathbf v(t_j);\Omega_j]}
{|\mathbf v(t_j)\cdot\hat{\mathbf r}_j|}}.
$$

Let f be in s³ m⁻⁶ and v_y,Δv_y in m/s. The integrated VDF f_xz is in s² m⁻⁵:

$$
f_{xz}(v_x,v_z)=\int f(v_x,v_y,v_z)dv_y
\simeq\sum_\ell f(v_x,v_{y,\ell},v_z)\Delta v_y.
$$

| Output | Meaning | Unit |
| --- | --- | --- |
| `f2d_volume` | Integrated volume contribution | s² m⁻⁵ |
| `f2d_ionosphere` | Integrated shell contribution | s² m⁻⁵ |
| `f2d_xz` | Sum of both terms | s² m⁻⁵ |
| `ionosphere_crossings` | Sum of crossing counts over sampled Vy at each (vx,vz) | Count, without Δv_y |
| `vx_km,vy_km,vz_km` | Detector velocity grid | km/s |
| `model` | `thin_shell_source_v1` | Identifier |
| `source_units` | Source and PSD unit dictionary | Strings |

Even a single Vy point is multiplied by `dvy_kms*1000`, representing a rectangular integral over that width rather than a three-dimensional slice. Density requires integration over vx and vz and checks of velocity coverage and resolution. The output is not differential flux per eV.

## 8. Discrete implementation

1. Trace the detector state backward in SI units with Boris or AdaptiveBoris.
2. Clip each step at the 200 km or 4R_M boundary.
3. Find all intersections of the valid position segment with the source sphere, in time order.
4. Synchronize the Boris half-step velocity at each intersection and evaluate the shell contribution.
5. Integrate the volume source over pieces separated by crossings and the valid endpoint.

For S_V in s² m⁻⁶ and t in s, the trapezoidal increment Δf_V is in s³ m⁻⁶:

$$
\Delta f_V\simeq\frac{S_{V,k}+S_{V,k+1}}2|t_{k+1}-t_k|.
$$

6. Continue to a boundary or time limit, accumulating volume production below the shell as well.
7. Accumulate the two three-dimensional contributions separately, then integrate over Vy.

Intersections use straight position segments and synchronized velocities from the same Boris step. Curved paths and multiple unresolved crossings introduce finite-step error. Fields are never evaluated outside the domain; 1e-8 m offsets near boundaries only avoid roundoff errors.

Status indices are 1=unused, 2=time limit, 3=200 km boundary, 4=nonfinite path, 5=outer boundary. A shell crossing is not termination. Partial accumulations from nonfinite paths may remain in output; check `status_counts[4]==0` before scientific use. Invalid source values or moments raise errors.

## 9. Zero-thickness limitations

Shell PSD can grow as radial speed approaches zero. Exact tangency violates the simple-root delta transformation. The implementation does not impose an arbitrary minimum physical radial speed. It rejects unresolved tangencies or radial speeds below `sqrt(eps(Float64))*max(norm(v),1 m/s)`.

A detector exactly on the source shell is rejected because of ambiguous initial-time attribution. At a source altitude of 200 km, the code uses the one-sided limit from inside the domain and counts one full crossing at the absorbing boundary.

A finite-thickness treatment would require a specified thickness and normalized radial profile W(r), in m⁻¹, satisfying an integral of one over radius. That additional physical model is not implemented.

Local surface production at 300 km is zero, while a backtraced path from there can cross 400 km and acquire a source contribution. This is accumulation along characteristics, not spatial spreading of the shell.

Some MHD pole moments are nearly zero but finite. Plot masks do not modify transport inputs. The model adds MAT volume production and the shell source without automatically removing physical overlap; this assumption needs independent assessment for total-production studies.

## 10. Run

```julia
using MarsTP, StaticArrays
cfg = BacktraceConfig(
    detector_Rm=SA[0.0,0.0,2.0],
    ionosphere_altitude_km=400.0,
    include_ionosphere=true,
    dt=-0.05,
    tspan=(0.0,-500.0),
    vx_min_kms=-20.0, vx_max_kms=20.0,
    vy_min_kms=-20.0, vy_max_kms=20.0,
    vz_min_kms=-20.0, vz_max_kms=20.0,
    dv_kms=10.0, dvy_kms=10.0,
)
result = run_backtrace_vdf(cfg)
```

This small velocity range illustrates the interface, not a converged scientific VDF. Prepare the Julia project, MHD, source-rate and GITM inputs. Use distinct output paths for different configurations.

```sh
julia --project=. test/runtests.jl
julia --project=. scripts/smoke_ionosphere.jl
```

Tests cover density/PDF distinctions, flux units, 200/400/800 km shells, the fixed inner boundary, sub-shell volume production, outward backtraced crossings, two-crossing geometry, shared endpoints, tangencies, synchronized velocities and step halving. They do not replace velocity-grid and lookback-time convergence checks.

The original Julia 1.12.6/TestParticle 0.23.3 validation passed 122 tests. A two-second MHD trajectory starting at (1.12 Rm,0,0) with velocity (10,0,0) km/s and dt=−0.05 s crossed the 400 km shell and continued to the time limit. Its shell contribution integrated over a single 1 km/s Vy width was about 1.17466e-5 s² m⁻⁵, an execution check rather than a converged prediction.
