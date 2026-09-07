# Maxwellian Monte Carlo source weights

Ported from [MarsASPEN.jl](https://github.com/ChiZhang123456/MarsASPEN.jl/blob/main/src/monte_carlo_weight.jl), blob `41d84d67602c5fa284f172a540cccc46531ca8e0`.
The original file supplies weights, not random draws or a transport solver.
MarsTP adds `sample_maxwellian_source`; default particle mass is O2+, not H.
Existing untracked shell Monte Carlo scripts are independent and unchanged.

For physical temperature T, sampling Maxwellian temperature Ts=c*T and bulk velocity U,
each Cartesian velocity component is sampled as U_k+sqrt(e*Ts/m)*randn().
Temperatures here mean kT/e in eV; convert input Kelvin using kB*T/e.
All three components are in the same Cartesian frame as the fields.
The returned initial state is [x,y,z,vx,vy,vz], in m and m/s.

The dimensionless importance ratio is

    w = g(v;U,T)/g(v;U,Ts)

where g is a normalized Maxwellian, in s^3 m^-3. With c=1 every w=1.
The sampled states plus weights constitute a discrete representation of the
source population, not a gridded VDF or a density measured at a detector.

The port preserves the original self-normalized density weights:

    Wn_i = n*w_i/sum(w)                      [m^-3]

Their sum is exactly the specified local density n. Self-normalized moments
have finite-sample bias. Zero source density selects the original unitless
mode, so `density_weights_m3` is then `nothing`; it does not denote a physical
zero-density ensemble. `macro_weights` is Wn, or unit_particle_weight*w in
unitless mode. A density weight is not a number of particles or a rate.

For an explicitly specified patch area A and outward unit normal er:

    rate_weights_s[i] = n*A*max(dot(v_i,er),0)*w_i/N  [s^-1]

N is the total number of Maxwellian draws, including inward ones. Inward draws
have zero rate; filter using their indices and retain the corresponding
weights. Do not renormalize by the number of retained particles. Sum these
weights to estimate the outward source rate of this patch. This estimator
represents a reservoir Maxwellian crossing a surface. At U=0 it has the
analytic total n*A*sqrt(e*T/m)/sqrt(2*pi), which is not zero.

This differs from the existing backtrace thin-sheet source F*g with
F=n*norm(U). Do not compare the two as identical physical sources. It also
differs from the existing shell example's outward-conditioned Maxwellian
with a separately prescribed injection flux. No such model is silently
replaced by this port. See `maxwellian_source.jl` for the explicit local MHD
moment adapter and forward-tracing connection. A shell simulation needs a
specified area discretization and independent reproducible RNG streams per
patch; one position is only a local patch approximation.

After propagation, sum rate weights for particles classified as escaped to
estimate escape rate [s^-1], counting each injected particle once. For a
steady source, cell density is sum(rate_weight*residence_time)/cell_volume
[m^-3]. Flux requires division by collecting area, and a velocity histogram
VDF additionally requires the velocity-bin volume. No detector or escape
quantity is produced by this sampler alone. Converge sample count, sampling
temperature, timestep and termination settings before interpreting them.
The sampler reports effective_sample_size = sum(w)^2/sum(w^2), dimensionless;
it diagnoses velocity importance sampling, not uncertainty of escape flux.

Validation: `julia --compiled-modules=existing --project=. test/runtests.jl`.
Tests cover SI thermal speed, fixed-seed reproducibility, normalization,
invalid inputs, weighted second moments, and analytic zero-drift outward
flux at sampling-temperature factors 1 and 4. No large trajectory simulation
is part of these tests. No new external dependency is needed; Random is a
Julia standard library.

## Other velocity distributions

Importance weighting uses the ratio of normalized physical and sampling velocity densities. The same principle and source-rate weighting can be applied to a κ distribution, provided the sampler covers its support and both densities are evaluated consistently. The temperature-rescaling exponential given above is specific to Maxwellians. A suitable heavy-tailed sampling distribution is generally preferable for κ tails. `sample_maxwellian_source` currently implements Maxwellian sampling only; a κ sampler requires additional implementation.
