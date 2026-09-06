# Adapted from ChiZhang123456/MarsASPEN.jl src/monte_carlo_weight.jl,
# Git blob 41d84d67602c5fa284f172a540cccc46531ca8e0.
# MarsTP defaults to O2+ rather than hydrogen. SI except temperature in eV.
Base.@kwdef struct MonteCarloWeight
    sampling_temperature_factor::Float64 = 1.0
    source_number_density_m3::Float64 = 0.0
    unit_particle_weight::Float64 = 1.0
end

_mc_positive(x) = isfinite(x) && x > 0
_mc_nonnegative(x) = isfinite(x) && x >= 0
const _MC_QE = 1.602176634e-19

"""Maxwellian thermal speed sqrt(2 kT/m) [m/s]; temperature_ev means kT/e."""
function thermal_speed_from_temperature_ev(temperature_ev::Real,
        mass_kg::Real=TP.SpeciesDict["O2+"].m)
    _mc_positive(temperature_ev) || throw(ArgumentError("temperature must be finite and positive"))
    _mc_positive(mass_kg) || throw(ArgumentError("mass must be finite and positive"))
    sqrt(2 * Float64(temperature_ev) * _MC_QE / Float64(mass_kg))
end

"""Dimensionless normalized Maxwellian ratio g(v;U,T)/gs(v;U,Ts)."""
function maxwellian_importance_weight_3d(U, v, temperature_ev::Real,
        sampled_temperature_ev::Real; mass_kg::Real=TP.SpeciesDict["O2+"].m)
    length(U) == length(v) == 3 && all(isfinite,U) && all(isfinite,v) ||
        throw(ArgumentError("velocities must contain three finite components"))
    a = thermal_speed_from_temperature_ev(temperature_ev,mass_kg)
    b = thermal_speed_from_temperature_ev(sampled_temperature_ev,mass_kg)
    d2 = sum((Float64(v[k])-Float64(U[k]))^2 for k in 1:3)
    exp(3log(b/a) + d2*(inv(b^2)-inv(a^2)))
end

"""Self-normalized density contribution [m^-3]; sum across one source equals n."""
function particle_density_weight(w::Real,n::Real,total::Real)
    _mc_nonnegative(w) && _mc_nonnegative(n) && _mc_positive(total) ||
        throw(ArgumentError("finite nonnegative weight/density and positive total required"))
    Float64(n)*(Float64(w)/Float64(total))
end

"""
    sample_maxwellian_source(N; position_m, bulk_velocity_m_s, temperature_ev,
        species="O2+", weights=MonteCarloWeight(), rng=Random.default_rng(),
        normal=nothing, area_m2=nothing)

Sample N velocities from an untruncated drifting Maxwellian at Ts=factor*T.
Return initial_states ([x,y,z,vx,vy,vz], m and m/s), dimensionless
importance_weights, density_weights_m3 (nothing for unitless mode), and
macro_weights (density weights or unit_particle_weight*importance_weights).
The Cartesian position and velocity must use the field coordinate frame.

If an outward normal and source patch area are supplied, also return
rate_weights_s = n*A*max(v dot normal,0)*importance_weight/N [s^-1].
This estimates outward crossing of a reservoir VDF, with NO self-normalization
or rejection of inward samples. Only positive-rate states need tracing.
Each patch needs its own N and local moments. This is not the prescribed
n*norm(U)*g thin-sheet injection model used by backtracing.
No propagation, collisions, escape classification, or files are produced here.
"""
function sample_maxwellian_source(N::Integer; position_m, bulk_velocity_m_s,
        temperature_ev, species="O2+", weights=MonteCarloWeight(),
        rng=Random.default_rng(), normal=nothing, area_m2=nothing)
    N > 0 || throw(ArgumentError("N must be positive"))
    x,U = SVector{3,Float64}(position_m),SVector{3,Float64}(bulk_velocity_m_s)
    all(isfinite,x) && all(isfinite,U) || throw(ArgumentError("nonfinite state"))
    c,n = weights.sampling_temperature_factor, weights.source_number_density_m3
    _mc_positive(c) && _mc_nonnegative(n) && _mc_nonnegative(weights.unit_particle_weight) ||
        throw(ArgumentError("invalid MonteCarloWeight settings"))
    mass = TP.SpeciesDict[species].m
    sigma = thermal_speed_from_temperature_ev(temperature_ev*c,mass)/sqrt(2)
    thermal_speed_from_temperature_ev(temperature_ev,mass)
    (normal === nothing) == (area_m2 === nothing) ||
        throw(ArgumentError("provide both normal and area_m2"))
    er = if normal === nothing
        nothing
    else
        e = SVector{3,Float64}(normal)
        all(isfinite,e) && _mc_positive(norm(e)) && _mc_positive(area_m2) && n>0 ||
            throw(ArgumentError("rate weights require positive density/area and finite nonzero normal"))
        normalize(e)
    end
    states = Vector{SVector{6,Float64}}(undef,N)
    w = Vector{Float64}(undef,N)
    rates = er === nothing ? nothing : zeros(N)
    for i in 1:N
        v = U + sigma*SVector{3,Float64}(randn(rng,3))
        states[i] = SVector{6,Float64}(x...,v...)
        w[i] = maxwellian_importance_weight_3d(U,v,temperature_ev,temperature_ev*c; mass_kg=mass)
        rates === nothing || (rates[i] = n*area_m2*max(dot(v,er),0)*w[i]/N)
    end
    total = sum(w)
    _mc_positive(total) && all(isfinite,w) || error("Importance weights underflowed/overflowed; adjust sampling temperature")
    density = n>0 ? particle_density_weight.(w,n,total) : nothing
    macro_values = density === nothing ? weights.unit_particle_weight*w : density
    return (; initial_states=states, importance_weights=w, density_weights_m3=density,
        macro_weights=macro_values, rate_weights_s=rates, effective_sample_size=total^2/sum(abs2,w),
        species, temperature_ev, sampling_temperature_ev=temperature_ev*c)
end
