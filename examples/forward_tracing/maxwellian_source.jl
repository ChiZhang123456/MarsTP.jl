using MarsTP, Random, LinearAlgebra

# One local O2+ reservoir patch, Cartesian field frame, area in m^2.
# Choose the patch geometry/area explicitly for your source discretization.
function sample_outflow_patch(source, position_m, area_m2; N=100, seed=42)
    p = ionosphere_properties(source,position_m)
    return sample_maxwellian_source(N; position_m=position_m/norm(position_m)*source.radius,
        bulk_velocity_m_s=p.Ui, temperature_ev=MarsTP.TP.kB*p.Ti/1.602176634e-19,
        weights=MonteCarloWeight(source_number_density_m3=p.n,sampling_temperature_factor=4.),
        normal=position_m,area_m2, rng=MersenneTwister(seed))
end

# source = load_ionosphere_source(;altitude_km=400.)
# sample = sample_outflow_patch(source,[Rm+400e3,0.,0.],1e6)
# ids = findall(>(0),sample.rate_weights_s)
# trajectories = trace_forward(sample.initial_states[ids];
#     config=ForwardTraceConfig(species=sample.species,tspan=(0.,1.)))
# Retain sample.rate_weights_s[ids] paired with those trajectories.
# Set and validate domain/impact/time-limit termination before escape studies.
