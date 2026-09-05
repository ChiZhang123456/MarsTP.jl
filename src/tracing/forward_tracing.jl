Base.@kwdef struct ForwardTraceConfig
    species::String = "O2+"
    mhd_file::String = data_path("mars_fields_spherical_from_dat.vts")
    tspan::Tuple{Float64, Float64} = (0.0, 500.0)
    solver::Symbol = :adaptive
    dt::Float64 = 0.2
    safety::Float64 = 1 / 32
end

function trace_forward(initial_states; config::ForwardTraceConfig = ForwardTraceConfig())
    fields = load_mhd_fields(config.mhd_file; electric_field = :total)
    param = mhd_param(fields; species = config.species)
    solver = config.solver == :boris ? TP.Boris() : TP.AdaptiveBoris(; safety = config.safety)
    return map(initial_states) do state
        prob = TP.TraceProblem(SVector{6, Float64}(state), config.tspan, param)
        config.solver == :boris ? TP.solve(prob, solver; dt = config.dt) : TP.solve(prob, solver)
    end
end
