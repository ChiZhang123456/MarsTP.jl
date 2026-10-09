Base.@kwdef struct ForwardTraceConfig
    species::String = "O2+"
    mhd_file::String = data_path("mars_fields_spherical_from_dat.vts")
    tspan::Tuple{Float64, Float64} = (0.0, 500.0)
    solver::Symbol = :adaptive
    dt::Float64 = 0.2
    safety::Float64 = 1 / 32
end

"""
    trace_forward(initial_states; config, fields=nothing, backend=nothing, saveat=())

Cartesian states are `[x,y,z,vx,vy,vz]` in m and m/s. Pass cached `MHDFields`
to avoid repeated VTK I/O. A hardware `backend` (e.g. `CUDA.CUDABackend()`)
uses the TestParticle 0.24 ensemble kernel and requires `solver=:boris`.
The return value preserves MarsTP's vector of single-member ensembles.

The backend integrates to the time limit without boundary callbacks. Choose
an in-domain interval; nonfinite saved states raise an error, never extrapolate.
Use `trace_forward_bounded` for per-particle spherical absorption and escape.
CPU backtracing/source accumulation is separate from this forward API.
"""
function trace_forward(initial_states; config::ForwardTraceConfig = ForwardTraceConfig(),
        fields=nothing, backend=nothing, saveat=())
    config.solver in (:boris, :adaptive) || throw(ArgumentError("Unknown solver"))
    backend === nothing || backend isa KA.Backend || throw(ArgumentError("Expected a KernelAbstractions backend"))
    backend === nothing || config.solver == :boris || throw(ArgumentError("GPU/backend tracing requires solver=:boris"))
    all(isfinite, config.tspan) && config.tspan[2] > config.tspan[1] ||
        throw(ArgumentError("Forward tspan must be finite and increasing"))
    config.solver != :boris || isfinite(config.dt) && config.dt > 0 ||
        throw(ArgumentError("Forward Boris dt must be positive and finite"))
    states = [SVector{6,Float64}(s) for s in initial_states]
    isempty(states) && return []
    all(s -> all(isfinite,s), states) || throw(ArgumentError("Nonfinite initial state"))
    fields = fields === nothing ? load_mhd_fields(config.mhd_file; electric_field = :total) : fields
    all(s -> first(fields.r) <= norm(s[1:3]) <= last(fields.r), states) ||
        throw(ArgumentError("Initial position outside radial field domain"))
    param = mhd_param(fields; species = config.species)
    solver = config.solver == :boris ? TP.Boris() : TP.AdaptiveBoris(; safety = config.safety)
    sols = if backend === nothing
        map(states) do state
            prob = TP.TraceProblem(state, config.tspan, param)
            config.solver == :boris ? TP.solve(prob, solver; dt = config.dt, saveat) : TP.solve(prob, solver; saveat)
        end
    else
        prob = TP.TraceProblem(first(states), config.tspan, param)
        # Upstream bulk input layout is (particles, 6), not (6, particles).
        u0 = Matrix{Float64}(undef, length(states), 6)
        for i in eachindex(states), j in 1:6
            u0[i,j] = states[i][j]
        end
        TP.solve(prob, solver, backend; dt=config.dt, trajectories=length(states),
            u0, saveat, save_everystep=isempty(saveat)).u
    end
    for (i,sol) in enumerate(sols)
        sol.retcode == TP.ReturnCode.Success || error("Forward particle $i failed: $(sol.retcode)")
        all(s -> all(isfinite,s), sol.u) || error("Forward particle $i has nonfinite states, possibly left the field domain")
    end
    return [TP.EnsembleSolution([s], 0.0, true) for s in sols]
end
