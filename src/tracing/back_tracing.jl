Base.@kwdef struct BacktraceConfig
    detector_Rm::SVector{3, Float64} = SA[-1.5, 0.0, 1.0]
    species::String = "O2+"
    mhd_file::String = data_path("mars_fields_spherical_from_dat.vts")
    source_file::String = data_path("O2plus_source_rates.mat")
    vx_min_kms::Float64 = -500.0
    vx_max_kms::Float64 = 500.0
    vy_min_kms::Float64 = -500.0
    vy_max_kms::Float64 = 500.0
    vz_min_kms::Float64 = -500.0
    vz_max_kms::Float64 = 500.0
    dv_kms::Float64 = 1.0
    dvy_kms::Float64 = 1.0
    tspan::Tuple{Float64, Float64} = (0.0, -500.0)
    solver::Symbol = :boris
    dt::Float64 = -0.2
    safety::Float64 = 0.005
    stream_vy::Bool = true
    checkpoint_file::String = project_path("output", "backtrace_checkpoint.jld2")
    output_file::String = project_path("output", "backtrace_vdf.jld2")
    progress_interval::Int = 10000
end

velocity_axes(c::BacktraceConfig) = (
    vx = collect((c.vx_min_kms * 1e3):(c.dv_kms * 1e3):(c.vx_max_kms * 1e3)),
    vy = collect((c.vy_min_kms * 1e3):(c.dvy_kms * 1e3):(c.vy_max_kms * 1e3)),
    vz = collect((c.vz_min_kms * 1e3):(c.dv_kms * 1e3):(c.vz_max_kms * 1e3)),
)

function _integrate_source(sol, q_itp, gitm, mass, r_inner, r_outer)
    f = 0.0
    status = 1
    traj = hasproperty(sol, :u) && sol.u isa AbstractVector ? sol.u[1] : sol
    for i in 1:(length(traj.t) - 1)
        u = traj.u[i]
        p = SVector{3, Float64}(u[1], u[2], u[3])
        v = SVector{3, Float64}(u[4], u[5], u[6])
        any(!isfinite, p) && return f, 3
        r = norm(p)
        r < r_inner && return f, 2
        r > r_outer && return f, 4
        any(!isfinite, v) && return f, 3

        q = q_itp(p)
        Tn = neutral_properties(gitm, p).Tn
        if isfinite(q) && q > 0 && isfinite(Tn) && Tn > 0
            vth = sqrt(2 * TP.kB * Tn / mass)
            g = VDF.Maxwellian(vth; u0 = SA[0.0, 0.0, 0.0])(v)
            f += q * g * abs(traj.t[i + 1] - traj.t[i])
        end
    end
    return f, status
end

function _solve_one_backtrace(position, velocity, tspan, param, solver, config)
    prob = TP.TraceProblem(vcat(position, velocity), tspan, param)
    config.solver == :boris && return TP.solve(prob, solver; dt = config.dt)
    return TP.solve(prob, solver)
end

function run_backtrace_vdf(config::BacktraceConfig = BacktraceConfig())
    config.species == "O2+" || error("Only O2+ source-rate backtracing is configured")
    mkpath(dirname(config.output_file))
    fields = load_mhd_fields(config.mhd_file; electric_field = :total)
    source = load_o2plus_source_rates(config.source_file)
    assert_same_grid(source, fields)
    gitm = load_gitm()
    param = mhd_param(fields; species = config.species)
    q_itp = TP.build_interpolator(TP.StructuredGrid, source.production_density,
        fields.r, fields.theta, fields.phi)
    axes = velocity_axes(config)
    nx, ny, nz = length(axes.vx), length(axes.vy), length(axes.vz)
    total = nx * ny * nz
    total > 100_000_000 && !config.stream_vy &&
        error("Use stream_vy=true for very large 3D velocity grids")

    position = config.detector_Rm .* Rm
    mass = TP.SpeciesDict[config.species].m
    solver = config.solver == :boris ? TP.Boris() : TP.AdaptiveBoris(; safety = config.safety)
    f2d = zeros(Float64, nx, nz)
    status_counts = zeros(Int, 5)
    done = Threads.Atomic{Int}(0)

    for (iy, vy) in pairs(axes.vy)
        slice = zeros(Float64, nx, nz)
        flags = zeros(Int, nx, nz)
        Threads.@threads for index in eachindex(slice)
            I = CartesianIndices(slice)[index]
            velocity = SA[axes.vx[I[1]], vy, axes.vz[I[2]]]
            sol = _solve_one_backtrace(position, velocity, config.tspan, param, solver, config)
            slice[I], flags[I] = _integrate_source(sol, q_itp, gitm, mass, Rinner, Router)
            n = Threads.atomic_add!(done, 1) + 1
            config.progress_interval > 0 && n % config.progress_interval == 0 &&
                println("backtrace progress: $n / $total")
        end
        f2d .+= slice .* (config.dvy_kms * 1e3)
        for flag in 1:4
            status_counts[flag + 1] += count(==(flag), flags)
        end
        if config.stream_vy
            jldsave(config.checkpoint_file; completed_iy = iy, f2d_xz = f2d,
                status_counts, vx_km = axes.vx ./ 1e3, vy_km = axes.vy ./ 1e3,
                vz_km = axes.vz ./ 1e3, config)
        end
    end

    result = (; vx_km = axes.vx ./ 1e3, vy_km = axes.vy ./ 1e3,
        vz_km = axes.vz ./ 1e3, f2d_xz = f2d, status_counts, config,
        units = "s^2 m^-5")
    jldsave(config.output_file; result...)
    return result
end
