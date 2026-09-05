Base.@kwdef struct BacktraceConfig
    detector_Rm::SVector{3, Float64} = SA[-1.5, 0.0, 1.0]
    species::String = "O2+"
    ionosphere_altitude_km::Float64 = 200.0
    include_ionosphere::Bool = true
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

# First entry into the inner sphere or exit from the outer sphere along a
# Boris position segment. Also catches a segment passing through the sphere.
function _backtrace_crossing(a, b, inner, outer)
    d = b-a
    A = dot(d,d)
    A == 0 && return (1.0, 1)
    best, status = 1.0, 1
    for (radius, flag) in ((inner,2), (outer,4))
        B, C = 2dot(a,d), dot(a,a)-radius^2
        disc = B^2-4A*C
        disc < 0 && continue
        for fraction in ((-B-sqrt(disc))/(2A), (-B+sqrt(disc))/(2A))
            0 <= fraction <= best || continue
            slope = dot(a+fraction*d,d)
            ((flag == 2 && slope < 0) || (flag == 4 && slope > 0)) || continue
            best, status = fraction, flag
        end
    end
    return best, status
end

function _volume_source(q_itp, gitm, mass, position, velocity)
    q = q_itp(position)
    isfinite(q) && q >= 0 || error("Invalid volume production at $position: $q")
    q == 0 && return 0.0
    Tn = neutral_properties(gitm, position).Tn
    isfinite(Tn) && Tn > 0 || error("Invalid neutral temperature at $position: $Tn")
    return q * VDF.Maxwellian(sqrt(2*TP.kB*Tn/mass); u0=SA[0.,0.,0.])(velocity)
end

# Stream quadrature in the callback. TestParticle supplies a staggered velocity;
# synchronize it at the accepted endpoint (or shell intersection) before using f.
function _trace_sources(position, velocity, param, config, volume_source, ionosphere;
        outer=Router)
    inner = config.include_ionosphere ? _ionosphere_radius(config.ionosphere_altitude_km) : Rinner
    inner < norm(position) < outer ||
        throw(ArgumentError("Detector must lie strictly between source/inner shell and outer boundary"))
    previous = Ref(position)
    previous_t = Ref(config.tspan[1])
    previous_q = Ref(volume_source(position, velocity))
    volume, boundary_f, flux = Ref(0.0), Ref(0.0), Ref(0.0)
    status = Ref(1)
    endpoint, end_velocity = Ref(position), Ref(velocity)
    end_time = Ref(config.tspan[1])
    function boundary(u, p, t)
        if !all(isfinite, u) || !isfinite(t)
            status[] = 3
            return true
        end
        step = t-previous_t[]
        step < 0 || error("Backtrace time must decrease")
        r, vhalf = SVector{3,Float64}(u[1:3]), SVector{3,Float64}(u[4:6])
        fraction, flag = _backtrace_crossing(previous[], r, inner, outer)
        x = previous[] + fraction*(r-previous[])
        hit_t = previous_t[] + fraction*step
        # Do not evaluate fields at the rejected out-of-domain endpoint.
        safe_x = x / norm(x) * clamp(norm(x), Rinner+1e-8, outer-1e-8)
        v = TP.update_velocity(vhalf, safe_x, (fraction-0.5)*step, hit_t, p)
        all(isfinite, v) || error("Invalid synchronized velocity at $safe_x")
        q = volume_source(safe_x, v)
        volume[] += 0.5*(previous_q[]+q)*abs(fraction*step)
        endpoint[], end_velocity[], end_time[] = x, v, hit_t
        if flag != 1
            status[] = flag
            if flag == 2 && config.include_ionosphere
                result = ionosphere(x,v)
                boundary_f[], flux[] = result.f, result.flux
            end
            return true
        end
        previous[], previous_t[], previous_q[] = r, t, q
        return false
    end
    prob = TP.TraceProblem(vcat(position, velocity), config.tspan, param)
    if config.solver == :boris
        # Fit an integer number of steps to the span, avoiding the solver's
        # nearest-integer step count leaving a short unintegrated final segment.
        steps = ceil(Int, abs((config.tspan[2]-config.tspan[1])/config.dt))
        dt = (config.tspan[2]-config.tspan[1])/steps
        TP.solve(prob, TP.Boris(); dt, isoutside=boundary, maxiters=steps+1,
            save_start=false, save_end=false, save_everystep=false)
    else
        TP.solve(prob, TP.AdaptiveBoris(; safety=config.safety); isoutside=boundary,
            save_start=false, save_end=false, save_everystep=false)
    end
    if status[] == 1 && !isapprox(end_time[], config.tspan[2]; atol=1e-7, rtol=1e-10)
        error("Backtrace stopped before the requested time without a boundary status")
    end
    return (; volume=volume[], ionosphere=boundary_f[], flux=flux[], status=status[],
        position=endpoint[], velocity=end_velocity[], time=end_time[])
end

function _validate_backtrace(config)
    radius = _ionosphere_radius(config.ionosphere_altitude_km)
    inner = config.include_ionosphere ? radius : Rinner
    all(isfinite, config.detector_Rm) && inner < norm(config.detector_Rm)*Rm < Router ||
        throw(ArgumentError("Detector must lie strictly above the source shell and below 4 Rm"))
    config.solver in (:boris,:adaptive) || throw(ArgumentError("solver must be :boris or :adaptive"))
    all(isfinite, config.tspan) && config.tspan[2] < config.tspan[1] ||
        throw(ArgumentError("tspan must be finite and decreasing"))
    isfinite(config.dt) && config.dt < 0 || throw(ArgumentError("dt must be finite and negative"))
    isfinite(config.safety) && config.safety > 0 || throw(ArgumentError("safety must be positive"))
    for (lo,hi,dv) in ((config.vx_min_kms,config.vx_max_kms,config.dv_kms),
            (config.vy_min_kms,config.vy_max_kms,config.dvy_kms),
            (config.vz_min_kms,config.vz_max_kms,config.dv_kms))
        all(isfinite,(lo,hi,dv)) && lo <= hi && dv > 0 ||
            throw(ArgumentError("Velocity axes require finite min <= max and positive spacing"))
    end
    return nothing
end

function run_backtrace_vdf(config::BacktraceConfig = BacktraceConfig())
    config.species == "O2+" || error("Only O2+ source-rate backtracing is configured")
    _validate_backtrace(config)
    mkpath(dirname(config.output_file))
    fields = load_mhd_fields(config.mhd_file; electric_field = :total)
    source = load_o2plus_source_rates(config.source_file; fill_invalid=false)
    ionosphere_source = config.include_ionosphere ? load_ionosphere_source(config.mhd_file;
        altitude_km=config.ionosphere_altitude_km) : nothing
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
    volume_source(p,v) = _volume_source(q_itp, gitm, mass, p,v)
    boundary_source(p,v) = ionosphere_distribution(ionosphere_source,p,v; mass)
    f2d = zeros(Float64, nx, nz)
    f2d_volume = zeros(Float64, nx, nz)
    f2d_ionosphere = zeros(Float64, nx, nz)
    status_counts = zeros(Int, 5)
    done = Threads.Atomic{Int}(0)

    for (iy, vy) in pairs(axes.vy)
        slice = zeros(Float64, nx, nz)
        boundary_slice = zeros(Float64, nx, nz)
        flags = zeros(Int, nx, nz)
        Threads.@threads for index in eachindex(slice)
            I = CartesianIndices(slice)[index]
            velocity = SA[axes.vx[I[1]], vy, axes.vz[I[2]]]
            result = _trace_sources(position, velocity, param, config, volume_source,
                boundary_source; outer=last(fields.r))
            slice[I], boundary_slice[I], flags[I] = result.volume, result.ionosphere, result.status
            n = Threads.atomic_add!(done, 1) + 1
            config.progress_interval > 0 && n % config.progress_interval == 0 &&
                println("backtrace progress: $n / $total")
        end
        f2d_volume .+= slice .* (config.dvy_kms * 1e3)
        f2d_ionosphere .+= boundary_slice .* (config.dvy_kms * 1e3)
        f2d .= f2d_volume .+ f2d_ionosphere
        for flag in 1:4
            status_counts[flag + 1] += count(==(flag), flags)
        end
        if config.stream_vy
            jldsave(config.checkpoint_file; completed_iy = iy, f2d_xz = f2d, f2d_volume, f2d_ionosphere,
                status_counts, vx_km = axes.vx ./ 1e3, vy_km = axes.vy ./ 1e3,
                vz_km = axes.vz ./ 1e3, config)
        end
    end

    result = (; vx_km = axes.vx ./ 1e3, vy_km = axes.vy ./ 1e3,
        vz_km = axes.vz ./ 1e3, f2d_xz = f2d, f2d_volume, f2d_ionosphere, status_counts, config,
        status_labels = ["unused", "time_limit", "inner_or_ionosphere", "nonfinite", "outer"],
        ionosphere_flux_definition = "n*norm(Ui) [m^-2 s^-1]; diagnostic, not a VDF multiplier",
        units = "s^2 m^-5")
    jldsave(config.output_file; result...)
    return result
end
