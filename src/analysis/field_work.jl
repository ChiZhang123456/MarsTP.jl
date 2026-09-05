struct FieldWorkInterpolators
    total::Any
    conv::Any
    hall::Any
end

function build_field_work_interpolators(mhd_file = data_path("mars_fields_spherical_from_dat.vts"))
    total = load_mhd_fields(mhd_file; electric_field = :total)
    conv = load_mhd_fields(mhd_file; electric_field = :conv)
    hall = load_mhd_fields(mhd_file; electric_field = :hall)
    return FieldWorkInterpolators(
        TP.build_interpolator(TP.StructuredGrid, total.E, total.r, total.theta, total.phi),
        TP.build_interpolator(TP.StructuredGrid, conv.E, conv.r, conv.theta, conv.phi),
        TP.build_interpolator(TP.StructuredGrid, hall.E, hall.r, hall.theta, hall.phi),
    )
end

function _cartesian_velocity_to_spherical(v, p)
    x, y, z = p
    rxy = hypot(x, y)
    r = hypot(rxy, z)
    sin_theta = rxy / r
    cos_theta = z / r
    sin_phi = rxy == 0 ? 0.0 : y / rxy
    cos_phi = rxy == 0 ? 1.0 : x / rxy
    vx, vy, vz = v
    return SA[
        vx * sin_theta * cos_phi + vy * sin_theta * sin_phi + vz * cos_theta,
        vx * cos_theta * cos_phi + vy * cos_theta * sin_phi - vz * sin_theta,
        -vx * sin_phi + vy * cos_phi,
    ]
end

function _work_one(itp, p, v, q, dt)
    E = itp(p)
    all(isfinite, E) || error("Nonfinite electric field in work integral")
    # StructuredGrid vector interpolators already return Cartesian components.
    return q * dot(E, v) * dt
end

function _work_trajectory(sol)
    if hasproperty(sol, :t)
        traj = sol
    elseif hasproperty(sol, :u) && length(sol.u) == 1
        return _work_trajectory(only(sol.u))
    else
        throw(ArgumentError("Expected one trajectory; use particle_field_work for multiple particles"))
    end
    length(traj.t) == length(traj.u) && !isempty(traj.t) ||
        throw(ArgumentError("Trajectory must have matching nonempty time and state arrays"))
    all(isfinite, traj.t) || throw(ArgumentError("Nonfinite trajectory times"))
    dt = diff(traj.t)
    (all(>(0), dt) || all(<(0), dt)) ||
        throw(ArgumentError("Trajectory times must be strictly monotonic"))
    all(u -> length(u) >= 6 && all(isfinite, u[1:6]), traj.u) ||
        throw(ArgumentError("States must contain finite Cartesian position and velocity"))
    return traj
end

"""
    field_work(sol, itp; species="O2+")

Signed total, convection and Hall work along one existing trajectory, in J and eV.
Uses midpoint quadrature of q E⋅v dt. Positions and velocities must be Cartesian
SI values; fields returned by `itp` must be Cartesian V/m. Decreasing times
retain their sign. A single-member ensemble is accepted; larger ensembles
must use `particle_field_work`. Does not write files.
"""
function field_work(sol, itp::FieldWorkInterpolators; species = "O2+")
    traj = _work_trajectory(sol)
    q = TP.SpeciesDict[species].q
    total_J = 0.0
    conv_J = 0.0
    hall_J = 0.0
    for i in 1:(length(traj.t) - 1)
        u0, u1 = traj.u[i], traj.u[i + 1]
        p = SA[0.5 * (u0[1] + u1[1]), 0.5 * (u0[2] + u1[2]), 0.5 * (u0[3] + u1[3])]
        v = SA[0.5 * (u0[4] + u1[4]), 0.5 * (u0[5] + u1[5]), 0.5 * (u0[6] + u1[6])]
        dt = traj.t[i + 1] - traj.t[i]
        all(isfinite, p) && all(isfinite, v) || error("Nonfinite trajectory in work integral")
        total_J += _work_one(itp.total, p, v, q, dt)
        conv_J += _work_one(itp.conv, p, v, q, dt)
        hall_J += _work_one(itp.hall, p, v, q, dt)
    end
    return (; total_J, conv_J, hall_J,
        total_eV = total_J / TP.eV, conv_eV = conv_J / TP.eV, hall_eV = hall_J / TP.eV)
end

"""
    Electric_field_work_profile(sol, itp; species="O2+")

Return time-resolved signed work and instantaneous power on the same trajectory.
`total_eV`, `conv_eV`, `hall_eV` are cumulative work, starting at zero.
`power_total_eV_s`, `power_conv_eV_s`, `power_hall_eV_s` are q E⋅v at saved states.
`kinetic_eV` and `delta_kinetic_eV` use nonrelativistic kinetic energy.
`energy_residual_eV = delta_kinetic_eV - total_eV` checks quadrature/integration
accuracy when only electromagnetic forces act. `field_sum_residual_eV` checks
whether total work equals convection plus Hall work. Totals are in `summary`.
Invalid fields/states raise an error instead of silently contributing zero.
Coarse saved trajectories can underestimate work; check saving cadence and dt.
"""
function Electric_field_work_profile(sol, itp::FieldWorkInterpolators; species = "O2+")
    traj = _work_trajectory(sol)
    sp = TP.SpeciesDict[species]
    n = length(traj.t)
    work = zeros(n, 3)
    power = zeros(n, 3)
    kinetic = zeros(n)
    fields = (itp.total, itp.conv, itp.hall)
    for i in 1:n
        u = traj.u[i]
        p, v = SA[u[1], u[2], u[3]], SA[u[4], u[5], u[6]]
        kinetic[i] = 0.5 * sp.m * dot(v, v) / TP.eV
        for j in 1:3
            power[i, j] = _work_one(fields[j], p, v, sp.q, 1.0) / TP.eV
        end
        if i > 1
            a = traj.u[i - 1]
            pm = 0.5 * (p + SA[a[1], a[2], a[3]])
            vm = 0.5 * (v + SA[a[4], a[5], a[6]])
            dt = traj.t[i] - traj.t[i - 1]
            for j in 1:3
                work[i, j] = work[i - 1, j] + _work_one(fields[j], pm, vm, sp.q, dt) / TP.eV
            end
        end
    end
    total_eV, conv_eV, hall_eV = work[:, 1], work[:, 2], work[:, 3]
    delta_kinetic_eV = kinetic .- first(kinetic)
    energy_residual_eV = delta_kinetic_eV .- total_eV
    field_sum_residual_eV = total_eV .- conv_eV .- hall_eV
    summary = (; total_eV = last(total_eV), conv_eV = last(conv_eV),
        hall_eV = last(hall_eV), delta_kinetic_eV = last(delta_kinetic_eV),
        energy_residual_eV = last(energy_residual_eV),
        field_sum_residual_eV = last(field_sum_residual_eV))
    return (; t = collect(traj.t), total_eV, conv_eV, hall_eV,
        power_total_eV_s = power[:, 1], power_conv_eV_s = power[:, 2],
        power_hall_eV_s = power[:, 3], kinetic_eV = kinetic, delta_kinetic_eV,
        energy_residual_eV, field_sum_residual_eV, summary)
end

"""
    particle_field_work(solutions, itp; species="O2+", profiles=false, threaded=false)

Analyze every particle in a TestParticle ensemble or a vector of trajectories
(including the vector of single-member ensembles returned by `trace_forward`).
Returns results in input order, using `field_work` or, if `profiles=true`,
`Electric_field_work_profile`. Build interpolators once and reuse them. Threaded mode
requires thread-safe, read-only interpolators. No files are saved.
"""
function particle_field_work(solutions, itp::FieldWorkInterpolators;
        species = "O2+", profiles = false, threaded = false)
    hasproperty(solutions, :t) &&
        throw(ArgumentError("Use field_work or Electric_field_work_profile for a single trajectory"))
    trajectories = solutions isa AbstractVector ? solutions : solutions.u
    analyze = profiles ? Electric_field_work_profile : field_work
    result = Vector{Any}(undef, length(trajectories))
    if threaded
        Threads.@threads for i in eachindex(trajectories)
            result[i] = analyze(trajectories[i], itp; species)
        end
    else
        for i in eachindex(trajectories)
            result[i] = analyze(trajectories[i], itp; species)
        end
    end
    return result
end

"""Compatibility alias for [`Electric_field_work_profile`](@ref)."""
const field_work_profile = Electric_field_work_profile
