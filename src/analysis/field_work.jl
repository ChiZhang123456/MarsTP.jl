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
    any(!isfinite, E) && return 0.0
    return q * dot(E, _cartesian_velocity_to_spherical(v, p)) * dt
end

function field_work(sol, itp::FieldWorkInterpolators; species = "O2+")
    traj = hasproperty(sol, :u) && sol.u isa AbstractVector ? sol.u[1] : sol
    q = TP.SpeciesDict[species].q
    total_J = 0.0
    conv_J = 0.0
    hall_J = 0.0
    for i in 1:(length(traj.t) - 1)
        u0, u1 = traj.u[i], traj.u[i + 1]
        p = SA[0.5 * (u0[1] + u1[1]), 0.5 * (u0[2] + u1[2]), 0.5 * (u0[3] + u1[3])]
        v = SA[0.5 * (u0[4] + u1[4]), 0.5 * (u0[5] + u1[5]), 0.5 * (u0[6] + u1[6])]
        dt = traj.t[i + 1] - traj.t[i]
        (any(!isfinite, p) || any(!isfinite, v)) && continue
        total_J += _work_one(itp.total, p, v, q, dt)
        conv_J += _work_one(itp.conv, p, v, q, dt)
        hall_J += _work_one(itp.hall, p, v, q, dt)
    end
    return (; total_J, conv_J, hall_J,
        total_eV = total_J / TP.eV, conv_eV = conv_J / TP.eV, hall_eV = hall_J / TP.eV)
end
