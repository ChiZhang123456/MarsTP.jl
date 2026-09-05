struct MHDFields
    r::Vector{Float64}
    theta::Vector{Float64}
    phi::Vector{Float64}
    E::Array{Float64, 4}
    B::Array{Float64, 4}
    electric_field::Symbol
    path::String
end

struct SpeciesMoments
    n::Array{Float64, 3}
    v::Array{Float64, 4}
    t::Array{Float64, 3}
end

function _mhd_axes(nr, ntheta, nphi)
    return (
        collect(exp.(range(log(Rinner), log(Router), length = nr))),
        collect(range(0.0, pi, length = ntheta)),
        collect(range(0.0, 2pi, length = nphi)),
    )
end

function _point_vector(point_data, names; scale = 1.0)
    keys_available = keys(point_data)
    for name in names
        if name in keys_available
            data = get_data(point_data[name])
            return scale == 1.0 ? data : data .* scale
        end
    end
    error("Missing vector field. Tried: $(join(names, ", "))")
end

function _electric_names(which)
    which in (:total, "total") && return ("E_Total [V/m]", "E_total [V/m]", "E [V/m]")
    which in (:conv, :convection, "conv", "convection") && return ("E_conv [V/m]", "Econv")
    which in (:hall, "hall") && return ("E_hall [V/m]",)
    which isa AbstractString && return (which,)
    error("Unknown electric field selector: $which")
end

function _rotate_vectors_to_spherical!(A, r, theta, phi)
    for k in eachindex(phi), j in eachindex(theta)
        sin_theta, cos_theta = sincos(theta[j])
        sin_phi, cos_phi = sincos(phi[k])
        for i in eachindex(r)
            ax, ay, az = A[:, i, j, k]
            A[1, i, j, k] = ax * sin_theta * cos_phi + ay * sin_theta * sin_phi + az * cos_theta
            A[2, i, j, k] = ax * cos_theta * cos_phi + ay * cos_theta * sin_phi - az * sin_theta
            A[3, i, j, k] = -ax * sin_phi + ay * cos_phi
        end
    end
    return A
end

function load_mhd_fields(
        path = data_path("mars_fields_spherical_from_dat.vts");
        electric_field = :total,
    )
    resolved = resolve_project_path(path)
    vtk = VTKFile(resolved)
    point_data = get_point_data(vtk)
    dims = ReadVTK.get_wholeextent(vtk.xml_file)[1]
    nr, ntheta, nphi = dims
    r, theta, phi = _mhd_axes(nr, ntheta, nphi)

    B_raw = if "B_Field [T]" in keys(point_data)
        get_data(point_data["B_Field [T]"])
    elseif "B [nT]" in keys(point_data)
        get_data(point_data["B [nT]"]) .* 1.0e-9
    elseif "B" in keys(point_data)
        get_data(point_data["B"]) .* 1.0e-9
    else
        error("Magnetic field not found in $resolved")
    end
    E_raw = _point_vector(point_data, _electric_names(electric_field))

    B = Float64.(reshape(B_raw, 3, nr, ntheta, nphi))
    E = Float64.(reshape(E_raw, 3, nr, ntheta, nphi))
    _rotate_vectors_to_spherical!(B, r, theta, phi)
    _rotate_vectors_to_spherical!(E, r, theta, phi)

    return MHDFields(r, theta, phi, E, B, Symbol(electric_field), resolved)
end

function _read_array(point_data, name, nr, ntheta, nphi)
    values = get_data(point_data[name])
    if length(values) == nr * ntheta * nphi
        return Float64.(reshape(values, nr, ntheta, nphi))
    elseif length(values) == 3 * nr * ntheta * nphi
        return Float64.(reshape(values, 3, nr, ntheta, nphi))
    end
    error("Unexpected array size for $name")
end

function load_mhd_moments(path = data_path("mars_fields_spherical_from_dat.vts"))
    resolved = resolve_project_path(path)
    vtk = VTKFile(resolved)
    point_data = get_point_data(vtk)
    nr, ntheta, nphi = ReadVTK.get_wholeextent(vtk.xml_file)[1]
    r, theta, phi = _mhd_axes(nr, ntheta, nphi)

    moments = (
        Hplus = SpeciesMoments(
            _read_array(point_data, "n_H^p [m^-3]", nr, ntheta, nphi),
            _read_array(point_data, "U_H^p [m/s]", nr, ntheta, nphi),
            _read_array(point_data, "T_H^p [K]", nr, ntheta, nphi),
        ),
        Oplus = SpeciesMoments(
            _read_array(point_data, "n_O^p^p [m^-3]", nr, ntheta, nphi),
            _read_array(point_data, "U_O^p^p [m/s]", nr, ntheta, nphi),
            _read_array(point_data, "T_O^p^p [K]", nr, ntheta, nphi),
        ),
        O2plus = SpeciesMoments(
            _read_array(point_data, "n_O^2^p [m^-3]", nr, ntheta, nphi),
            _read_array(point_data, "U_O^2^p [m/s]", nr, ntheta, nphi),
            _read_array(point_data, "T_O^2^p [K]", nr, ntheta, nphi),
        ),
        CO2plus = SpeciesMoments(
            _read_array(point_data, "n_C^O^2^p [m^-3]", nr, ntheta, nphi),
            _read_array(point_data, "U_C^O^2^p [m/s]", nr, ntheta, nphi),
            _read_array(point_data, "T_C^O^2^p [K]", nr, ntheta, nphi),
        ),
        Te = _read_array(point_data, "T_e [K]", nr, ntheta, nphi),
    )
    return (; r, theta, phi, moments...)
end

function mhd_param(fields::MHDFields; species = "O2+")
    return TP.prepare(
        fields.r,
        fields.theta,
        fields.phi,
        fields.E,
        fields.B;
        species = TP.SpeciesDict[species],
        gridtype = TP.StructuredGrid,
    )
end
