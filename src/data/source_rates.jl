struct O2plusSourceRates
    path::String
    r::Vector{Float64}
    theta::Vector{Float64}
    phi::Vector{Float64}
    production_density::Array{Float64, 3}
    production_rate::Array{Float64, 3}
    outflow_rate::Array{Float64, 2}
    cell_volume::Array{Float64, 3}
    surface_area::Array{Float64, 2}
    metadata::Dict{String, Any}
end

function _mat_value(data, name)
    haskey(data, name) || error("Missing variable $name")
    return data[name]
end

function load_o2plus_source_rates(path = data_path("O2plus_source_rates.mat"); fill_invalid = true)
    resolved = resolve_project_path(path)
    data = matread(resolved)
    production_density = Float64.(_mat_value(data, "O2plus_production_density_m3_inv_s_inv"))
    production_rate = Float64.(_mat_value(data, "O2plus_production_rate_s_inv"))
    outflow_rate = Float64.(_mat_value(data, "O2plus_outflow_rate_s_inv"))
    cell_volume = Float64.(_mat_value(data, "cell_volume_m3"))
    surface_area = Float64.(_mat_value(data, "surface_area_m2"))
    r = Float64.(vec(_mat_value(data, "r_m")))
    theta = Float64.(vec(_mat_value(data, "theta_rad")))
    phi = Float64.(vec(_mat_value(data, "phi_rad")))

    fill_invalid && foreach(A -> A[.!isfinite.(A)] .= 0.0,
        (production_density, production_rate, outflow_rate))

    size(production_density) == (length(r), length(theta), length(phi)) ||
        error("O2+ production density grid does not match coordinates")
    size(production_rate) == size(production_density) ||
        error("O2+ production rate grid does not match density grid")
    size(cell_volume) == size(production_density) ||
        error("Source cell volume grid does not match density grid")
    size(outflow_rate) == (length(theta), length(phi)) ||
        error("O2+ outflow grid does not match angular coordinates")
    size(surface_area) == size(outflow_rate) ||
        error("Surface area grid does not match outflow grid")

    metadata = Dict{String, Any}()
    for key in keys(data)
        value = data[key]
        if !(value isa AbstractArray) || length(value) <= 10
            metadata[String(key)] = value
        end
    end

    return O2plusSourceRates(
        resolved, r, theta, phi, production_density, production_rate,
        outflow_rate, cell_volume, surface_area, metadata,
    )
end

function assert_same_grid(source::O2plusSourceRates, fields::MHDFields)
    isapprox(source.r, fields.r; rtol = 1.0e-8, atol = 1.0e-6) ||
        error("Source-rate radial grid does not match MHD grid")
    isapprox(source.theta, fields.theta; rtol = 1.0e-12, atol = 1.0e-12) ||
        error("Source-rate theta grid does not match MHD grid")
    isapprox(source.phi, fields.phi; rtol = 1.0e-12, atol = 1.0e-12) ||
        error("Source-rate phi grid does not match MHD grid")
    return true
end
