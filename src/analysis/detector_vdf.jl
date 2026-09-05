function load_detector_vdf(path)
    data = load(resolve_project_path(path))
    f = data["f2d_xz"]
    positive = f[isfinite.(f) .& (f .> 0)]
    return (; vx_km = data["vx_km"], vz_km = data["vz_km"], f2d_xz = f,
        max_f = isempty(positive) ? NaN : maximum(positive),
        positive_count = length(positive), units = get(data, "units", "unknown"))
end

function plot_vx_vz(vdf; output = nothing, color_min = nothing, color_max = nothing)
    f = copy(vdf.f2d_xz)
    valid = isfinite.(f) .& (f .> 0)
    logs = fill(isnothing(color_min) ? -20.0 : color_min, size(f))
    logs[valid] .= log10.(f[valid])
    cmax = isnothing(color_max) ? maximum(logs[valid]; init = -8.0) : color_max
    cmin = isnothing(color_min) ? cmax - 8.0 : color_min
    fig = Figure(size = (900, 700))
    ax = Axis(fig[1, 1], xlabel = "Vx [km/s]", ylabel = "Vz [km/s]",
        title = "Backtraced O2+ f(Vx,Vz)")
    hm = heatmap!(ax, vdf.vx_km, vdf.vz_km, logs; colormap = :turbo,
        colorrange = (cmin, cmax))
    Colorbar(fig[1, 2], hm, label = "log10 f")
    isnothing(output) || save(output, fig, px_per_unit = 2)
    return fig
end
