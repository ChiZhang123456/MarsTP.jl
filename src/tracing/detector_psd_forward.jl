"""
    forward_psd(solutions; detector_m, side_m, vlim, vgrid, species="O2+",
                rate_weights_s, option="3D", velocity_unit=:m_s,
                coordinate_system="unspecified")

Steady-source, finite-volume velocity PSD from saved forward trajectories.
`solutions` accepts `trace_forward` output, a TestParticle ensemble, a vector of
trajectories, or one trajectory with `t` and `u` ([x,y,z,vx,vy,vz,...]). Saved
positions/velocities/times must be synchronized Cartesian m, m/s, s. The cube
is axis-aligned in that same frame; `detector_m` is its center, `side_m` its side.

`rate_weights_s` is required: one nonnegative physical particle rate per input
trajectory, in particles/s, already including source Monte Carlo normalization.
Density or dimensionless weights cannot substitute for rates. All trajectories
must belong to `species`; known TestParticle mass/charge parameters are checked.

`vlim` is a positive symmetric limit, `(min,max)` shared by all axes, or three
`(min,max)` pairs. `vgrid` is the number of bins (positive integer), or three
bin counts. `velocity_unit=:m_s` (default) or `:km_s` applies ONLY to `vlim`.
Returned velocity centers/edges are always m/s; returned PSD is always SI.
`option` is "3D", "Vx-Vy", "Vy-Vz", or "Vx-Vz" (also :xyz/:xy/:yz/:xz).
Two-dimensional output integrates the omitted axis over its specified range.

Returns `psd`, `axes`, `velocity_centers_m_s`, `velocity_edges_m_s`, `units`,
`density_in_range_m3`, `density_outside_vlim_m3`, `density_total_m3`, and
per-trajectory residence times and saved termination codes. Array dimension
order follows `axes`. 3D units are s^3 m^-6, 2D units are s^2 m^-5.

Accumulates Q*dt within the cube AND velocity bin. Each saved segment is linear
in position and velocity and is split at cube faces and velocity-bin edges.
Bins and cube use [lower,upper), including the final velocity upper edge.
This handles crossings with both endpoints outside. Curvature lost between
saved samples cannot be recovered: check timestep/save-cadence convergence.
Only saved time intervals contribute; no extrapolation or extra time/sample
normalization is performed. Failed solver return codes and invalid data error.
No tracing, field loading, plotting, or file writing is performed.
"""
function forward_psd(solutions; detector_m, side_m, vlim, vgrid,
        species="O2+", rate_weights_s, option="3D", velocity_unit=:m_s,
        coordinate_system="unspecified")
    haskey(TP.SpeciesDict, species) || throw(ArgumentError("Unknown species: $species"))
    sp = TP.SpeciesDict[species]
    length(detector_m) == 3 && all(isfinite, detector_m) ||
        throw(ArgumentError("detector_m must contain three finite coordinates"))
    isfinite(side_m) && side_m > 0 || throw(ArgumentError("side_m must be positive"))
    center = SVector{3,Float64}(detector_m)
    side = Float64(side_m)
    lower, upper = center .- side/2, center .+ side/2
    all(isfinite, lower) && all(isfinite, upper) && all(upper .> lower) &&
        isfinite(side^3) && side^3 > 0 || throw(ArgumentError("Unrepresentable cube"))
    edges = _psd_edges(vlim, vgrid, velocity_unit)
    widths = diff.(edges)
    dims = length.(widths)
    key = lowercase(replace(string(option), "v"=>"", "V"=>"", "-"=>""))
    kept = key in ("3d", "xyz") ? (1,2,3) : key == "xy" ? (1,2) :
        key == "yz" ? (2,3) : key == "xz" ? (1,3) :
        throw(ArgumentError("option must be 3D, Vx-Vy, Vy-Vz, or Vx-Vz"))
    trajectories = hasproperty(solutions, :t) ? (solutions,) :
        solutions isa AbstractVector ? solutions : solutions.u
    rate_weights_s isa AbstractVector || rate_weights_s isa Tuple ||
        throw(ArgumentError("rate_weights_s must contain one rate per trajectory"))
    length(trajectories) == length(rate_weights_s) ||
        throw(ArgumentError("One rate_weights_s entry is required per trajectory"))
    rates = Float64[q for q in rate_weights_s]
    all(q -> isfinite(q) && q >= 0, rates) || throw(ArgumentError("Invalid particle rates"))
    occupancy = zeros(ntuple(j -> dims[kept[j]], length(kept)))
    residence = zeros(length(trajectories))
    outside = zeros(length(trajectories))
    retcodes = String[]
    cuts = Float64[]
    for (i, wrapped) in enumerate(trajectories)
        traj = _psd_trajectory(wrapped, sp)
        push!(retcodes, hasproperty(traj, :retcode) ? string(traj.retcode) : "unavailable")
        for j in 1:(length(traj.t)-1)
            a, b = traj.u[j], traj.u[j+1]
            dt = Float64(traj.t[j+1]) - Float64(traj.t[j])
            isfinite(dt) && dt > 0 || throw(ArgumentError("Forward times must strictly increase"))
            x = SVector{3,Float64}(a[1],a[2],a[3])
            dx = SVector{3,Float64}(b[1],b[2],b[3]) - x
            v = SVector{3,Float64}(a[4],a[5],a[6])
            dv = SVector{3,Float64}(b[4],b[5],b[6]) - v
            all(isfinite, dx) && all(isfinite, dv) || throw(ArgumentError("Segment overflow"))
            lo, hi = 0.0, 1.0
            for k in 1:3
                if dx[k] == 0
                    if !(lower[k] <= x[k] < upper[k])
                        hi = lo
                        break
                    end
                else
                    p, q = (lower[k]-x[k])/dx[k], (upper[k]-x[k])/dx[k]
                    lo, hi = max(lo, min(p,q)), min(hi, max(p,q))
                end
            end
            hi > lo || continue
            residence[i] += dt*(hi-lo)
            empty!(cuts)
            push!(cuts, lo, hi)
            for k in 1:3
                dv[k] == 0 && continue
                v0, v1 = v[k]+lo*dv[k], v[k]+hi*dv[k]
                firstedge = searchsortedfirst(edges[k], min(v0,v1))
                lastedge = searchsortedlast(edges[k], max(v0,v1))
                for n in firstedge:lastedge
                    alpha = (edges[k][n]-v[k])/dv[k]
                    lo < alpha < hi && push!(cuts, alpha)
                end
            end
            sort!(cuts)
            for n in 1:(length(cuts)-1)
                duration = dt*(cuts[n+1]-cuts[n])
                duration > 0 || continue
                vm = v + ((cuts[n]+cuts[n+1])/2)*dv
                bins = ntuple(3) do k
                    value = vm[k]
                    value < first(edges[k]) || value > last(edges[k]) ? 0 :
                        min(searchsortedlast(edges[k], value), dims[k])
                end
                if any(==(0), bins)
                    outside[i] += duration
                else
                    index = ntuple(j -> bins[kept[j]], length(kept))
                    occupancy[index...] += rates[i]*duration
                end
            end
        end
    end
    total_density = sum(rates .* residence) / side^3
    outside_density = sum(rates .* outside) / side^3
    density_in_range = sum(occupancy) / side^3
    for index in CartesianIndices(occupancy)
        volume = side^3 * prod(widths[kept[j]][index[j]] for j in eachindex(kept))
        isfinite(volume) && volume > 0 || throw(ArgumentError("Unrepresentable phase-space volume"))
        occupancy[index] /= volume
    end
    all(isfinite, occupancy) && isfinite(total_density) && isfinite(outside_density) ||
        throw(ArgumentError("PSD accumulation overflow"))
    return (; psd=occupancy, axes=map(k -> (:vx,:vy,:vz)[k], kept),
        velocity_centers_m_s=map(k -> (edges[k][1:end-1]+edges[k][2:end])/2, kept),
        velocity_edges_m_s=map(k -> edges[k], kept), all_velocity_edges_m_s=edges,
        units=length(kept)==3 ? "s^3 m^-6" : "s^2 m^-5",
        density_in_range_m3=density_in_range, density_outside_vlim_m3=outside_density,
        density_total_m3=total_density, residence_s=residence, outside_vlim_residence_s=outside,
        retcodes, species, detector_m=center, side_m=side, coordinate_system,
        interpolation=:piecewise_linear)
end

function _psd_edges(vlim, vgrid, velocity_unit)
    scale = velocity_unit == :m_s ? 1.0 : velocity_unit == :km_s ? 1000.0 :
        throw(ArgumentError("velocity_unit must be :m_s or :km_s"))
    limits = vlim isa Real ? ntuple(_ -> (-vlim,vlim), 3) :
        length(vlim)==2 && all(x -> x isa Real, vlim) ? ntuple(_ -> vlim, 3) : vlim
    counts = vgrid isa Integer ? (vgrid,vgrid,vgrid) : vgrid
    length(limits)==3 && length(counts)==3 || throw(ArgumentError("Expected three axes"))
    return ntuple(3) do k
        counts[k] isa Integer && !(counts[k] isa Bool) && counts[k] > 0 ||
            throw(ArgumentError("vgrid must specify positive integer bin counts"))
        lim = limits[k]
        length(lim)==2 || throw(ArgumentError("Each vlim must be (min,max)"))
        lo, hi = Float64(lim[1])*scale, Float64(lim[2])*scale
        isfinite(lo) && isfinite(hi) && hi > lo && isfinite(hi-lo) ||
            throw(ArgumentError("Velocity limits must be finite and increasing"))
        e = collect(range(lo,hi; length=counts[k]+1))
        all(>(0), diff(e)) || throw(ArgumentError("Velocity bins are too narrow"))
        e
    end
end

function _psd_trajectory(sol, sp)
    while !hasproperty(sol, :t)
        hasproperty(sol, :u) && length(sol.u)==1 ||
            throw(ArgumentError("Expected a trajectory or single-member ensemble"))
        sol = only(sol.u)
    end
    length(sol.t)==length(sol.u) && !isempty(sol.t) ||
        throw(ArgumentError("Trajectory needs matching nonempty times and states"))
    all(isfinite, sol.t) && all(u -> length(u)>=6 && all(isfinite,u[1:6]), sol.u) ||
        throw(ArgumentError("Trajectory has invalid times or Cartesian states"))
    if hasproperty(sol, :retcode)
        sol.retcode in (TP.ReturnCode.Success, TP.ReturnCode.Terminated) ||
            throw(ArgumentError("Failed trajectory: $(sol.retcode)"))
    end
    if hasproperty(sol, :prob) && sol.prob isa TP.TraceProblem
        p = sol.prob.p
        isapprox(p[1], sp.q/sp.m; rtol=1e-10) && isapprox(p[2], sp.m; rtol=1e-10) ||
            throw(ArgumentError("Trajectory species does not match requested species"))
    end
    return sol
end
